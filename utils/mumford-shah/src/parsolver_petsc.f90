!AUTHORS:
!Stefano Zhampini
!Hom Nath Gharti
!REFERENCE:
!PETSC documentation
! -----------------------------------------------------------------------
module parsolver_petsc
use ksp_constants
use global
use math_library_mpi,only:maxvec,minvec
use mpi_library,only:check_allocate,sync_process
use ghost_library_mpi,only:ngpart,gpart
implicit none
!------------------------------------------------------------------------------
!                    Include files
!------------------------------------------------------------------------------
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscviewer.h90"
PetscBool      flg,flg_ch,flg_lu,flg_ilu,mat_symmetry
PetscInt       petsc_solver_type
integer,parameter :: SUPERLU=2,MUMPS=3
PetscInt       ival,icntl
PetscReal      val

! Level-1 solver
Vec              xvec,bvec,local_vec,gxvec
Mat              Amat,Amatc,Fmat,Umat!,Tmat
KSP              ksp
PC               pc
PetscReal        atol,dtol,rtol
PetscInt         iter,maxiter
! For communications from local to global   
VecScatter             pscat,vscat
! Stores l2g map info 
ISLocalToGlobalMapping l2gmap                    
!PetscBool        flg

PetscInt :: nzeros_max,nzeros_min,nzerosoff_max
PetscInt,allocatable :: nnzero_diag(:),nnzero_offdiag(:)
PetscInt :: ngdof_part
PetscErrorCode   ierr
!integer :: ierr
character(len=250),private :: myfname=" => parsolver_petsc.f90"
character(len=500),private :: errsrc
integer :: mypart 
contains

subroutine petsc_initialize()
implicit none
errsrc=trim(myfname)//' => petsc_initialize'
!------------------------------------------------------------------------------
! initialize petsc
!------------------------------------------------------------------------------
call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
end subroutine petsc_initialize
!===============================================================================

subroutine petsc_create_vector()
implicit none
IS global_is,local_is

errsrc=trim(myfname)//' => petsc_create_vector'

!------------------------------------------------------------------------------
! create vector objects
!------------------------------------------------------------------------------
call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof,xvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,bvec,ierr)
CHKERRQ(ierr)

! local vector
call VecCreateSeq(PETSC_COMM_SELF,neq,local_vec,ierr)
CHKERRQ(ierr)

! objects needed for global vector scattering to local vector
! create local and global IS (index set) objects from the array of local and
! global indices
call ISCreateGeneral(PETSC_COMM_WORLD,neq,l2gdof(1:),PETSC_COPY_VALUES,global_is,ierr)
CHKERRQ(ierr)
call ISCreateStride(PETSC_COMM_SELF,neq,0,1,local_is,ierr);
CHKERRQ(ierr)
! create VecScatter object which is needed to scatter PETSc parallel vectors
call VecScatterCreate(bvec,global_is,local_vec,local_is,vscat,ierr)
CHKERRQ(ierr)
call ISDestroy(global_is,ierr) ! no longer necessary
call ISDestroy(local_is,ierr)  ! no longer necessary

end subroutine petsc_create_vector
!===============================================================================

subroutine petsc_matrix_preallocate_size()
implicit none
Vec         nnzv,nzeror_gvec,nzeror_dvec,nzeror_ovec,iproc_gvec,               &
            interface_gvec,ninterface_dvec,ninterface_ovec,nself_gvec
PetscInt :: i,istart,iend,n,n1,ncol_part,nrow_part
PetscInt :: nnzmax,lsize,idxinsert(neq),ldof(neq)
PetscInt,allocatable :: nzeros(:),ig_array(:) 
PetscScalar rval,valinsert(neq),nnzv_v(1)
PetscOffset nnzv_i
PetscInt, allocatable :: nnz(:)

PetscInt :: ig0,ig1
PetscInt :: icount,igdof,ind,maxrank0,ng,ng0,ng1,np0,ngrow
PetscInt,allocatable :: inzeror_array(:),iproc_array(:)
PetscInt,allocatable :: nnzero_diagr(:),nnzero_offdiagr(:)
PetscScalar,pointer :: nzeror_array(:),rproc_array(:)
PetscScalar,pointer :: nzeror_darray(:),nzeror_oarray(:),rnself_array(:)
PetscReal :: fac_ni,max_ni,pmax,pmin,rnid,rnioffd,rnd,rnoffd,rproc,ZERO

PetscInt :: ir,ic,igr,igc,ir0,ic0,igr0,igc0
PetscInt :: nd,noffd,nid,nioffd
PetscInt :: i_bool,i_ndof

PetscInt :: nibool,ng_interface
PetscInt,allocatable :: ibool_interface(:),ig_interface(:),isg_interface(:),   &
                        nself_array(:)
PetscInt, allocatable :: ninterface_darray(:),ninterface_oarray(:)
PetscScalar,allocatable :: rg_interface(:),rnself_lgarray(:)
PetscScalar,pointer :: rninterface_darray(:),rninterface_oarray(:)

character(len=10) :: ptail
character(len=60) :: outf_name
logical :: is_file
PetscInt :: count_diag,count_nsparse


errsrc=trim(myfname)//' => petsc_matrix_preallocate_size'

mypart= myid-1 !consistency with myrank numbering
write(ptail,'(i4)')mypart
ptail=adjustl(ptail)

!------------------------------------------------------------------------------
! predetermine or read compressed size of the sparse matrix
!------------------------------------------------------------------------------
!if(RECYCLE_STEP_VAL.gt.0 .and. exist_recycle)then
!  outf_name='recycle_petsc'//trim(ptail)
!  open(100,file=outf_name,access='stream',form='unformatted',action='read',status='old')
!  read(100)ngrow
!  allocate(nnzero_diag(ngrow),nnzero_offdiag(ngrow),stat=ierr)
!  call check_allocate(ierr,errsrc)
!  read(100)nnzero_diag
!  read(100)nnzero_offdiag
!  close(100)
!else
! count number of nonzeros per row
allocate(nzeros(neq),stat=ierr)
call check_allocate(ierr,errsrc)

nzeros=0;
do i=1,nsparse
  nzeros(krow_sparse(i))=nzeros(krow_sparse(i))+1
enddo
nzeros_max=maxvec(nzeros)
nzeros_min=minvec(nzeros)
nzerosoff_max=nzeros_max
call sync_process

! precompute ownership range OR partion layout
!
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! ng0=ng1
! all processors have equal gdofs
  ng=ng0
  ig0=mypart*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! mypart is 0-based
  if(mypart.le.maxrank0)then
    ng=ng0
    ig0=mypart*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !mypart.gt.maxrank0
    ng=ng1
    ig0=np0*ng0+(mypart-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(*,*)'ERROR: illegal value of "np0"!'
  stop
endif
call VecDuplicate(xvec,nzeror_gvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,nzeror_dvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,nzeror_ovec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,iproc_gvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,interface_gvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,nself_gvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,ninterface_dvec,ierr)
CHKERRQ(ierr)
call VecDuplicate(xvec,ninterface_ovec,ierr)
CHKERRQ(ierr)

! assign owner processor ID to each gdof (or row)
allocate(ig_array(ng),rproc_array(ng),stat=ierr)
call check_allocate(ierr,errsrc)
ig_array=(/ (i,i=ig0,ig1) /)
!rproc=real(mypart)
rproc_array=real(mypart)
call VecSetValues(iproc_gvec,ng,ig_array,rproc_array,INSERT_VALUES,ierr);
CHKERRQ(ierr)
deallocate(ig_array,rproc_array)
call VecAssemblyBegin(iproc_gvec,ierr)
CHKERRQ(ierr)
call VecAssemblyEnd(iproc_gvec,ierr)
CHKERRQ(ierr)
call VecMin(iproc_gvec,PETSC_NULL_INTEGER,pmin,ierr)
call VecMax(iproc_gvec,PETSC_NULL_INTEGER,pmax,ierr)
! copy solution to local array
allocate(iproc_array(neq),rproc_array(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(iproc_gvec,rproc_array)
iproc_array=int(rproc_array)

! assign interface ID to each gdofs
rval=1.0
! inner core
do i=1,ngpart
   nibool=gpart(i)%nnode
    allocate(ibool_interface(nibool),stat=ierr)
    call check_allocate(ierr,errsrc)
    ibool_interface=gpart(i)%node
    do i_bool=1,nibool
      do i_ndof=1,NNDOF
        igdof=ggdof(i_ndof,ibool_interface(i_bool))-1
        if(igdof.ge.0)call VecSetValues(interface_gvec,1,igdof,rval,INSERT_VALUES,ierr);
      enddo
    enddo
    deallocate(ibool_interface)
enddo
call VecAssemblyBegin(interface_gvec,ierr)
CHKERRQ(ierr)
call VecAssemblyEnd(interface_gvec,ierr)
CHKERRQ(ierr)

!call sync_process
!! stop all the MPI processes, and exit
!call MPI_FINALIZE(ierr)

! copy solution to local array
allocate(isg_interface(neq),rg_interface(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(interface_gvec,rg_interface)
isg_interface=int(rg_interface)

! estimate correction for the number of nonzero entries in the diagonal and
! nondiagonal portion
! self interface
!rval=-1.0
!call VecSet(nself_gvec,rval,ierr) ! subtract self
rval=1.0
do i=1,neq
  if(isg_interface(i).eq.1)then    
    call VecSetValues(nself_gvec,1,l2gdof(i),rval,ADD_VALUES,ierr);
  endif
enddo
call VecAssemblyBegin(nself_gvec,ierr)
CHKERRQ(ierr)
call VecAssemblyEnd(nself_gvec,ierr)
CHKERRQ(ierr)
call VecGetLocalSize(nself_gvec,n,ierr)

allocate(rnself_lgarray(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(nself_gvec,rnself_lgarray)
call VecGetArrayF90(nself_gvec,rnself_array,ierr)
allocate(nself_array(n),stat=ierr)
call check_allocate(ierr,errsrc)
nself_array=int(rnself_array(1:n))
where(nself_array.gt.0)nself_array=nself_array-1 ! subtract self
call VecRestoreArrayF90(nself_gvec,rnself_array,ierr)
call VecDestroy(nself_gvec,ierr)

!outf_name='isg_interface'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')isg_interface
!close(1)

! factor for maximum number of interfaces for each nondiagonal entry of the
! stiffness matrix
! the factor below is valid ONLY for rectagular partitioning of the global model
max_ni=8.0
fac_ni=0.0

! first element
igr0=kgrow_sparse(1)-1
igc0=kgcol_sparse(1)-1
ir0=krow_sparse(1)
ic0=kcol_sparse(1)
nd=0; noffd=0
rnid=0.; rnioffd=0.
if(iproc_array(ir0).eq.iproc_array(ic0))then
  nd=1;
  if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.1.0 .and. rnself_lgarray(ic0).gt.1.0)then
    fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
    rnid=1.0/fac_ni
  endif
else
  noffd=1
  if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.1.0 .and. rnself_lgarray(ic0).gt.1.0)then
    fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
    rnioffd=1.0/fac_ni
  endif
endif
count_diag=0
count_nsparse=nd+noffd
do i=2,nsparse
  igr=kgrow_sparse(i)-1
  igc=kgcol_sparse(i)-1
  ir=krow_sparse(i)
  ic=kcol_sparse(i)
  if(l2gdof(ir).ne.igr.or.l2gdof(ic).ne.igc)then
    print*,'strange:',l2gdof(ir),igr,l2gdof(ic),igc
    stop
  endif
  if(igr.ne.igr0)then
    ! new row starts
    ! set values computed so far
    !if(igr0==6943)count_diag=count_diag+nd
    count_nsparse=count_nsparse+nd+noffd
    rnd=real(nd)
    rnoffd=real(noffd)
    call VecSetValues(nzeror_dvec,1,igr0,rnd,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    call VecSetValues(nzeror_ovec,1,igr0,rnoffd,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    
    call VecSetValues(ninterface_dvec,1,igr0,rnid,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    call VecSetValues(ninterface_ovec,1,igr0,rnioffd,ADD_VALUES,ierr)
    CHKERRQ(ierr)

    ! reset
    nd=0; noffd=0
    rnid=0.; rnioffd=0.
    igr0=igr !kgrow_sparse(i)-1
    igc0=igc !kgcol_sparse(i)-1
    ir0=ir !krow_sparse(i)
    ic0=ic !kcol_sparse(i)

    if(iproc_array(ir0).eq.iproc_array(ic0))then
      nd=1;
      if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.0.0 .and. rnself_lgarray(ic0).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
        rnid=1.0/fac_ni
      endif
    else
      noffd=1
      if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.0.0 .and. rnself_lgarray(ic0).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
        rnioffd=1.0/fac_ni
      endif
    endif
  else
    ! count
    if(iproc_array(ir).eq.iproc_array(ic))then
      nd=nd+1;
      if(igr.ne.igc .and. rnself_lgarray(ir).gt.0.0 .and. rnself_lgarray(ic).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir),rnself_lgarray(ic)))
        rnid=rnid+(1.0/fac_ni)
      endif
    else
      noffd=noffd+1
      if(igr.ne.igc .and. rnself_lgarray(ir).gt.0.0 .and. rnself_lgarray(ic).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir),rnself_lgarray(ic)))
        rnioffd=rnioffd+(1.0/fac_ni)
      endif
    endif
  endif
  if(i.eq.nsparse)then
    ! for last
    !if(igr0==6943)count_diag=count_diag+nd
    count_nsparse=count_nsparse+nd+noffd
    rnd=real(nd)
    rnoffd=real(noffd)
    call VecSetValues(nzeror_dvec,1,igr0,rnd,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    call VecSetValues(nzeror_ovec,1,igr0,rnoffd,ADD_VALUES,ierr)
    CHKERRQ(ierr)

    call VecSetValues(ninterface_dvec,1,igr0,rnid,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    call VecSetValues(ninterface_ovec,1,igr0,rnioffd,ADD_VALUES,ierr)
    CHKERRQ(ierr)
  endif
enddo
deallocate(krow_sparse,kcol_sparse)
call sync_process
call VecAssemblyBegin(nzeror_dvec,ierr)
call VecAssemblyEnd(nzeror_dvec,ierr)
call VecAssemblyBegin(nzeror_ovec,ierr)
call VecAssemblyEnd(nzeror_ovec,ierr)

call VecAssemblyBegin(ninterface_dvec,ierr)
call VecAssemblyEnd(ninterface_dvec,ierr)
call VecAssemblyBegin(ninterface_ovec,ierr)
call VecAssemblyEnd(ninterface_ovec,ierr)

call VecGetLocalSize(nzeror_dvec,n,ierr)

call VecGetArrayF90(nzeror_dvec,nzeror_darray,ierr)
allocate(nnzero_diag(n),stat=ierr)
call check_allocate(ierr,errsrc)
nnzero_diag=int(nzeror_darray(1:n))

!outf_name='ninterface_self'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nself_array
!close(1)
!! apply correction for repeatition due to interfaces
where(nnzero_diag.gt.0)nnzero_diag=nnzero_diag-nself_array

if(myid.eq.1)print*,n,minval(nzeror_darray),maxval(nzeror_darray),minval(nnzero_diag),maxval(nnzero_diag)
call sync_process

!outf_name='nzeror_diagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nnzero_diag
!close(1)
call VecRestoreArrayF90(nzeror_dvec,nzeror_darray,ierr)
call VecDestroy(nzeror_dvec,ierr)

call VecGetArrayF90(nzeror_ovec,nzeror_oarray,ierr)
allocate(nnzero_offdiag(n))
nnzero_offdiag=int(nzeror_oarray(1:n))
!outf_name='nzeror_offdiagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nnzero_offdiag
!close(1)
call VecRestoreArrayF90(nzeror_ovec,nzeror_oarray,ierr)
call VecDestroy(nzeror_ovec,ierr)

! correction
! I do not know why but there are some DOFs where the correction exceeds by 4 or
! 8 therefore to be safe we need to subtract this from all
call VecGetArrayF90(ninterface_dvec,rninterface_darray,ierr)
!where(rninterface_darray.gt.0.0 .and. rninterface_darray.lt.1.0)rninterface_darray=1.0
allocate(ninterface_darray(n),stat=ierr)
call check_allocate(ierr,errsrc)
ninterface_darray=int(rninterface_darray(1:n))
call VecRestoreArrayF90(ninterface_dvec,rninterface_darray,ierr)
call VecDestroy(ninterface_dvec,ierr)
where(ninterface_darray.gt.0)ninterface_darray=ninterface_darray-4
where(ninterface_darray.lt.0)ninterface_darray=0
!outf_name='ninterface_diagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')ninterface_darray
!close(1)

call VecGetArrayF90(ninterface_ovec,rninterface_oarray,ierr)
!where(rninterface_oarray.gt.0.0 .and. rninterface_oarray.lt.1.0)rninterface_oarray=1.0
allocate(ninterface_oarray(n),stat=ierr)
call check_allocate(ierr,errsrc)
ninterface_oarray=int(rninterface_oarray(1:n))
call VecRestoreArrayF90(ninterface_ovec,rninterface_oarray,ierr)
call VecDestroy(ninterface_ovec,ierr)
where(ninterface_oarray.gt.0)ninterface_oarray=ninterface_oarray-8
where(ninterface_oarray.lt.0)ninterface_oarray=0
!outf_name='ninterface_offdiagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')ninterface_oarray
!close(1)

call sync_process

do i=1,nsparse
  rval=1.
  igdof=kgrow_sparse(i)-1 ! fortran index
  call VecSetValues(nzeror_gvec,1,igdof,rval,ADD_VALUES,ierr);
  CHKERRQ(ierr)
enddo
call VecAssemblyBegin(nzeror_gvec,ierr)
CHKERRQ(ierr)
call VecAssemblyEnd(nzeror_gvec,ierr)
CHKERRQ(ierr)
call VecGetLocalSize(nzeror_gvec,n,ierr)
CHKERRQ(ierr)
deallocate(kgrow_sparse,kgcol_sparse)
call VecGetArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRQ(ierr)
allocate(inzeror_array(n),stat=ierr)
call check_allocate(ierr,errsrc)
inzeror_array=int(nzeror_array(1:n))
!outf_name='nzeror'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')inzeror_array
!close(1)

call VecRestoreArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRQ(ierr)
call VecDestroy(nzeror_gvec,ierr)
CHKERRQ(ierr)
call sync_process
where(nnzero_diag.lt.0)nnzero_diag=0
where(nnzero_offdiag.lt.0)nnzero_offdiag=0
where(nnzero_diag.gt.ng)nnzero_diag=ng
call sync_process
if(myid.eq.1)print*,'success!',nzeros_max
deallocate(nzeros)
!endif

end subroutine petsc_matrix_preallocate_size
!===============================================================================

subroutine petsc_create_matrix()
implicit none
PetscInt :: istart,iend

errsrc=trim(myfname)//' => petsc_create_matrix'
!------------------------------------------------------------------------------
! create the matrix and preallocate
!------------------------------------------------------------------------------
call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
call MatSetType(Amat,MATMPIAIJ,ierr)
CHKERRQ(ierr)
call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr)
CHKERRQ(ierr)

!does not work before preallocation
!call MatGetLocalSize(Amat,nrow_part,ncol_part,ierr)

! preallocation
call MatMPIAIJSetPreallocation(Amat,PETSC_NULL_INTEGER,nnzero_diag,         &
PETSC_NULL_INTEGER,nnzero_offdiag,ierr)
CHKERRQ(ierr)

if(myid.eq.1)print*,'ngdof:',ngdof
!if(myid.eq.1)print*,'nnzero_offdiag:',nnzero_offdiag

!call MatSetOption(Amat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
!call MatSetOption(Amat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
call MatSetFromOptions(Amat,ierr)
CHKERRQ(ierr)

call MatGetOwnershipRange(Amat,istart,iend,ierr)
CHKERRQ(ierr)
call sync_process

!! convoutional stiffness matrix-------------------------------------------------
!call MatCreate(PETSC_COMM_WORLD,Amatc,ierr)
!call MatSetType(Amatc,MATMPIAIJ,ierr)
!CHKERRQ(ierr)
!call MatSetSizes(Amatc,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr)
!CHKERRQ(ierr)
!
!!does not work before preallocation
!!call MatGetLocalSize(Amat,nrow_part,ncol_part,ierr)
!
!! preallocation
!call MatMPIAIJSetPreallocation(Amatc,PETSC_NULL_INTEGER,nnzero_diagconv,         &
!PETSC_NULL_INTEGER,nnzero_offdiag,ierr)
!CHKERRQ(ierr)
!
!if(myid.eq.1)print*,'ngdof:',ngdof,minval(nnzero_diagconv),maxval(nnzero_diagconv),sum(nnzero_diagconv)
!!if(myid.eq.1)print*,'nnzero_offdiag:',nnzero_offdiag
!
!!call MatSetOption(Amat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
!!call MatSetOption(Amat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
!call MatSetFromOptions(Amatc,ierr)
!CHKERRQ(ierr)
end subroutine petsc_create_matrix
!===============================================================================

subroutine petsc_create_solver()
implicit none

!------------------------------------------------------------------------------
! Create the linear solver and set various options
!------------------------------------------------------------------------------
! define solver type
petsc_solver_type=0

! Create linear solver context

call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
CHKERRQ(ierr)
!call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr) ! PETSc version < 3.5
call KSPSetOperators(ksp,Amat,Amat,ierr) !version >= 3.5.0
CHKERRQ(ierr)
call KSPSetInitialGuessNonzero(ksp,PETSC_false,ierr)
!since the euqutions are noNDIMensionalized, the scaling is not necessary?
call KSPSetDiagonalScale(ksp,PETSC_TRUE,ierr)
CHKERRQ(ierr)
call KSPSetReusePreconditioner(ksp,PETSC_FALSE,ierr)

if(petsc_solver_type==0)then
  if(myid.eq.1)print*,'Solver type: GMRES'
  !call KSPSetType(ksp,KSPMINRES,ierr);
  !call KSPSetType(ksp,KSPBCGSL,ierr);
  call KSPSetType(ksp,KSPLGMRES,ierr);
  call KSPGetPC(ksp,pc,ierr)
  CHKERRQ(ierr)
  !call PCSetType(pc,PCNONE,ierr)
  !call PCSetType(pc,PCGAMG,ierr)
  call PCSetType(pc,PCPBJACOBI,ierr)
  !call PCSetType(pc,PCPBJACOBI,ierr)
  !call KSPSetType(ksp,KSPGCR,ierr);
  !call KSPSetType(ksp,KSPCG,ierr);
  !call KSPSetType(ksp,KSPPREONLY,ierr);
  !CHKERRQ(ierr)
  !call PCSetType(pc,PCGAMG,ierr)
  CHKERRQ(ierr)
elseif(petsc_solver_type==1)then
  if(myid.eq.1)print*,'Solver type: CG'
  call KSPSetType(ksp,KSPCG,ierr);
  !call KSPSetType(ksp,KSPPREONLY,ierr);
  CHKERRQ(ierr)
  call KSPGetPC(ksp,pc,ierr)
  CHKERRQ(ierr)
  !call PCSetType(pc,PCNONE,ierr)
  call PCSetType(pc,PCHYPRE,ierr)
  !call PCSetType(pc,PCBJACOBI,ierr)
  !call PCSetType(pc,PCGAMG,ierr)
  CHKERRQ(ierr)
  call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
  CHKERRQ(ierr)
elseif(petsc_solver_type.eq.SUPERLU)then
  if(myid.eq.1)print*,'Solver type: SUPERLU'
  flg_ilu = PETSC_FALSE;
  flg_lu     = PETSC_FALSE;
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_lu",flg_lu,flg,ierr);
  CHKERRQ(ierr)
  !PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_ilu",flg_ilu,flg,ierr);
  if(flg_lu .or. flg_ilu)then
    call KSPSetType(ksp,KSPPREONLY,ierr);
    CHKERRQ(ierr)
    call KSPGetPC(ksp,pc,ierr);
    CHKERRQ(ierr)
    if(flg_lu)then
      call PCSetType(pc,PCLU,ierr);
      CHKERRQ(ierr)
    elseif(flg_ilu)then
      call PCSetType(pc,PCILU,ierr);
      CHKERRQ(ierr)
    endif
    call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
    CHKERRQ(ierr)
    call PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU,ierr);
    CHKERRQ(ierr)
    call PCFactorSetUpMatSolverPackage(pc,ierr); ! call MatGetFactor() to create F
    CHKERRQ(ierr)
 
    call PCFactorGetMatrix(pc,Fmat,ierr);
    CHKERRQ(ierr)
    !call MatSuperluSetILUDropTol(Fmat,1.e-8,ierr);
    !CHKERRQ(ierr)
  endif
elseif(petsc_solver_type.eq.MUMPS)then
  if(myid.eq.1)print*,'Solver type: MUMPS'
  flg_lu    = PETSC_FALSE;
  flg_ch = PETSC_FALSE;
  !call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_mumps_lu",flg_lu,flg,ierr);
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_mumps_ch",flg_ch,flg,ierr);
  if(flg_lu .or. flg_ch)then
    call KSPSetType(ksp,KSPPREONLY,ierr);
    call KSPGetPC(ksp,pc,ierr);
    if(flg_lu)then
      call PCSetType(pc,PCLU,ierr);
    elseif(flg_ch)then
      call MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr); ! set MUMPS id%SYM=1
      call PCSetType(pc,PCCHOLESKY,ierr);
    endif
    call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
    CHKERRQ(ierr)
    call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,ierr);
    call PCFactorSetUpMatSolverPackage(pc,ierr); ! call MatGetFactor() to create F
  
    call PCFactorGetMatrix(pc,Fmat,ierr);
    icntl = 7; ival = 2;
    call MatMumpsSetIcntl(Fmat,icntl,ival,ierr);
    icntl = 1; val = 0.0;
    call MatMumpsSetCntl(Fmat,icntl,val,ierr);
  endif
endif

call KSPSetTolerances(ksp,KSP_RTOL,KSP_ATOL,KSP_DTOL,KSP_MAXITER,ierr)
CHKERRQ(ierr)

!  Set runtime options, e.g.,
!    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.
call KSPSetFromOptions(ksp,ierr)

end subroutine petsc_create_solver
!===============================================================================

subroutine petsc_set_ksp_operator
implicit none
errsrc=trim(myfname)//' => petsc_set_ksp_operator'                               
                                                                                 
! set ksp operators                                                              
!call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr) ! PETSc version < 3.5
call KSPSetOperators(ksp,Amat,Amat,ierr) !version >= 3.5.0                       
CHKERRQ(ierr)
end subroutine petsc_set_ksp_operator
!===============================================================================

subroutine petsc_set_ksp_operatorconv
implicit none
errsrc=trim(myfname)//' => petsc_set_ksp_operator'                               
                                                                                 
! set ksp operators                                                              
!call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr) ! PETSc version < 3.5
call KSPSetOperators(ksp,Amatc,Amatc,ierr) !version >= 3.5.0                       
CHKERRQ(ierr)
end subroutine petsc_set_ksp_operatorconv
!===============================================================================

subroutine petsc_set_stiffness_matrix(storekmat)
use math_library_mpi,only:sumscal
use ieee_arithmetic
implicit none
real(kind=kreal),intent(in) :: storekmat(:,:,:)
integer :: i,i_elmt,ielmt,j,n
integer :: ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)

PetscInt irow,jcol,istart,iend,ndiag,noffdiag,ng,ngrow
integer :: ncols,ncols_vali,nsparse1
PetscInt,allocatable :: cols(:),nnzero_diag(:),nnzero_offdiag(:)
real(kind=8),allocatable :: vals(:)
character(len=10) :: ptail
character(len=60) :: outf_name
Vec   vdiag
PetscScalar rval
PetscScalar,pointer :: diag_array(:)
PetscInt,allocatable :: ibarray(:)
PetscScalar,pointer :: barray(:)
PetscReal :: stol,mnorm
logical :: isu,isp,isphi
logical,allocatable :: iselmt(:)

real(kind=8) :: kmat(NEDOF,NEDOF),xval

interface
  subroutine issymmetric(mypart,mat)
  implicit none
  integer,parameter :: kreal=selected_real_kind(15)
  integer,intent(in) :: mypart
  real(kind=kreal),intent(in) :: mat(:,:)
  end subroutine issymmetric
end interface
! Set and assemble matrix.
!  - Note that MatSetValues() uses 0-based row and column numbers
!  in Fortran as well as in C (as set here in the array "col").
! stage 0: store all elements

mypart= myid-1 !consistency with myrank numbering
call MatZeroEntries(Amat,ierr)
CHKERRQ(ierr)
call sync_process
rval=1.0

! entirely in solid
do i_elmt=1,nelmt
  ielmt=i_elmt
  ggdof_elmt=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))
  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  do i=1,NEDOF
    do j=1,NEDOF
    irow=i; jcol=j
    if(ggdof_elmt(irow).ge.0.and.ggdof_elmt(jcol).ge.0)then !.and.storekmat_intact_ic(i,j,i_elmt).ne.0.0_kreal)then
      xval=storekmat(i,j,ielmt)
      if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then
        write(*,*)'ERROR: stiffness matrix has nonfinite value/s!',mypart,ielmt,&
        mat_id(ielmt),xval,minval(abs(storekmat)),maxval(abs(storekmat))
        stop
      endif
      call MatSetValues(Amat,1,ggdof_elmt(irow),1,ggdof_elmt(jcol),           &
      storekmat(i,j,ielmt),ADD_VALUES,ierr)
      CHKERRQ(ierr)
    endif
    enddo
  enddo
enddo

call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)

call MatCreateVecs(Amat,vdiag,PETSC_NULL_OBJECT,ierr)
call MatGetDiagonal(Amat,vdiag,ierr)
call VecGetLocalSize(vdiag,n,ierr)
CHKERRQ(ierr)
call VecGetArrayF90(vdiag,diag_array,ierr)
CHKERRQ(ierr)
print*,'NZEROs in diagonal:',mypart,n,count(diag_array==0.),minval(abs(diag_array)),maxval(abs(diag_array))
call VecRestoreArrayF90(vdiag,diag_array,ierr)
call sync_process
call VecDestroy(vdiag,ierr)
!endif
end subroutine petsc_set_stiffness_matrix
!=======================================================

subroutine petsc_set_stiffness_matrixconv(nup,nvalency,interpfgll,gdof_elmt,nodalnu0,nodalgauss)
use math_constants,only:gtol
use math_library,only:upind_i2j,upinds_i2js
use math_library_mpi,only:sumscal
use global,only:storederiv,storejw
use ieee_arithmetic
implicit none
integer,intent(in) :: nup
integer,intent(in) :: nvalency(:)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
integer,intent(in) :: gdof_elmt(:,:)
real(kind=kreal),intent(in) :: nodalnu0(:,:),nodalgauss(0:)

integer :: i,i_elmt,ielmt,inum,istat,iup,j,j_node,n,nnzero
integer :: egdof(nedof),ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF),num(ngll)
integer :: i_node,inode,iups(nnode),jnodes(nnode),node2gll(nnode)
integer :: ncols,ncols_vali,nsparse1
real(kind=kreal) :: elmtjw(ngll),nugll(nedof)
real(kind=kreal) :: dinterpf(ndim,ngll),interpf(ngll,1)
real(kind=kreal) :: kmat(nedof,nedof),kmat_gradm(nedof,nedof),kmat_m(nedof,nedof)
PetscInt irow,jcol,istart,iend,ndiag,noffdiag,ng,ngrow
PetscInt,allocatable :: cols(:),nnzero_diag(:),nnzero_offdiag(:)
PetscInt,allocatable :: icols(:)
real(kind=8),allocatable :: vals(:)
real(kind=8),allocatable :: gauss(:)
PetscScalar,allocatable :: gkmat(:,:),kvec(:)
character(len=10) :: ptail
character(len=60) :: outf_name
Vec   vdiag
PetscScalar rval
PetscScalar,pointer :: diag_array(:)
PetscInt,allocatable :: ibarray(:)
PetscScalar,pointer :: barray(:)
PetscReal :: stol,mnorm,zero
logical :: isu,isp,isphi
logical,allocatable :: iselmt(:)

real(kind=8) :: xval

! Set and assemble matrix.
!  - Note that MatSetValues() uses 0-based row and column numbers
!  in both Fortran/C.
! stage 0: store all elements

mypart= myid-1 !consistency with myrank numbering
call MatZeroEntries(Amatc,ierr)
CHKERRQ(ierr)
call sync_process
rval=1.0
zero=0.0
allocate(gkmat(nedof,nnode),stat=istat)
if(istat.ne.0)then
  write(*,*)'ERROR: not enough memory for "gkmat"!'
  stop
endif
! for image problem all elements are identical
elmtjw=storejw(:,1)
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  node2gll(num)=(/ (i,i=1,ngll) /)
enddo

allocate(gauss(nnode))
jnodes=(/ (i_node,i_node=1,nnode) /)
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,i_elmt)
  ggdof_elmt=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))
  if(any(num.ne.ggdof_elmt))stop 'strange!'
  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0

  ! compute stiffness

  ! constant part of the convolutional term
  gkmat=zero
  do i=1,ngll
    inode=num(i)
    do j_node=1,nnode
      iup=upind_i2j(nup,nnode,inode,j_node)
      if(nodalgauss(iup).le.gtol)cycle
      gkmat(i,j_node)=nodalgauss(iup)
    enddo
    !iups=upinds_i2js(nup,nnode,inode,nnode,jnodes)
    !gauss=nodalgauss(iups)
    !where(gauss.le.gtol)gauss=zero ! to make martix sparse
    !gkmat(i,:)=gauss
    gkmat(i,:)=elmtjw(node2gll(jnodes))*(gkmat(i,:)*nvalency)
  enddo
  egdof=gdof_elmt(:,i_elmt) 
  nugll=nodalnu0(1,num)
  kmat=zero; kmat_gradm=zero; kmat_m=zero
  do i=1,ngll
    iups=upinds_i2js(nup,nnode,num(i),ngll,num)
    interpf(:,1)=interpfgll(i,:)
    dinterpf=storederiv(:,:,i,i_elmt)
     
    kmat_gradm = kmat_gradm+gam_img*nugll(i)*nugll(i)* &
                 matmul(transpose(dinterpf),dinterpf)*storejw(i,i_elmt)
    kmat_m = kmat_m+matmul(interpf,transpose(interpf))*storejw(i,i_elmt)
    !kmat = kmat + (kmat_gradm+kmat_m)*storejw(i,i_elmt)
  end do ! i=1,ngll
  gkmat=matmul(kmat_m,gkmat)
  gkmat(:,egdof)=gkmat(:,egdof)+kmat_gradm

  do i=1,ngll
    ! count non zeros
    nnzero=count(abs(gkmat(i,:)).gt.zero)
    !print*,'hello0:',nnzero
    allocate(icols(nnzero),kvec(nnzero))
    kvec=zero
    icols=-1
    inum=0
    do i_node=1,nnode
      if(abs(gkmat(i,i_node)).gt.zero)then
        inum=inum+1
        icols(inum)=i_node
        kvec(inum)=gkmat(i,i_node)
      endif
      !call MatSetValues(Amatc,ngll,num,nnzero,icols,kvec,ADD_VALUES,ierr)
    enddo
    ! petsc uses 0 based indices
    irow=num(i)-1
    icols=icols-1
    call MatSetValues(Amatc,1,irow,nnzero,icols,kvec,ADD_VALUES,ierr)
    CHKERRQ(ierr)
    deallocate(icols,kvec)
  enddo
enddo
print*,'success setting up the matrix!'

call MatAssemblyBegin(Amatc,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(Amatc,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)

call MatCreateVecs(Amatc,vdiag,PETSC_NULL_OBJECT,ierr)
call MatGetDiagonal(Amatc,vdiag,ierr)
call VecGetLocalSize(vdiag,n,ierr)
CHKERRQ(ierr)
call VecGetArrayF90(vdiag,diag_array,ierr)
CHKERRQ(ierr)
print*,'NZEROs in diagonal:',mypart,n,count(diag_array==0.),minval(abs(diag_array)),maxval(abs(diag_array))
call VecRestoreArrayF90(vdiag,diag_array,ierr)
call sync_process
call VecDestroy(vdiag,ierr)
!endif
end subroutine petsc_set_stiffness_matrixconv
!=======================================================

subroutine petsc_set_vector(rload)
!use global,only:l2gdof,nelmt,NEDOF
use ieee_arithmetic
implicit none
PetscScalar,intent(in) :: rload(0:)
PetscScalar fval,ZERO !,none,ONE
PetscInt    i,gdof
PetscInt    istart,iend
PetscScalar,pointer :: larray(:)
integer :: i_elmt,ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)
real(kind=kreal) :: xval

call VecGetOwnershipRange(bvec,istart,iend,ierr)
CHKERRQ(ierr)

mypart= myid-1 !consistency with myrank numbering

ZERO=0.0
call VecSet(bvec,ZERO,ierr)
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)
if(myid.eq.1)print*,'vector setting & assembly complete!'
end subroutine petsc_set_vector
!=======================================================

subroutine petsc_set_initialguess(rload)
!use global,only:l2gdof,NEDOF
implicit none
PetscScalar,intent(in) :: rload(0:)
PetscScalar fval,ZERO !,none,ONE
PetscInt    i,gdof
PetscInt    istart,iend
PetscScalar,pointer :: larray(:)
integer :: i_elmt,ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)
real(kind=kreal) :: xval

call VecGetOwnershipRange(bvec,istart,iend,ierr)
CHKERRQ(ierr)

mypart= myid-1 !consistency with myrank numbering
ZERO=0.0
call VecSet(bvec,ZERO,ierr)
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)
if(myid.eq.1)print*,'vector setting & assembly complete!'

end subroutine petsc_set_initialguess
!=======================================================

subroutine petsc_solve(sdata,cg_iter)
implicit none
MatNullSpace nullspace
PetscScalar sdata(:)
PetscInt    cg_iter
PetscInt    ireason
PetscScalar,pointer :: larray(:)

mypart= myid-1 !consistency with myrank numbering

!! null space
!call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!!call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_OBJECT, nullspace,ierr);
!call MatSetNullSpace(Amat,nullspace,ierr)
!call MatNullSpaceRemove(nullspace,bvec,PETSC_NULL_OBJECT,ierr)
!call MatNullSpaceDestroy(nullspace,ierr)
!TMP 
!TMP !call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_OBJECT, nullspace,ierr);
!TMP !call KSPSetNullSpace(ksp, nullspace,ierr);
!TMP !call MatNullSpaceDestroy(nullspace,ierr);

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call KSPSolve(ksp,bvec,xvec,ierr)

! View solver info; we could instead use the option -ksp_view
call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Check the error
!call VecAXPY(x,none,u,ierr)
!call VecNorm(x,NORM_2,norm,ierr)
call KSPGetConvergedReason(ksp,ireason,ierr)
call KSPGetIterationNumber(ksp,cg_iter,ierr)
if (mypart.lt.1)then
  print*,'converged reason',ireason
  print*,'Iterations:',cg_iter
endif

! copy solution to local array
call scatter_globalvec(xvec,sdata)
end subroutine petsc_solve
!=======================================================

subroutine scatter_globalvec(global_vec,larray)
implicit none

Vec         global_vec
PetscScalar larray(:)
PetscScalar a_v(1)
PetscOffset a_i
PetscInt    n

PetscScalar,pointer :: array_data(:)
!call VecScatterBegin(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_REVERSE,ierr)
!CHKERRQ(ierr)
!call VecScatterEnd(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_REVERSE,ierr)
!CHKERRQ(ierr)
!call VecGetSize(local_vec,n,ierr)
!CHKERRQ(ierr)
!call VecGetArray(local_vec,a_v,a_i,ierr)
!CHKERRQ(ierr)
!larray(1:n)=a_v(a_i+1:a_i+n) ! TODO: use FORTRAN POINTER
!if(myid.eq.1)print*,'larray:',minval(larray(1:n)),maxval(larray(1:n))
!call VecRestoreArray(local_vec,a_v,a_i,ierr)
!CHKERRQ(ierr)

call VecScatterBegin(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
CHKERRQ(ierr)
call VecScatterEnd(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
CHKERRQ(ierr)
call VecGetSize(local_vec,n,ierr)
call VecGetArrayF90(local_vec,array_data,ierr)
CHKERRQ(ierr)
larray(1:n)=array_data(1:n)
call VecRestoreArrayF90(local_vec,array_data,ierr)
CHKERRQ(ierr)
return
end subroutine scatter_globalvec
!=======================================================

subroutine petsc_load()
PetscViewer viewer

call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"PetscMatVecStore",FILE_MODE_READ,viewer,ierr)
call MatLoad(Amat,viewer,ierr)
call VecLoad(bvec,viewer,ierr)
call PetscViewerDestroy(viewer,ierr)
end subroutine petsc_load
!=======================================================

subroutine petsc_save()
PetscViewer viewer

call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"PetscMatVecStore",FILE_MODE_WRITE,viewer,ierr)
call MatView(Amat,viewer,ierr)
call VecView(bvec,viewer,ierr)
call PetscViewerDestroy(viewer,ierr)
end subroutine petsc_save
!=======================================================

subroutine petsc_finalize()
implicit none

! Free work space.  All PETSc objects should be destroyed when they
! are no longer needed.

call VecDestroy(xvec,ierr)
call VecDestroy(bvec,ierr)
call MatDestroy(Amat,ierr)
call KSPDestroy(ksp,ierr)
call VecScatterDestroy(vscat,ierr)
call PetscFinalize(ierr)

end subroutine petsc_finalize
!=======================================================

end module parsolver_petsc
!=======================================================
