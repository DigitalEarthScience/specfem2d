!AUTHORS:
!Stefano Zhampini
!Hom Nath Gharti
!REFERENCE:
!PETSC documentation
! -----------------------------------------------------------------------
module solver_petsc
use ksp_constants
use global
use math_library_mpi,only:maxvec,minvec
use mpi_library,only:check_allocate,sync_process
use ghost_library_mpi,only:ngpart,gpart
implicit none
!------------------------------------------------------------------------------
!                    Include files
!------------------------------------------------------------------------------
!
! This program uses CPP for preprocessing, as indicated by the use of
! PETSc include files in the directory petsc/include/finclude.  This
! convention enables use of the CPP preprocessor, which allows the use
! of the #include statements that define PETSc objects and variables.
!
! Use of the conventional Fortran include statements is also supported
! In this case, the PETsc include files are located in the directory
! petsc/include/foldinclude.
!
! Since ONE must be very careful to include each file no more than once
! in a Fortran routine, application programmers must exlicitly list
! each file needed for the various PETSc compONEnts within their
! program (unlike the C/C++ interface).
!
! See the Fortran section of the PETSc users manual for details.
!
! The following include statements are required for KSP Fortran programs:
!   petscsys.h       - base PETSc routines
!   petscvec.h    - vectors
!   petscmat.h    - matrices
!   petscksp.h    - Krylov subspace methods
!   petscpc.h     - preconditiONErs
! Other include statements may be needed if using additional PETSc
! routines in a Fortran program, e.g.,
!   petscviewer.h - viewers
!   petscis.h     - index sets
!
!   include "petscsys.h"
!   include "petscvec.h"
!   include "petscmat.h"
!   include "petscksp.h"
!   include "petscpc.h"
!#include <finclude/petsckspdef.h>
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscviewer.h90"
PetscBool      flg,flg_ch,flg_lu,flg_ilu,mat_symmetry
PetscInt       solver_type
integer,parameter :: SUPERLU=2,MUMPS=3
PetscInt       ival,icntl
PetscReal      val

! Level-1 solver
Vec              xvec,bvec,local_vec,gxvec
Mat              Amat,Tmat,Fmat,Umat
KSP              ksp
PC               pc
PetscReal        atol,dtol,rtol
PetscInt         iter,maxiter
! For communications from local to global   
VecScatter             pscat,vscat,vscat_all
! Stores l2g map info 
ISLocalToGlobalMapping l2gmap                    
!PetscBool        flg

!! Level-2 solver
!Vec              x,b,u
!Mat              A
!KSP              ksp
!PC               pc
!PetscReal        norm,atol,dtol,rtol
!PetscInt         iter,maxiter
PetscInt :: nzeros_max,nzeros_min,nzerosoff_max
PetscInt :: ngdof_part
PetscInt :: ig0,ig1
PetscErrorCode   ierr
!integer :: ierr
character(len=250),private :: myfname=" => solver_petsc.f90"
character(len=500),private :: errsrc
contains
!=======================================================
! Level-1 solver
!=======================================================
subroutine petsc_initialize()
implicit none
Vec         nnzv,nzeror_gvec,nzeror_dvec,nzeror_ovec,iproc_gvec,               &
            interface_gvec,ninterface_dvec,ninterface_ovec,nself_gvec
PetscInt :: i,istart,iend,n,n1,ncol_part,nrow_part
PetscInt :: nnzmax,lsize,idxinsert(neq),ldof(neq)
PetscInt,allocatable :: nzeros(:),ig_array(:) 
PetscScalar rval,valinsert(neq),nnzv_v(1)
PetscOffset nnzv_i
PetscInt, allocatable :: nnz(:)
IS global_is,local_is

PetscInt :: icount,igdof,ind,maxrank0,ng,ng0,ng1,np0,ngrow
PetscInt,allocatable :: inzeror_array(:),iproc_array(:),nzeros_row(:)
PetscInt,allocatable :: nnzero_diag(:),nnzero_offdiag(:)
PetscInt,allocatable :: nnzero_diagr(:),nnzero_offdiagr(:)
PetscScalar,pointer :: nzeror_array(:),rproc_array(:)
PetscScalar,pointer :: nzeror_darray(:),nzeror_oarray(:),rnself_array(:)
PetscReal :: fac_ni,max_ni,pmax,pmin,rnid,rnioffd,rnd,rnoffd,rproc,ZERO

PetscInt :: ir,ic,igr,igc,ir0,ic0,igr0,igc0
PetscInt :: nd,noffd,nid,nioffd
PetscInt :: i_bool,i_ndof,ncount

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
errsrc=trim(myfname)//' => petsc_initialize'
!------------------------------------------------------------------------------
! initialize petsc
!------------------------------------------------------------------------------
call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',ngdof,flg,ierr)
!if(myrank==0)print*,'hi0!'
!call sync_process
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute the matrix and right-hand-side vector that define
! the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Create matrix. When using MatCreate(), the matrix format can
! be specified at runtime.

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

call VecScatterCreateToAll(xvec,vscat_all,gxvec,ierr)

write(ptail,'(i4)')myrank
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
ncount=0
do i=1,nsparse
  nzeros(krow_sparse(i))=nzeros(krow_sparse(i))+1
  if(kgrow_sparse(i)==6944)ncount=ncount+1
enddo
nzeros_max=maxvec(nzeros)
nzeros_min=minvec(nzeros)
nzerosoff_max=nzeros_max
!nzeros_max=4*nzeros_max
!nzeros=nzeros
!nzeros=5*nzeros
!print*,'nzeros 6943rd row:',myrank,ncount
call sync_process
if(myrank==0)print*,'nzeros in 1th index:',nzeros(1)
if(myrank==0)print*,'ngdof:',ngdof,' nzeros_max:',nzeros_max,' nzeros_min:',nzeros_min,count(krow_sparse==1)

! precompute ownership range OR partion layout
!
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! ng0=ng1
! all processors have equal gdofs
  ng=ng0
  ig0=myrank*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! myrank is 0-based
  if(myrank.le.maxrank0)then
    ng=ng0
    ig0=myrank*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !myrank.gt.maxrank0
    ng=ng1
    ig0=np0*ng0+(myrank-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(*,*)'ERROR: illegal value of "np0"!'
  stop
endif
if(myrank==0)print*,'OK0:',ng0,ng1,ng,ig0,ig1
call sync_process
allocate(nzeros_row(ng),stat=ierr)
call check_allocate(ierr,errsrc)
nzeros_row=0
if(myrank==0)then
  open(1,file='test_file_proc',action='write',status='replace')
  write(1,*)ng,ig0,ig1
endif
do i=1,nsparse
 if(kgrow_sparse(i)-1.ge.ig0 .and. kgrow_sparse(i)-1.le.ig1)then
   ind=kgrow_sparse(i)-ig0 ! fortran indexing
   if(myrank==0)write(1,*)ind,kgrow_sparse(i)
   !if(myrank==0)print*,'ind:',ind,ng,kgrow_sparse(i),ig0,ig1,maxval(l2gdof) 
   nzeros_row(ind)=nzeros_row(ind)+1
 endif
enddo
if(myrank==0)close(1)
!nzeros_row=2*nzeros_row
if(myrank==0)print*,'OK1:',nzeros_row(1),minval(nzeros_row),maxval(nzeros_row)
call sync_process
outf_name='precomp_nonzeros'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')nzeros_row
close(1)
deallocate(nzeros_row)
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
if(myrank==0)print*,'vector'

! assign owner processor ID to each gdof (or row)
allocate(ig_array(ng),rproc_array(ng),stat=ierr)
call check_allocate(ierr,errsrc)
ig_array=(/ (i,i=ig0,ig1) /)
!rproc=real(myrank)
rproc_array=real(myrank)
call VecSetValues(iproc_gvec,ng,ig_array,rproc_array,INSERT_VALUES,ierr);
CHKERRQ(ierr)
deallocate(ig_array,rproc_array)
call VecAssemblyBegin(iproc_gvec,ierr)
CHKERRQ(ierr)
call VecAssemblyEnd(iproc_gvec,ierr)
CHKERRQ(ierr)
call VecMin(iproc_gvec,PETSC_NULL_INTEGER,pmin,ierr)
call VecMax(iproc_gvec,PETSC_NULL_INTEGER,pmax,ierr)
if(myrank==0)print*,'iproc range',pmin,pmax; call sync_process
!call VecGetArrayF90(iproc_gvec,rproc_array,ierr)
!CHKERRQ(ierr)
!allocate(iproc_array(ng))
!iproc_array=int(rproc_array(1:n))
!call VecRestoreArrayF90(iproc_gvec,rproc_array,ierr)
!CHKERRQ(ierr)
! copy solution to local array
allocate(iproc_array(neq),rproc_array(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(iproc_gvec,rproc_array)
iproc_array=int(rproc_array)
if(myrank==0)print*,'vector3 iproc',minval(iproc_array),maxval(iproc_array); call sync_process
!!TODO: use local scatter
!call VecScatterCreateToAll(iproc_gvec,pscat,iproc_garray,ierr);
!call VecScatterBegin(pscat,iproc_gvec,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
!call VecScatterEnd(pscat,iproc_gvec,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
!call VecScatterDestroy(pscat);

! assign interface ID to each gdofs
rval=1.0
! inner core
do i=1,ngpart
   nibool=gpart(i)%nnode
    allocate(ibool_interface(nibool),stat=ierr)
    call check_allocate(ierr,errsrc)
    ibool_interface=gpart(i)%node
    !ng_interface=nibool*NNDOF
    !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
    !ig_interface=reshape(ggdof_ic(1,ibool_interface), (/ ng_interface /) )
    !ig_interface=ig_interface-1
    !rg_interface=1.0
    !call VecSetValues(interface_gvec,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
    !deallocate(ibool_interface,ig_interface,rg_interface)
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

if(myrank==0)print*,'OH0'
call sync_process
allocate(rnself_lgarray(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(nself_gvec,rnself_lgarray)
if(myrank==0)print*,'OH1'
call sync_process
call VecGetArrayF90(nself_gvec,rnself_array,ierr)
if(myrank==0)print*,'OH2'
call sync_process
allocate(nself_array(n),stat=ierr)
call check_allocate(ierr,errsrc)
nself_array=int(rnself_array(1:n))
where(nself_array.gt.0)nself_array=nself_array-1 ! subtract self
call VecRestoreArrayF90(nself_gvec,rnself_array,ierr)
call VecDestroy(nself_gvec,ierr)

if(myrank==0)print*,'maximum value of nself:',maxval(nself_array)
!call sync_process
!! stop all the MPI processes, and exit
!call MPI_FINALIZE(ierr)
!stop
!! count nonzero entries in the diagonal and nondiagonal portion
!ZERO=0.
!call VecSet(nzeror_dvec,ZERO,ierr)
!call VecSet(nzeror_ovec,ZERO,ierr)

outf_name='isg_interface'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')isg_interface
close(1)

if(myrank==0)open(11,file='test_interface',action='write',status='replace')
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
    if(myrank==0)write(11,*)ir,ic,isg_interface(ir),isg_interface(ic)
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
call sync_process
print*,'counted nsparse:',myrank,nsparse,count_nsparse,count_nsparse-nsparse
print*,'counted ndiag:',myrank,count_diag
call sync_process
deallocate(krow_sparse,kcol_sparse)
if(myrank==0)close(11)
call sync_process
if(myrank==0)print*,'vector4:',NNDOF; call sync_process
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

outf_name='ninterface_self'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')nself_array
close(1)
!! apply correction for repeatition due to interfaces
where(nnzero_diag.gt.0)nnzero_diag=nnzero_diag-nself_array

if(myrank==0)print*,n,minval(nzeror_darray),maxval(nzeror_darray),minval(nnzero_diag),maxval(nnzero_diag)
call sync_process

outf_name='nzeror_diagonal'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')nnzero_diag
close(1)
call VecRestoreArrayF90(nzeror_dvec,nzeror_darray,ierr)
call VecDestroy(nzeror_dvec,ierr)

call VecGetArrayF90(nzeror_ovec,nzeror_oarray,ierr)
allocate(nnzero_offdiag(n))
nnzero_offdiag=int(nzeror_oarray(1:n))
outf_name='nzeror_offdiagonal'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')nnzero_offdiag
close(1)
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
outf_name='ninterface_diagonal'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')ninterface_darray
close(1)

call VecGetArrayF90(ninterface_ovec,rninterface_oarray,ierr)
!where(rninterface_oarray.gt.0.0 .and. rninterface_oarray.lt.1.0)rninterface_oarray=1.0
allocate(ninterface_oarray(n),stat=ierr)
call check_allocate(ierr,errsrc)
ninterface_oarray=int(rninterface_oarray(1:n))
call VecRestoreArrayF90(ninterface_ovec,rninterface_oarray,ierr)
call VecDestroy(ninterface_ovec,ierr)
where(ninterface_oarray.gt.0)ninterface_oarray=ninterface_oarray-8
where(ninterface_oarray.lt.0)ninterface_oarray=0
outf_name='ninterface_offdiagonal'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')ninterface_oarray
close(1)

!nnzero_diag=nnzero_diag-ninterface_darray
!nnzero_offdiag=nnzero_offdiag-ninterface_oarray
if(myrank.eq.0)print*,'corrections:',minval(ninterface_darray),                &
maxval(ninterface_darray),minval(ninterface_oarray),                           &
maxval(ninterface_oarray),minval(nnzero_diag),maxval(nnzero_diag),       &
minval(nnzero_offdiag),maxval(nnzero_offdiag),count(nnzero_offdiag.eq.-1)
call sync_process

!where(nnzero_offdiag==0)nnzero_offdiag=1
!call PetscSplitOwnership(PETSC_COMM_WORLD,ngdof_part,ngdof,ierr)
!print*,'partition:',myrank,ngdof,ngdof_part
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
if(myrank==0)print*,'size of vector:',ng,n,minval(kgrow_sparse),ig0
deallocate(kgrow_sparse,kgcol_sparse)
call VecGetArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRQ(ierr)
allocate(inzeror_array(n),stat=ierr)
call check_allocate(ierr,errsrc)
inzeror_array=int(nzeror_array(1:n))
outf_name='nzeror'//trim(ptail)
open(1,file=outf_name,action='write',status='replace')
write(1,'(i4)')inzeror_array
close(1)

call VecRestoreArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRQ(ierr)
call VecDestroy(nzeror_gvec,ierr)
CHKERRQ(ierr)


!!!allocate(nnzero_diagr(ng),nnzero_offdiagr(ng))
!!!outf_name='nonzeros'//trim(ptail)
!!!inquire(file=outf_name,exist=is_file)
!!!if(is_file)then
!!!  open(1,file=outf_name,action='read',status='old')
!!!  do i=1,ng
!!!    read(1,*)nnzero_diagr(i),nnzero_offdiagr(i),inzeror_array(i)
!!!  enddo
!!!  close(1)
!!!
!!!  if(any(nnzero_diagr.gt.nnzero_diag))then
!!!    print*,'ohhhh diagonal:',myrank,maxval(nnzero_diagr-nnzero_diag)
!!!    stop
!!!  endif 
!!!  if(any(nnzero_offdiagr.gt.nnzero_offdiag))then
!!!    print*,'ohhhh offdiagonal:',myrank,maxval(nnzero_offdiagr-nnzero_offdiag)
!!!    stop
!!!  endif 
!!!  nnzero_diag=nnzero_diagr
!!!  nnzero_offdiag=nnzero_offdiagr
!!!  if(myrank==0)print*,'WARNING: matrix size is read from file!'
!!!endif
if(myrank==0)print*,'WHAT!'
call sync_process
where(nnzero_diag.lt.0)nnzero_diag=0
where(nnzero_offdiag.lt.0)nnzero_offdiag=0
where(nnzero_diag.gt.ng)nnzero_diag=ng
if(myrank==0)print*,'OH0!',nzeros_max
call sync_process
print*,'preallocation size range:',myrank,minval(nnzero_diag),              &
maxval(nnzero_diag),minval(nnzero_offdiag),maxval(nnzero_offdiag)
if(myrank==0)print*,'preallocation size ends!',nzeros_max
call sync_process
if(myrank==0)print*,'success!',nzeros_max
deallocate(nzeros)
!endif

!------------------------------------------------------------------------------
! create the matrix and preallocate
!------------------------------------------------------------------------------
call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
call MatSetType(Amat,MATMPIAIJ,ierr)
CHKERRQ(ierr)
call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr)
CHKERRQ(ierr)
!call MatSetOption(Amat,MAT_SYMMETRIC,PETSC_TRUE,ierr);
!CHKERRQ(ierr)

!does not work before preallocation
!call MatGetLocalSize(Amat,nrow_part,ncol_part,ierr)
!print*,'local size:',myrank,nrow_part,ncol_part

! preallocation
!nzeros_max=5*nzeros_max
!call MatMPIAIJSetPreallocation(Amat,nzeros_max,PETSC_NULL_INTEGER,nzeros_max,     &
!PETSC_NULL_INTEGER,ierr)
!CHKERRQ(ierr)

!call MatMPIAIJSetPreallocation(Amat,nzeros_max,inzeror_array,nzerosoff_max,   &
!PETSC_NULL_INTEGER,ierr)
!CHKERRQ(ierr)
!nnzero_diag=5*nnzero_diag
!nnzero_offdiag=5*nnzero_offdiag
print*,'sizes:',ng,n,size(nnzero_diag),size(nnzero_offdiag)
call MatMPIAIJSetPreallocation(Amat,PETSC_NULL_INTEGER,nnzero_diag,         &
PETSC_NULL_INTEGER,nnzero_offdiag,ierr)
CHKERRQ(ierr)

!call MatMPIAIJSetPreallocation(Amat,nzeros_max,inzeror_array,nzeros_max,     &
!inzeror_array,ierr)
!CHKERRQ(ierr)

!call MatMPIAIJSetPreallocation(Amat,nzeros_max,nzeros,nzeros_max,     &
!PETSC_NULL_INTEGER,ierr)
!call MatMPIAIJSetPreallocation(Amat,PETSC_NULL_INTEGER,nzeros,     &
!PETSC_NULL_INTEGER,nzeros,ierr)
!call MatMPIAIJSetPreallocation(Amat,nnzmax,PETSC_NULL_INTEGER,nnzmax,     &
!PETSC_NULL_INTEGER,ierr)
!CHKERRQ(ierr)
!call MatSeqAIJSetPreallocation(Amat,nzeros_max,nzeros,ierr)
!CHKERRQ(ierr)

!call MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,   &
!nzeros_max,PETSC_NULL_INTEGER,nzeros_max,PETSC_NULL_INTEGER,A,ierr)

if(myrank==0)print*,'ngdof:',ngdof
!if(myrank==0)print*,'Matrix size:',size(Amat,1),size(Amat,2)
!call MatSetOption(Amat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
!call MatSetOption(Amat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
call MatSetFromOptions(Amat,ierr)
CHKERRQ(ierr)
!call MatSetUp(Amat,ierr)
if(myrank==0)print*,'ierr1:',ierr


call MatGetOwnershipRange(Amat,istart,iend,ierr)
CHKERRQ(ierr)
call sync_process
print*,'ownership range:',myrank,istart,ig0,iend-1,ig1
!!call MatZeroEntries(Amat,ierr)
if(myrank==0)print*,'ierr4:',ierr
CHKERRQ(ierr)
call sync_process
if(myrank==0)print*,'matrix'
!call sync_process

!! Null space matrix for rigid body motion
!call MatCreate(PETSC_COMM_WORLD,Umat,ierr)
!call MatSetType(Umat,MATMPIAIJ,ierr)
!CHKERRQ(ierr)
!call MatSetSizes(Umat,PETSC_DECIDE,PETSC_DECIDE,ngdof,6,ierr)
!CHKERRQ(ierr)
!call MatMPIAIJSetPreallocation(Umat,PETSC_NULL_INTEGER,3,         &
!PETSC_NULL_INTEGER,3,ierr)
!CHKERRQ(ierr)
!call MatSetFromOptions(Umat,ierr)
!CHKERRQ(ierr)
if(myrank==0)print*,'ierr1:',ierr
!------------------------------------------------------------------------------
! Create the linear solver and set various options
!------------------------------------------------------------------------------
! define solver type
solver_type=0

! Create linear solver context

call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
CHKERRQ(ierr)
!call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr) ! PETSc version < 3.5
call KSPSetOperators(ksp,Amat,Amat,ierr) !version >= 3.5.0
CHKERRQ(ierr)
call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
!since the euqutions are noNDIMensionalized, the scaling is not necessary?
!call KSPSetDiagonalScale(ksp,PETSC_TRUE,ierr)
!CHKERRQ(ierr)

if(solver_type==0)then
  if(myrank==0)print*,'Solver type: GMRES'
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
elseif(solver_type==1)then
  if(myrank==0)print*,'Solver type: CG'
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
elseif(solver_type.eq.SUPERLU)then
  if(myrank==0)print*,'Solver type: SUPERLU'
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
elseif(solver_type.eq.MUMPS)then
  if(myrank==0)print*,'Solver type: MUMPS'
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

!rtol = 1.0d-4
!atol = 1.0d-16
!dtol = 1.0d46
!maxiter = 500
!call KSPSetTolerances(ksp,rtol,atol,dtol,maxiter,ierr)
call KSPSetTolerances(ksp,KSP_RTOL,KSP_ATOL,KSP_DTOL,KSP_MAXITER,ierr)
CHKERRQ(ierr)
if(myrank==0)print*,'ksp2'

!  Set runtime options, e.g.,
!    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.
call KSPSetFromOptions(ksp,ierr)

end subroutine petsc_initialize
!=======================================================

subroutine petsc_set_stiffness_matrix(storekmat)
use math_library_mpi,only:sumscal
use ieee_arithmetic
implicit none
real(kind=kreal),intent(in) :: storekmat(:,:,:)
integer :: i,i_elmt,ielmt,j,ncount,n
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
integer :: imapu(NEDOFU),imapphi(NEDOFPHI)
logical :: isu,isp,isphi
logical,allocatable :: iselmt(:)

real(kind=8) :: kmat(NEDOF,NEDOF),xval
interface
  subroutine issymmetric(myrank,mat)
  implicit none
  integer,parameter :: kreal=selected_real_kind(15)
  integer,intent(in) :: myrank
  real(kind=kreal),intent(in) :: mat(:,:)
  end subroutine issymmetric
end interface
! Set and assemble matrix.
!  - Note that MatSetValues() uses 0-based row and column numbers
!  in Fortran as well as in C (as set here in the array "col").
! stage 0: store all elements

!call MatZeroEntries(Amat,ierr)
!CHKERRQ(ierr)
if(myrank==0)print*,'hi homnath1!'
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
      !print*,'STRANGE IC!!!'
      xval=storekmat(i,j,ielmt)
      if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then
        write(*,*)'ERROR: stiffness matrix has nonfinite value/s!',myrank,ielmt,&
        mat_id(ielmt),xval,minval(abs(storekmat)),maxval(abs(storekmat))
        stop
      endif
      call MatSetValues(Amat,1,ggdof_elmt(irow),1,ggdof_elmt(jcol),           &
      storekmat(i,j,ielmt),ADD_VALUES,ierr)
      CHKERRQ(ierr)
      !print*,'size of kmat:',shape(storekmat_intact_ic)
      !kmat_solid=storekmat_intact_ic(:,:,i_elmt)
      !call issymmetric(myrank,kmat_solid)
    endif
    enddo
  enddo
enddo

call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRQ(ierr)


! check symmetry
call MatDuplicate(Amat,MAT_DO_NOT_COPY_VALUES,Tmat,ierr)
CHKERRQ(ierr)
call MatTranspose(Amat,MAT_INITIAL_MATRIX,Tmat,ierr)
CHKERRQ(ierr)
rval=-1.0
call MatAXPY(Tmat,rval,Amat,SAME_NONZERO_PATTERN,ierr)
call MatNorm(Amat,NORM_FROBENIUS,mnorm,ierr)
if(myrank==0)print*,'Matrix norm:',mnorm
call MatNorm(Tmat,NORM_FROBENIUS,mnorm,ierr)
if(myrank==0)print*,'Symmetry norm:',mnorm
call MatDestroy(Tmat,ierr)
!!!test matrix symmetry is not supported for mataij format
!!!stol=1d-6
!!!call MatIsSymmetric(Amat,stol,mat_symmetry,ierr)
!!!CHKERRQ(ierr)
!!!if(myrank==0)print*,'Matrix symmetry:',mat_symmetry

!call MatSetOption(Amat,MAT_SYMMETRIC,PETSC_TRUE,ierr);
!CHKERRQ(ierr)
!call MatSetOption(Amat,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
!CHKERRQ(ierr)
if(myrank==0)print*,'matrix setting & assembly complete11!'
call sync_process
!call close_process
!if(myrank==0)then
!allocate(cols(2*nzeros_max),vals(2*nzeros_max))
!TMP!if(RECYCLE_STEP_VAL>0.and. .not.exist_recycle)then
!TMP!  call MatGetOwnershipRange(Amat,istart,iend,ierr)
!TMP!  CHKERRQ(ierr)
!TMP!  allocate(cols(2*nzeros_max))
!TMP!  write(ptail,'(i4)')myrank
!TMP!  ptail=adjustl(ptail)
!TMP!  outf_name='recycle_petsc'//trim(ptail)
!TMP!  call sync_process
!TMP!  open(100,file=outf_name,access='stream',form='unformatted',action='write',status='replace')
!TMP!
!TMP!  ngrow=iend-istart
!TMP!  write(100)ngrow
!TMP!  allocate(nnzero_diag(ngrow),nnzero_offdiag(ngrow),stat=ierr)
!TMP!  call check_allocate(ierr,errsrc)
!TMP!  ncount=0
!TMP!  irow=0
!TMP!  do i=istart,iend-1
!TMP!    cols=-1
!TMP!    !call MatGetRow(Amat,i,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);
!TMP!    call MatGetRow(Amat,i,ncols,cols,PETSC_NULL_SCALAR,ierr);
!TMP!    CHKERRQ(ierr)
!TMP!    
!TMP!    ncount=ncount+ncols
!TMP!    !if(myrank==0)print*,'nzeros in',i,'th row:',myrank,ncols,nzeros_max
!TMP!    !if(ncols.gt.nzeros_max)print*,'Very strange!',myrank,ncols,nzeros_max
!TMP!    irow=irow+1
!TMP!    nnzero_diag(irow)=count(cols.ge.ig0.and.cols.le.ig1)
!TMP!    nnzero_offdiag(irow)=ncols-nnzero_diag(irow)
!TMP!    !!endif
!TMP!    call MatRestoreRow(Amat,i,ncols,cols,PETSC_NULL_SCALAR,ierr);
!TMP!    CHKERRQ(ierr)
!TMP!    !call MatRestoreRow(Amat,istart,ncols,cols,vals,ierr);
!TMP!    !call sync_process
!TMP!    !call MatRestoreRow(Amat,i,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);
!TMP!    !CHKERRQ(ierr)
!TMP!  enddo
!TMP!  write(100)nnzero_diag
!TMP!  write(100)nnzero_offdiag
!TMP!  close(100)
!TMP!  deallocate(cols,nnzero_diag,nnzero_offdiag)
!TMP!  !
!TMP!  call sync_process
!TMP!  nsparse1=sumscal(ncount)
!TMP!  if(myrank==0)print*,'total sparse elements:',nsparse1
!TMP!endif
!call close_process
!call MatGetVecs(Amat,vdiag,PETSC_NULL_OBJECT,ierr)
!version >= 3.6.0
call MatCreateVecs(Amat,vdiag,PETSC_NULL_OBJECT,ierr)
call MatGetDiagonal(Amat,vdiag,ierr)
!call MatGetDiagonal(Amat,vdiag,ierr)
call VecGetLocalSize(vdiag,n,ierr)
CHKERRQ(ierr)
!call sync_process
!call close_process
!allocate(diag_array(n))
call VecGetArrayF90(vdiag,diag_array,ierr)
CHKERRQ(ierr)
print*,'NZEROs in diagonal:',myrank,n,count(diag_array==0.),minval(abs(diag_array)),maxval(abs(diag_array))
!!open(1,file='diagonalPETSC',action='write',status='replace')
!!write(1,'(f10.6)')diag_array
!!close(1)
call VecRestoreArrayF90(vdiag,diag_array,ierr)
call sync_process
call VecDestroy(vdiag,ierr)
!endif
end subroutine petsc_set_stiffness_matrix
!=======================================================

!subroutine petsc_set_rigidbody_matrix()
!use ieee_arithmetic
!use math_constants,only:ONE,ZERO
!use global,only:g_coord,ggdof,g_num,nnode,NNDOF
!implicit none
!integer,parameter,dimension(6) :: null_icol=(/1,2,3,4,5,6/)
!integer :: i,i_dof,i_node,j,ncount,n
!integer :: ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)
!
!PetscInt irow,jcol,istart,iend,ndiag,noffdiag,ng,ngrow
!integer :: ncols,ncols_vali,nsparse1
!real(kind=kreal) :: xp,yp,zp
!real(kind=kreal),dimension(6) :: null_row
!character(len=10) :: ptail
!character(len=60) :: outf_name
!Vec   vdiag
!PetscScalar rval
!PetscScalar,pointer :: diag_array(:)
!PetscInt,allocatable :: ibarray(:)
!PetscScalar,pointer :: barray(:)
!PetscReal :: stol,mnorm
!integer :: imapu(NEDOFU),imapphi(NEDOFPHI)
!logical :: isu,isp,isphi
!logical,allocatable :: iselmt(:)
!
!real(kind=8) :: kmat(NEDOF,NEDOF),xval
!! Set and assemble matrix.
!!  - Note that MatSetValues() uses 0-based row and column numbers
!!  in Fortran as well as in C (as set here in the array "col").
!! stage 0: store all elements
!
!!call MatZeroEntries(Amat,ierr)
!!CHKERRQ(ierr)
!
!do i_node=1,nnode
!  xp=g_coord(1,i_node)
!  yp=g_coord(2,i_node)
!  zp=g_coord(3,i_node)
!  do i_dof=1,NNDOF
!    if(ggdof(i_dof,i_node).gt.0)then
!      null_row=ZERO
!      if(i_dof.eq.1)then ! ux
!        null_row(1)= ONE
!        null_row(5)= zp
!        null_row(6)=-yp
!      elseif(i_dof.eq.2)then ! uy
!        null_row(2)= ONE
!        null_row(4)=-zp
!        null_row(6)= xp
!      elseif(i_dof.eq.3)then ! uz
!        null_row(3)= ONE
!        null_row(4)= yp
!        null_row(6)=-xp
!      endif
!      call MatSetValues(Umat,1,ggdof(i_dof,i_node),6,null_icol,null_row,INSERT_VALUES,ierr)
!      CHKERRQ(ierr)
!    endif
!  enddo
!enddo
!call MatAssemblyBegin(Umat,MAT_FINAL_ASSEMBLY,ierr)
!CHKERRQ(ierr)
!call MatAssemblyEnd(Umat,MAT_FINAL_ASSEMBLY,ierr)
!CHKERRQ(ierr)
!
!end subroutine petsc_set_rigidbody_matrix
!!=======================================================
!
!! remove null space for rigid body motion
!subroutine petsc_remove_nullspace
!
!call MatMatTrabsposeMult(Umat,Umat,:wq
!
!
!
!
!end subroutine petsc_remove_nullspace
!!=======================================================

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

!print*,'globalIND:',myrank,istart,iend!,'actualIND:',minval(kgrow_sparse),maxval(kgrow_sparse)
!if(myrank==0)print*,'rload:',maxval(abs(rload))
!call sync_process

! Set exact solution; then compute right-hand-side vector.
!none=-1.0
!ONE=1.0
ZERO=0.0
call VecSet(bvec,ZERO,ierr)
!!do i=1,neq
!!  xval=rload(i)
!!  if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then
!!    stop 'OH GOD!'
!!  endif
!!enddo
!!open(1,file='rload',action='write',status='replace')
!!write(1,'(e12.6)')rload(1:)*1e10
!!close(1)
!if(myrank==0)print*,'OK0',minval(l2gdof(1:)),maxval(l2gdof(1:)),size(load(1:)),neq; call sync_process
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);
!do i=1,neq
!  !call VecSetValues(bvec,neq,l2gdof(1:)-1,load(1:),ADD_VALUES);
!  gdof=l2gdof(i)-1
!  call VecSetValues(bvec,1,gdof,load(i),ADD_VALUES,ierr);
!enddo
!!if(myrank==0)print*,'OK1'; call sync_process
!!! inner core
!!do i_elmt=1,NSPEC_INNER_CORE
!!  if(idoubling_inner_core(i_elmt)==IFLAG_IN_FICTITIOUS_CUBE)cycle
!!  ggdof_elmt=reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
!!  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
!!  !ncount=0; idof=-1; igdof=-1
!!  do i=1,NEDOF
!!    if(ggdof_elmt(i).ge.0)then
!!      !ncount=ncount+1
!!      !idof(ncount)=i
!!      !igdof(ncount)=ggdof_elmt(i)
!!      !if(myrank==0)print*,'hello in IC:',i_elmt,i,j,ggdof_elmt(i),ggdof_elmt(j)
!!      fval=storefvec_inner_core1(i,i_elmt)
!!      call VecSetValues(bvec,1,ggdof_elmt(i),fval,ADD_VALUES,ierr)
!!      CHKERRQ(ierr)
!!    endif
!!  enddo
!!enddo
!!if(myrank==0)print*,'IC:ok'; call sync_process
!!! outer core
!!do i_elmt=1,NSPEC_OUTER_CORE
!!  ggdof_elmt=reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
!!  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
!!  !ncount=0; idof=-1; igdof=-1
!!  do i=1,NEDOF
!!    if(ggdof_elmt(i).ge.0)then
!!      !ncount=ncount+1
!!      !idof(ncount)=i
!!      !igdof(ncount)=ggdof_elmt(i)
!!      !if(myrank==0)print*,'hello in IC:',i_elmt,i,j,ggdof_elmt(i),ggdof_elmt(j)
!!      fval=storefvec_outer_core1(i,i_elmt)
!!      call VecSetValues(bvec,1,ggdof_elmt(i),fval,ADD_VALUES,ierr)
!!      CHKERRQ(ierr)
!!    endif
!!  enddo
!!enddo
!!if(myrank==0)print*,'OC:ok'; call sync_process
!!! crust mantle
!!do i_elmt=1,NSPEC_CRUST_MANTLE
!!  ggdof_elmt=reshape(ggdof_cm(:,inode_elmt_cm(:,i_elmt)),(/NEDOF/))
!!  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
!!  !ncount=0; idof=-1; igdof=-1
!!  do i=1,NEDOF
!!    if(ggdof_elmt(i).ge.0)then
!!      !ncount=ncount+1
!!      !idof(ncount)=i
!!      !igdof(ncount)=ggdof_elmt(i)
!!      !if(myrank==0)print*,'hello in IC:',i_elmt,i,j,ggdof_elmt(i),ggdof_elmt(j)
!!      fval=storefvec_crust_mantle1(i,i_elmt)
!!      call VecSetValues(bvec,1,ggdof_elmt(i),fval,ADD_VALUES,ierr)
!!      CHKERRQ(ierr)
!!    endif
!!  enddo
!!enddo
!!if(myrank==0)print*,'CM:ok'; call sync_process
!!! transition infinite
!!! no load (right-hand side)
!!do i_elmt=1,NSPEC_TRINFINITE
!!  ggdof_elmt=reshape(ggdof_trinf(:,inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
!!  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
!!  !ncount=0; idof=-1; igdof=-1
!!  do i=1,NEDOF
!!    if(ggdof_elmt(i).ge.0)then
!!      !ncount=ncount+1
!!      !idof(ncount)=i
!!      !igdof(ncount)=ggdof_elmt(i)
!!      !if(myrank==0)print*,'hello in IC:',i_elmt,i,j,ggdof_elmt(i),ggdof_elmt(j)
!!      !fval=storefvec_trinfinite1(i,i_elmt)
!!      call VecSetValues(bvec,1,ggdof_elmt(i),ZERO,ADD_VALUES,ierr)
!!      CHKERRQ(ierr)
!!    endif
!!  enddo
!!enddo
!!if(myrank==0)print*,'TRINF:ok'; call sync_process
!!! infinite
!!! no load (right-hand side)
!!do i_elmt=1,NSPEC_INFINITE
!!  ggdof_elmt=reshape(ggdof_inf(:,inode_elmt_inf(:,i_elmt)),(/NEDOF/))
!!  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
!!  !ncount=0; idof=-1; igdof=-1
!!  do i=1,NEDOF
!!    if(ggdof_elmt(i).ge.0)then
!!      !ncount=ncount+1
!!      !idof(ncount)=i
!!      !igdof(ncount)=ggdof_elmt(i)
!!      !if(myrank==0)print*,'hello in IC:',i_elmt,i,j,ggdof_elmt(i),ggdof_elmt(j)
!!      !fval=storefvec_infinite1(i,i_elmt)
!!      call VecSetValues(bvec,1,ggdof_elmt(i),ZERO,ADD_VALUES,ierr)
!!      CHKERRQ(ierr)
!!    endif
!!  enddo
!!enddo
!!if(myrank==0)print*,'INF:ok'; call sync_process

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)
if(myrank==0)print*,'vector setting & assembly complete!'
!!call VecGetArrayF90(bvec,larray,ierr)
!!CHKERRQ(ierr)
!!open(1,file='loadPETSC',action='write',status='replace')
!!write(1,'(e12.6)')larray*1e10
!!close(1)
!!print*,'loadPETSC:',maxval(abs(larray))
!!call VecRestoreArrayF90(bvec,larray,ierr)

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

ZERO=0.0
call VecSet(bvec,ZERO,ierr)
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)
if(myrank==0)print*,'vector setting & assembly complete!'

end subroutine petsc_set_initialguess
!=======================================================

subroutine petsc_solve(sdata,cg_iter)
implicit none
MatNullSpace nullspace
PetscScalar sdata(:)
PetscInt    cg_iter
PetscInt    ireason
PetscScalar,pointer :: larray(:)
!!call MatMult(Amat,bvec,xvec,ierr)
!!call VecGetArrayF90(xvec,larray,ierr)
!!CHKERRQ(ierr)
!!open(1,file='kpPETSC',action='write',status='replace')
!!write(1,'(e12.6)')larray*1e10
!!close(1)
!!call VecRestoreArrayF90(xvec,larray,ierr)
!! copy solution to local array
!call scatter_globalvec(bvec,sdata)

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
if (myrank.lt.1)then
  print*,'converged reason',ireason
  print*,'Iterations:',cg_iter
endif
!if(norm .gt. 1.e-12)then
!  write(*,'(a,e11.4,a,i5)')'Norm of error:',norm,', Iterations:',its
!else
!  write(*,'(a,i5,a)')'Norm of error < 1.e-12, Iterations:',its
!endif

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
!if(myrank==0)print*,'larray:',minval(larray(1:n)),maxval(larray(1:n))
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

end module solver_petsc
!=======================================================
