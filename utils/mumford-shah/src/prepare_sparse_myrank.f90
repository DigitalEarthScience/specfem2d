module sparse
contains
subroutine prepare_sparse()
use math_library,only:i_uniinv,i8_uniinv
use math_library_mpi,only:maxscal,minvec,maxvec
use mpi_library,only:check_allocate,sync_process
use global
implicit none
integer,parameter :: kint8=selected_int_kind(13)
logical :: ismpi
integer :: errcode
character(len=250) :: errtag

integer :: i,j,i_elmt,ielmt,i_count,n,ncount
integer :: igdof,jgdof
integer :: nmax,ndof
integer :: gdof_elmt(NEDOF),ggdof_elmt(NEDOF)
integer(kind=kint8),allocatable :: ind0(:),iorder(:)
integer,allocatable :: row0(:),col0(:),grow0(:),gcol0(:)

integer :: nglob_ic,nglob_oc,nglob_cm,nglob_trinf,nglob_inf

character(len=12) :: spm
character(len=60) :: fname

integer :: nx,ny,nz,nbyte,off0,gmin,gmax
integer :: i_bool,ibool,i_s,i0,i1,ier,j_proc
logical :: isbig

! counting nonZERO elements in offdiagonal portion
integer :: grow,ig0,ig1,ind,neq_part1
integer,allocatable :: nZERO_rowoff1(:)


integer :: irow,jcol,neq_actual
integer :: idgdof(NEDOF),idggdof(NEDOF)
logical :: isu,isphi
integer,allocatable :: gdof_read(:,:),gnum_read(:,:)
integer :: itest(6),otest(6)
integer :: ngz,ng0,ng1,np0,maxrank0,nnode_read
logical,allocatable :: iseq(:)
integer,allocatable :: igorder(:)
integer :: ierr, mypart
character(len=250) :: myfname=" => prepare_runtime.f90"
character(len=500) :: errsrc

character(len=20) :: proc_str
character(len=80) :: ofname
errsrc=trim(myfname)//' => prepare_timerun_sparse'

mypart = myid-1 !for consistency with myrank file numbering
ismpi=.true.
if(mypart==0) then
  write(stdout,*) 'preparing sparse matrix.'
  write(stdout,*)
endif
 
nmax=nelmt*(NEDOF*NEDOF)
allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax),stat=ierr)
call check_allocate(ierr,errsrc)
if(mypart==0)print*,NEDOF,nmax
call sync_process
! read global degrees of freedoms from DATABASE files
! inner core
write(spm,'(i10)')mypart
call sync_process
fname='partition/ggdof_proc'//trim(adjustl(spm))
if(myrank==2)print*,trim(fname)
open(10,file=fname,action='read',status='old')
read(10,*)nnode_read
if(nnode.ne.nnode_read)then
  ! Error
  write(*,*)'ERROR: nnode & nnode_read mismatch!',nnode,nnode_read
  stop
endif
allocate(ggdof(NNDOF,nnode_read),stat=ierr)
call check_allocate(ierr,errsrc)
read(10,*)ggdof
close(10)

allocate(gdof_read(nndof,nnode),gnum_read(nenod,nelmt))
fname='partition/gdof'//trim(adjustl(spm))
open(10,file=fname,action='read',status='old')
read(10,*)nnode_read
read(10,*)gdof_read
close(10)
if(any(.not.(gdof==gdof_read)))print*,myrank,'gdof NOT_EQUALT_TO gdof_read!'

fname='partition/gnum'//trim(adjustl(spm))
open(10,file=fname,action='read',status='old')
read(10,*)nnode_read
read(10,*)gnum_read
close(10)
if(any(.not.(g_num==gnum_read)))print*,myrank,'gnum NOT_EQUALT_TO gnum_read!',count(g_num-gnum_read.ne.0)
! total degrees of freedoms
ngdof=maxscal(maxval(ggdof))

if(myrank==0)write(*,'(a,i12)')'Total global degrees of freedom:',ngdof

if(myrank==0)print*,'nedofs:',nedof

! precompute ownership range OR partion layout
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! ng0=ng1
! all processors have equal gdofs
  ngz=ng0
  ig0=myrank*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! myrank is 0-based
  if(myrank.le.maxrank0)then
    ngz=ng0
    ig0=myrank*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !myrank.gt.maxrank0
    ngz=ng1
    ig0=np0*ng0+(myrank-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(*,*)'ERROR: illegal value of "np0"!'
  stop
endif

if(myrank==0)print*,'prepare_runtime OK0:',ng0,ng1,ngz,ig0,ig1


write(proc_str,'(i10)')myrank
ofname='partition/gdofsp'//trim(adjustl(proc_str))
open(22,file=ofname,action='write',status='replace')
write(22,*)nnode
write(22,*)gdof
close(22)

allocate(iseq(0:neq))
iseq=.false.
! stage 0: store all elements
ncount=0

do i_elmt=1,nelmt
  ielmt=i_elmt
  gdof_elmt=reshape(gdof(:,g_num(:,ielmt)),(/NEDOF/))
  ggdof_elmt=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))
  if(myrank==0)then
    if(i_elmt.eq.23)print*,'in prepare_sparse:',myrank,i_elmt,neq,gdof_elmt
    !if(any(gdof_elmt.eq.1365))stop 'OMG!'
  endif
  iseq(gdof_elmt)=.true.
  if(myrank==0.and.i_elmt==1)print*,'hello gdof_elmt:',gdof_elmt
  if(myrank==0.and.i_elmt==1)print*,'hello ggdof_elmt:',ggdof_elmt
  if(myrank==2.and.ielmt==23)then
    print*,'HC:',gdof_elmt
    print*,'HC1:',ggdof_elmt
  endif
  idgdof=gdof_elmt; idggdof=ggdof_elmt
  where(idgdof.gt.0)idgdof=1
  where(idggdof.gt.0)idggdof=1
  if(myrank==0)print*,'difference:',idgdof-idggdof
  if(any((idgdof-idggdof).ne.0))then
    print*,myrank,' gdof mismatch!',ielmt
    print*,'OK:',idgdof-idggdof
    print*,'gdof:',gdof_elmt
    print*,'ggdof:',ggdof_elmt
    stop
  endif
  if(myrank==0)print*,'count dof ic:',count(idgdof==1),count(idggdof==1)
  !call sync_process
  !call close_process
  do i=1,NEDOF
    do j=1,NEDOF
      igdof=gdof_elmt(i)
      jgdof=gdof_elmt(j)
      if(igdof.gt.0.and.jgdof.gt.0)then !.and.storekmat(i,j,i_elmt).ne.0.0_kreal)then
        ncount=ncount+1
        row0(ncount)=igdof
        col0(ncount)=jgdof
        grow0(ncount)=ggdof_elmt(i)
        gcol0(ncount)=ggdof_elmt(j)
      endif
    enddo
  enddo
enddo
if(myrank==0)print*,'kmat dONE!',ncount,NEDOF
print*,'oh GODSP!',myrank,neq,maxval(gdof)
call sync_process
do i=1,neq
  if(.not.iseq(i))then
    print*,'wow:',myrank,i
  endif
enddo

if(count(.not.iseq).gt.1)then
  write(*,*)'ERRORSP: some degrees of freedoms missing!',myrank,count(.not.iseq),maxval(gdof)
  !stop
endif
deallocate(iseq)
neq_actual=maxval(gdof)
print*,myrank,'neq:',neq,neq_actual,maxval(gdof),maxval(ggdof)
call sync_process
! stage 1: assemble duplicates
! sort global indices
allocate(ind0(ncount),iorder(ncount),stat=ierr)
call check_allocate(ierr,errsrc)
ind0=0
!ind0=neq*(row0(1:ncount)-1)+col0(1:ncount)
do i=1,ncount
  ind0(i)=int(neq,kint8)*(int(row0(i),kint8)-1_kint8)+int(col0(i),kint8)
  if(ind0(i).lt.0)print*,'IMPOSSIBLE:',myrank,neq,row0(i),col0(i),ind0(i)
enddo

call i8_uniinv(ind0,iorder)
nsparse=maxval(iorder)
if(myrank==0)print*,'test point:',minval(ind0),maxval(ind0),maxval(gdof),minval(row0),minval(col0),maxval(row0),maxval(col0)
if(myrank==0)print*,'neq:',neq,'Nsparse:',nsparse,ncount,minval(row0(1:ncount)),minval(col0(1:ncount)),count(row0>0),count(col0>0)
call sync_process
!!  allocate(kmat_sparse(nsparse),krow_sparse(nsparse),kcol_sparse(nsparse))
allocate(krow_sparse(nsparse),kcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)
allocate(kgrow_sparse(nsparse),kgcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)
!if(myrank==0)print*,'OK02',nmax,n
!call sync_process

!kmat_sparse1=0.0_kreal
krow_sparse=-1
kcol_sparse=-1
kgrow_sparse=-1
kgcol_sparse=-1
do i_count=1,ncount!nmax
  krow_sparse(iorder(i_count))=row0(i_count)
  kcol_sparse(iorder(i_count))=col0(i_count)
  kgrow_sparse(iorder(i_count))=grow0(i_count)
  kgcol_sparse(iorder(i_count))=gcol0(i_count)
enddo
!if(myrank==0)print*,'OK03'
!call sync_process
if(minval(krow_sparse).lt.1.or.minval(kcol_sparse).lt.1.or.                    &
minval(kgrow_sparse).lt.1.or.minval(kgcol_sparse).lt.1)then
  write(*,*)'ERROR: local and global indices are less than 1!',                &
  minval(krow_sparse),minval(kcol_sparse),minval(kgrow_sparse),                &
  minval(kgcol_sparse)
  stop
endif

deallocate(row0,col0,grow0,gcol0,ind0,iorder)

! local DOF to global DOF mapping
allocate(l2gdof(0:neq),stat=ierr)
call check_allocate(ierr,errsrc)
l2gdof=-1
!l2gdof(gdof_ic)=ggdof_ic(1,:)
!l2gdof(gdof_oc)=ggdof_oc(1,:)
!l2gdof(gdof_cm)=ggdof_cm(1,:)
!l2gdof(gdof_trinf)=ggdof_trinf(1,:)
!l2gdof(gdof_inf)=ggdof_inf(1,:)
ndof=nnode*NNDOF
l2gdof(reshape(gdof, (/ndof/)))=reshape(ggdof, (/ndof/))
if(myrank==0)then
  do i=1,nsparse
    !print*,'testing:',kgrow_sparse(i),l2gdof(krow_sparse(i)),kgcol_sparse(i),l2gdof(kcol_sparse(i))
    if(kgrow_sparse(i).ne.l2gdof(krow_sparse(i)).or.kgcol_sparse(i).ne.l2gdof(kcol_sparse(i)))then
      print*,'VERY STRANGE!!!!!'
      stop
    endif
  enddo
endif
l2gdof=l2gdof-1 ! PETSC uses 0 indexing
gmin=minvec(l2gdof(1:))
gmax=maxvec(l2gdof(1:))
if(myrank==0)write(*,*)'l2gdof range:',gmin,gmax
call sync_process
if(minval(l2gdof(1:)).lt.0)then
  write(*,*)'ERROR: local-to-global indices are less than 1!'
  stop
endif

allocate(igorder(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call i_uniinv(l2gdof(1:),igorder)

end subroutine prepare_sparse

!TMP!!
!TMP!!-------------------------------------------------------------------------------------------------
!TMP!!
!TMP!subroutine prepare_timerun_freesurface()
!TMP!use specfem_par,only:kreal=>CUSTOM_REAL,NGLLX,NGLLY,NGLLXY
!TMP!use specfem_par_crustmantle
!TMP!implicit none
!TMP!integer :: i,j,i_spec,inum,ispec,nsnode,nsnode_all
!TMP!integer :: i1,i2,j1,j2,iquad
!TMP!logical,allocatable :: isnode(:)
!TMP!integer,allocatable :: inode_order(:),nodelist(:),nodelist_spec(:,:,:)
!TMP!!real(kind=kreal),allocatable :: xstore_freesurf(:),ystore_freesurf(:),zstore_freesurf(:)
!TMP!
!TMP!! inteface for routines with assumed shaped arrays not in the module
!TMP!interface
!TMP!subroutine i_uniinv(XDONT, IGOEST)
!TMP!  implicit none
!TMP!  integer,intent(in)  :: XDONT(:)
!TMP!  integer,intent(out) :: IGOEST(:)
!TMP!end subroutine i_uniinv
!TMP!end interface
!TMP!
!TMP!! surface node list    
!TMP!nsnode_all=NSPEC2D_TOP_CM*NGLLXY
!TMP!allocate(nodelist(nsnode_all),inode_order(nsnode_all))
!TMP!allocate(nodelist_spec(NGLLX,NGLLY,NSPEC2D_TOP_CM))
!TMP!inum=0    
!TMP!do i_spec=1,NSPEC2D_TOP_CM      
!TMP!  ispec=ibelm_top_crust_mantle(i_spec)
!TMP!  do j=1,NGLLY
!TMP!    do i=1,NGLLX
!TMP!      inum=inum+1
!TMP!      nodelist_spec(i,j,i_spec)=inum
!TMP!      nodelist(inum)=ibool_crust_mantle(i,j,NGLLZ,ispec)
!TMP!    enddo
!TMP!  enddo        
!TMP!enddo
!TMP!!print*,'hi0'
!TMP!!call sync_process
!TMP!!call close_process!exit_MPI(myrank,'what!')
!TMP!if(nsnode_all.ne.inum)then
!TMP!  write(*,*)'ERROR: total number of surface nodes mismatch!'
!TMP!  stop
!TMP!endif
!TMP!
!TMP!call i_uniinv(nodelist,inode_order)
!TMP!!call close_process!exit_MPI(myrank,'what!')
!TMP!nsnode=maxval(inode_order)
!TMP!nnode_freesurf=nsnode
!TMP!!print*,'hi0 OK',nsnode_all,nsnode
!TMP!!call sync_process
!TMP!allocate(ibool_freesurf(nsnode),isnode(nsnode),coord_freesurf(NDIM,nsnode))
!TMP!!!,mirxs(NDIM,nsnode,NGLLZ),mirxs1(NDIM,nsnode),mirxs2(NDIM,nsnode))
!TMP!!print*,'hi0'
!TMP!!call sync_process
!TMP!!call exit_MPI(myrank,'what!')
!TMP!isnode=.false.
!TMP!! assign surface nodes: xs
!TMP!do i=1,nsnode_all
!TMP!  if(.not.isnode(inode_order(i)))then
!TMP!    coord_freesurf(1,inode_order(i))=xstore_crust_mantle(nodelist(i))
!TMP!    coord_freesurf(2,inode_order(i))=ystore_crust_mantle(nodelist(i))
!TMP!    coord_freesurf(3,inode_order(i))=zstore_crust_mantle(nodelist(i))
!TMP!    isnode(inode_order(i))=.true.
!TMP!    ibool_freesurf(inode_order(i))=nodelist(i)
!TMP!  endif
!TMP!enddo
!TMP!deallocate(isnode)
!TMP!NSUBQUAD_FREESURF=NSPEC2D_TOP_CM*(NGLLX-1)*(NGLLY-1)
!TMP!allocate(connect_freesurf(4,NSUBQUAD_FREESURF))
!TMP!iquad=0
!TMP!do i_spec=1,NSPEC2D_TOP_CM      
!TMP!  do j=1,NGLLY-1
!TMP!    do i=1,NGLLX-1
!TMP!      iquad=iquad+1
!TMP!      i1=i; i2=i+1
!TMP!      j1=j; j2=j+1
!TMP!      connect_freesurf(1,iquad)=inode_order(nodelist_spec(i1,j1,i_spec))
!TMP!      connect_freesurf(2,iquad)=inode_order(nodelist_spec(i2,j1,i_spec))
!TMP!      connect_freesurf(3,iquad)=inode_order(nodelist_spec(i2,j2,i_spec))
!TMP!      connect_freesurf(4,iquad)=inode_order(nodelist_spec(i1,j2,i_spec))
!TMP!    enddo
!TMP!  enddo        
!TMP!enddo
!TMP!deallocate(nodelist_spec)
!TMP!end subroutine prepare_timerun_freesurface
!TMP!
!TMP!!
!TMP!!-------------------------------------------------------------------------------------------------
!TMP!!

end module sparse
