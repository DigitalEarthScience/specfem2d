module sparse
contains
subroutine prepare_sparse(nup)!,nodalgauss)
use math_constants,only:gtol
use math_library,only:i_uniinv,i8_uniinv,upind_i2j
use math_library_mpi,only:maxscal,minvec,maxvec
use mpi_library,only:check_allocate,sync_process
use global
implicit none
integer,intent(in) :: nup
!real(kind=kreal),intent(in) :: nodalgauss(0:)
integer,parameter :: kint8=selected_int_kind(13)
logical :: ismpi
integer :: errcode
character(len=250) :: errtag

integer :: i,j,i_elmt,ielmt,i_count,n,ncount
integer :: i_node,j_node,iup
integer :: igdof,jgdof
integer :: nmax,ndof
integer :: gdof_elmt(NEDOF),ggdof_elmt(NEDOF)
integer(kind=kint8),allocatable :: ind0(:),iorder(:)
integer,allocatable :: row0(:),col0(:),grow0(:),gcol0(:)

integer :: nglob_ic,nglob_oc,nglob_cm,nglob_trinf,nglob_inf

character(len=12) :: spm
character(len=60) :: fname

integer :: nx,ny,nz,nbyte,off0,gmin,gmax
integer :: i_bool,ibool,i_s,i0,i1,ier,j_proc,istat
logical :: isbig

! counting nonzero elements in offdiagonal portion
integer :: grow,ig0,ig1,ind,neq_part1
integer,allocatable :: nzero_row(:)
logical,allocatable :: iskmat(:,:)


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
if(myid.eq.1) then
  write(stdout,*) 'preparing sparse matrix.'
  write(stdout,*)
endif
 
nmax=nelmt*(NEDOF*NEDOF)
allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax),stat=ierr)
call check_allocate(ierr,errsrc)

allocate(ggdof(NNDOF,nnode),stat=ierr)
call check_allocate(ierr,errsrc)

if(nproc.eq.1)then
  ggdof=gdof
else
  ! read global degrees of freedoms from DATABASE files
  write(spm,'(i10)')mypart
  fname='partition/ggdof_proc'//trim(adjustl(spm))
  open(10,file=fname,action='read',status='old')
  read(10,*)nnode_read
  if(nnode.ne.nnode_read)then
    ! Error
    write(*,*)'ERROR: nnode & nnode_read mismatch!',nnode,nnode_read
    stop
  endif
  read(10,*)ggdof
  close(10)

  allocate(gdof_read(nndof,nnode),gnum_read(nenod,nelmt))
  fname='partition/gdof'//trim(adjustl(spm))
  open(10,file=fname,action='read',status='old')
  read(10,*)nnode_read
  read(10,*)gdof_read
  close(10)
  if(any(.not.(gdof==gdof_read)))print*,mypart,'gdof NOT_EQUAL_TO gdof_read!'

  fname='partition/gnum'//trim(adjustl(spm))
  open(10,file=fname,action='read',status='old')
  read(10,*)nnode_read
  read(10,*)gnum_read
  close(10)
  if(any(.not.(g_num==gnum_read)))print*,mypart,'gnum NOT_EQUAL_TO gnum_read!',count(g_num-gnum_read.ne.0)
endif

! total degrees of freedoms
ngdof=maxscal(maxval(ggdof))

if(myid.eq.1)write(*,'(a,i12)')'Total global degrees of freedom:',ngdof

! precompute ownership range OR partion layout
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! all processors have equal gdofs
  ngz=ng0
  ig0=mypart*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! mypart is 0-based
  if(mypart.le.maxrank0)then
    ngz=ng0
    ig0=mypart*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !mypart.gt.maxrank0
    ngz=ng1
    ig0=np0*ng0+(mypart-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(*,*)'ERROR: illegal value of "np0"!'
  stop
endif

!write(proc_str,'(i10)')mypart
!ofname='partition/gdofsp'//trim(adjustl(proc_str))
!open(22,file=ofname,action='write',status='replace')
!write(22,*)nnode
!write(22,*)gdof
!close(22)

allocate(iseq(0:neq))
iseq=.false.
! stage 0: store all elements
ncount=0

do i_elmt=1,nelmt
  ielmt=i_elmt
  gdof_elmt=reshape(gdof(:,g_num(:,ielmt)),(/NEDOF/))
  ggdof_elmt=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))
  iseq(gdof_elmt)=.true.
  idgdof=gdof_elmt; idggdof=ggdof_elmt
  where(idgdof.gt.0)idgdof=1
  where(idggdof.gt.0)idggdof=1
  if(any((idgdof-idggdof).ne.0))then
    print*,mypart,' gdof mismatch!',ielmt
    print*,'OK:',idgdof-idggdof
    print*,'gdof:',gdof_elmt
    print*,'ggdof:',ggdof_elmt
    stop
  endif
  do i=1,NEDOF
    do j=1,NEDOF
      igdof=gdof_elmt(i)
      jgdof=gdof_elmt(j)
      if(igdof.gt.0.and.jgdof.gt.0)then 
        ncount=ncount+1
        row0(ncount)=igdof
        col0(ncount)=jgdof
        grow0(ncount)=ggdof_elmt(i)
        gcol0(ncount)=ggdof_elmt(j)
      endif
    enddo
  enddo
enddo
call sync_process

if(count(.not.iseq).gt.1)then
  write(*,*)'ERRORSP: some degrees of freedoms missing!',mypart,count(.not.iseq),maxval(gdof)
endif
deallocate(iseq)
neq_actual=maxval(gdof)
call sync_process
! stage 1: assemble duplicates
! sort global indices
allocate(ind0(ncount),iorder(ncount),stat=ierr)
call check_allocate(ierr,errsrc)
ind0=0
do i=1,ncount
  ind0(i)=int(neq,kint8)*(int(row0(i),kint8)-1_kint8)+int(col0(i),kint8)
  if(ind0(i).lt.0)print*,'IMPOSSIBLE:',mypart,neq,row0(i),col0(i),ind0(i)
enddo

call i8_uniinv(ind0,iorder)
nsparse=maxval(iorder)
if(myid.eq.1)print*,'neq:',neq,'Nsparse:',nsparse
call sync_process
allocate(krow_sparse(nsparse),kcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)
allocate(kgrow_sparse(nsparse),kgcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)

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
ndof=nnode*NNDOF
l2gdof(reshape(gdof, (/ndof/)))=reshape(ggdof, (/ndof/))
if(myid.eq.1)then
  do i=1,nsparse
    if(kgrow_sparse(i).ne.l2gdof(krow_sparse(i)).or.kgcol_sparse(i).ne.l2gdof(kcol_sparse(i)))then
      print*,'VERY STRANGE!!!!!'
      stop
    endif
  enddo
endif
l2gdof=l2gdof-1 ! PETSC uses 0 indexing
gmin=minvec(l2gdof(1:))
gmax=maxvec(l2gdof(1:))
if(myid.eq.1)write(*,*)'l2gdof range:',gmin,gmax
call sync_process
if(minval(l2gdof(1:)).lt.0)then
  write(*,*)'ERROR: local-to-global indices are less than 1!'
  stop
endif

allocate(igorder(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call i_uniinv(l2gdof(1:),igorder)

!!sparsity for the convoluted stiffness matrix
!!WARNING: ONLY WORKS FOR SERIAL AND VERY SMALL MODEL
!! count sparse element in all rows for convoluted matrix
!allocate(iskmat(neq,neq),stat=istat)
!if(istat.ne.0)then
!  write(*,*)'Not enough memory for "iskmat"!'
!  stop
!endif
!iskmat=.false.
!do i_node=1,nnode
!  iskmat(i_node,i_node)=.true.
!enddo
!!do i_node=1,nnode-1
!!  do j_node=i_node+1,nnode
!!  ! get the iup index
!!  iup=upind_i2j(nup,nnode,i_node,j_node)
!!  ! compare and set
!!  if(nodalgauss(iup).le.gtol)cycle
!!  iskmat(i_node,j_node)=.true.
!!  iskmat(j_node,i_node)=.true.
!!  enddo
!!enddo
!do i_node=1,nnode
!  do j_node=1,nnode
!    ! get the iup index
!    iup=upind_i2j(nup,nnode,i_node,j_node)
!    ! compare and set
!    if(nodalgauss(iup).le.gtol)cycle
!    iskmat(i_node,j_node)=.true.
!    iskmat(j_node,i_node)=.true.
!  enddo
!enddo
!print*,'nsparse convolution:',count(iskmat),count(iskmat(1,:)),iskmat(1,64),iskmat(1,291),iskmat(1,1111),iskmat(1,1432)
!allocate(nnzero_diagconv(nnode))
!do i=1,nnode
!  nnzero_diagconv(i)=count(iskmat(i,:))
!enddo
!deallocate(iskmat)

end subroutine prepare_sparse
end module sparse
