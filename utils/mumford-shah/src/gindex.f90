! DEVELOPER
!   Hom Nath Gharti
! HISTORY
!   Feb 07,2014; Sep 27,2013
!TODO: numbe of processor can be automatically determined from
!OUTPUT_FILES/values_from_mesher.h. Which can also be used to prompt a message
!if xmeshfem3D hasn't been run yet.
module global_dof
contains
subroutine gindex()
use set_precision
!use earth_constants
use global,only:myrank,nproc,NNDOF,IDOFU
implicit none

! local parameters
integer :: nnode
integer :: i_proc,j_proc
integer :: i_p,j_p
character(len=20) :: snproc

integer :: i,i_dof,ielmt,ix,j,k,i_elmt,i_node
integer :: ig,ispec

! local
logical,allocatable :: isnode(:)

! global
integer,allocatable :: ignode(:)
logical,allocatable :: isgnode(:)

integer,allocatable :: gnf(:,:),gnf_read(:,:)

integer :: ig_end,neq_ext

logical,allocatable :: isgnf(:,:),iseq(:)

integer :: igdof,gnf_end,nibool,nnode_ext
integer,allocatable :: gghost(:,:),ighost(:)
character(len=20) :: fhead
character(len=20) :: spm,spn
character(len=80) :: fname,ifname,pfname

integer,allocatable :: tmpvec(:)
integer,allocatable :: tmpmat(:,:)
integer :: ibool_center
real(kind=kreal) :: xp,yp,zp
logical :: is_center

! ghost partitions                                                               
integer :: ngpart
type ghost_partition                                                             
  integer :: id,nnode                                                     
  integer,dimension(:),allocatable :: inode                                  
end type ghost_partition                                                         
type(ghost_partition),dimension(:),allocatable :: gpart
if(myrank.ne.0)then
  write(*,*)'ERROR: gindex routine must be run on a single processor!'
endif
if(myrank==0)write(*,'(a)')'<< xgindex3D...'

! initialize global indices
ig_end=0 ! global node
gnf_end=0    ! global gdof

! loop through the processors
do i_proc=0,nproc-1
  write(spm,'(i10)')i_proc
  ! read partition information
  pfname='partition/partitioninfo'//trim(adjustl(spm))
  open(22,file=pfname,action='read',status='old')
  read(22,*)nnode
  read(22,*)ngpart
  allocate(gpart(ngpart))
  do i=1,ngpart
    read(22,*)gpart(i)%id
    read(22,*)gpart(i)%nnode
    allocate(gpart(i)%inode(gpart(i)%nnode))
    read(22,*)gpart(i)%inode
  enddo
  close(22)
  
  print*,'Processor:',i_proc,' Neighbours:',ngpart
  
  allocate(isnode(nnode),ignode(nnode),isgnode(nnode))
  isnode=.false.
  ignode=-1;    isgnode=.false.

  !============================================================

  ! global nodal indexing

  ! WARNING: is it correct to put these statements here?
  isgnode=.false.

  fhead='gnode'
  ! copy global indices from preceeding partitions
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc<i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(ighost(nibool))
      read(10,*)ighost(1:nibool)
      close(10,status='delete')
      isgnode(gpart(i)%inode)=.true.
      ignode(gpart(i)%inode)=ighost
      deallocate(ighost)
    endif
  enddo
  print*,'Previous largest node ID:',ig_end

  ! indexify global nodes and store in a region array
  ! inner core
  ig=ig_end 
  do i_node=1,nnode
    if(.not.isgnode(i_node))then
      ig=ig+1
      isgnode(i_node)=.true.
      ignode(i_node)=ig
    endif
  enddo
  
  fhead='gnode'
  ! save global indices for neighbouring partitions
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc>i_proc)then
    write(spn,'(i10)')j_proc
    fname='tmp/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
    open(10,file=fname,action='write',status='replace')
    write(10,*)gpart(i)%nnode
    allocate(tmpvec(gpart(i)%nnode))
    tmpvec=ignode(gpart(i)%inode)
    write(10,*)tmpvec
    deallocate(tmpvec)
    close(10)
    endif
  enddo
  ig_end=maxval(ignode)
  print*,'Largest node ID:',ig_end
  write(spm,'(i10)')i_proc
  fname='partition/gnode_proc'//trim(adjustl(spm))
  open(10,file=fname,action='write',status='replace')
  write(10,*)nnode
  write(10,*)ignode
  close(10)


  ! global indexing of degrees of freedoms
  allocate(gnf(NNDOF,nnode),isgnf(NNDOF,nnode))
  gnf=0
  isgnf=.false.

  ! externally defined Dirichlet BCs
  write(spm,'(i10)')i_proc
  ifname='partition/external_gdof'//trim(adjustl(spm))
  open(22,file=ifname,action='read',status='old')
  read(22,*)nnode_ext
  read(22,*)neq_ext
  if(nnode.ne.nnode_ext)then
    write(*,*)'ERROR: nnode & nnode_ext mismatch!'
    stop
  endif
  allocate(gnf_read(NNDOF,nnode_ext),iseq(0:neq_ext))
  read(22,*)gnf_read
  close(22)

  gnf=gnf_read
  ! copy global indices from preceeding partitions
  fhead='gdof'
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc<i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*)gghost
      close(10)
      isgnf(:,gpart(i)%inode)=.true.
      gnf(:,gpart(i)%inode)=gghost
      deallocate(gghost)
    endif
  enddo
  print*,'Previous largest gnf ID:',gnf_end

  ! test only
  iseq=.false.
  do i_node=1,nnode
    do i_dof=1,NNDOF
      iseq(gnf_read(i_dof,i_node))=.true.
    enddo
  enddo
  if(count(.not.iseq).gt.1)then
    write(*,*)'ERRORSP gindex: some degrees of freedoms missing!',i_proc,neq_ext,nnode,&
    count(iseq),count(.not.iseq),neq_ext,minval(gnf_read),maxval(gnf_read)
  endif

  where(gnf_read>0)gnf_read=1
  igdof=gnf_end ! gdof
  do i_node=1,nnode
    do i_dof=1,NNDOF
      if(gnf(i_dof,i_node).gt.0 .and. .not.isgnf(i_dof,i_node))then
        isgnf(i_dof,i_node)=.true.
        igdof=igdof+1
        gnf(i_dof,i_node)=igdof
      endif
    enddo
  enddo
  deallocate(gnf_read,iseq)
  ! save global degrees of freedom for neighbouring partitions
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc>i_proc)then
    write(spn,'(i10)')j_proc
    fname='tmp/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
    open(10,file=fname,action='write',status='replace')
    write(10,*)gpart(i)%nnode
    allocate(tmpmat(NNDOF,gpart(i)%nnode))
    tmpmat=gnf(:,gpart(i)%inode)
    write(10,*)tmpmat
    deallocate(tmpmat)
    close(10)
    endif
  enddo

  gnf_end=maxval(gnf)
  print*,'Largest gnf ID:',gnf_end
  write(spm,'(i10)')i_proc
  fname='partition/ggdof_proc'//trim(adjustl(spm))
  open(10,file=fname,action='write',status='replace')
  write(10,*)nnode
  write(10,*)gnf
  close(10)

  deallocate(gnf,isgnf)
  do i=1,ngpart                                                                    
    deallocate(gpart(i)%inode)!,gpart(i)%gdof)                                        
  enddo                                                                            
  deallocate(gpart) 
  deallocate(isnode,ignode,isgnode)
enddo !i_proc=0,nproc 
  
if(myrank==0)write(*,'(a)')' SUCCESS!'
  
end subroutine gindex
end module global_dof

