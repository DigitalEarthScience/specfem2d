module mesh
contains
! this subroutine reads the input information from a structured ASCII text file
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO:
!   - prompt warning or error for unknown argument/s
subroutine read_mesh(ismpi,inp_fname,errcode,errtag,ispartmesh)
use global
use math_constants,only:zero,zerotol
use string_library
implicit none

integer :: i
character(len=*),intent(in) :: inp_fname
logical,intent(in) :: ismpi
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
logical,optional,intent(in) :: ispartmesh
character(len=250) :: line
character(len=800) ::tag
character(len=80) :: strval,token
character(len=1) :: tmp_char
character(len=80),dimension(50) :: args
character(len=20) :: ptail_inp
integer :: id,ind,ios,narg,slen

integer :: i_elmt,ielmt,i_node,inode,i_mat,imat,tmp_nelmt,tmp_nnode !,mat_domain

character(len=20) :: format_str
character(len=250) :: fname
character(len=150) :: data_path

integer :: ipart,nproc_inp ! partition ID
real(kind=4),allocatable :: coordgll(:,:)

errtag="ERROR: unknown!"
errcode=-1

ipart=myid-1 ! partition ID starts from 0

! reading mesh information
if(myid==1)write(*,'(a)',advance='no')'reading mesh information...'

if(ismpi.and.nproc.gt.1)then
  ptail_inp=trim(ptail)
else
  ptail_inp=""
endif
! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

! read connectivity
fname=trim(data_path)//trim(confile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',form='unformatted',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11)nelmt
print*,nelmt,nenod,ngll
allocate(g_num(nenod,nelmt))

read(11)g_num
close(11)

nnode=maxval(g_num)
! read coordinates information
allocate(g_coord(ndim,nnode))
allocate(coordgll(ngll,nelmt))
do i=1,ndim
  fname=trim(data_path)//trim(coordfile(i))//trim(ptail_inp)
  !print*,fname
  open(unit=11,file=trim(fname),status='old',action='read',form='unformatted',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  read(11)coordgll
  do i_elmt=1,nelmt
    g_coord(i,g_num(:,i_elmt))=coordgll(:,i_elmt)
  enddo
enddo
deallocate(coordgll)
close(11)

errcode=0
if(myid==1)write(*,*)'complete!'

end subroutine read_mesh
end module mesh
