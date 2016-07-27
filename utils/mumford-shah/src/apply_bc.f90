module bc
contains
subroutine modify_gdof(lgdof,neq,errcode,errtag)
use global,only:nndof,nnode,g_num,gdof,myid
implicit none
integer,dimension(nndof,nnode),intent(inout) :: lgdof
integer,intent(out) :: neq
integer,intent(out) :: errcode
character(len=20) :: proc_str
character(len=80) :: ofname
character(len=250),intent(out) :: errtag
integer :: i,j,ipart

ipart=myid-1
write(proc_str,'(i10)')ipart

errtag="ERROR: unknown!"
errcode=-1
! compute modified lgdof
neq=0
do j=1,ubound(lgdof,2)
  do i=1,ubound(lgdof,1)
    if(lgdof(i,j)/=0)then
      neq=neq+1
      lgdof(i,j)=neq
    endif
  enddo
enddo

! compute nodal to global
errcode=0
return
end subroutine modify_gdof
!===============================================================================
end module bc

