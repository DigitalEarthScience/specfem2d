module blur
contains
! this subroutine computed image blur kernel, derivative and gradien of the
! kernel
! REVISION
!   HNG, Jul 12,2011; ; HNG, Apr 09,2010
subroutine blur_kernel(ismpi,s,x,h,dh_ds,dgradh2_ds,errcode,errtag)
use global
use math_library,only:norm
use math_constants, only : one,half,two,pi
implicit none
logical,intent(in) :: ismpi
real(kind=kreal),intent(in) :: s
real(kind=kreal),dimension(ndim),intent(in) :: x
real(kind=kreal),intent(out) :: h,dh_ds,dgradh2_ds
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

real(kind=kreal) :: s2,x2
integer :: ipart ! partition ID

errtag="ERROR: unknown!"
errcode=-1

ipart=myid-1 ! partition ID starts from 0
!if(ismpi)then
!  write(format_str,*)ceiling(log10(real(nproc)+1))
!  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
!  write(ptail, fmt=format_str)'_proc',ipart
!else
!  ptail=""
!endif

x2=sum(x*x)
s2=s*s

! kernel
h=(half/(pi*s2))*exp(-half*x2/s2)

! derivative
dh_ds=-(one/(pi*s**3))*(one-half*x2/s2)*exp(-half*x2/s2)

! derivative of the kernel
dgradh2_ds=-two*x2/(pi**2*s**9)*(one-0.25_kreal*x2/s2)*exp(-x2/s2)

errcode=0
return


end subroutine blur_kernel
!===============================================================================

! compute gauss kernel factor for convolution
subroutine compute_gaussian_kernel(sigma,nodalgauss)
use math_constants,only:half,pi,zero
use global
implicit none
real(kind=kreal),intent(in) :: sigma
real(kind=kreal),intent(out) :: nodalgauss(0:)

integer :: i,j,iup
real(kind=kreal) :: halfinvs2,fac,xsq
real(kind=kreal) :: xi(ndim),xj(ndim)

halfinvs2=half/(sigma*sigma)
fac=halfinvs2/pi

! diagonal elements have index 0
nodalgauss(0)=fac
iup=0
do i=1,nnode-1
  xi=g_coord(:,i)
  do j=i+1,nnode
    iup=iup+1
    xj=g_coord(:,j)
    xsq=sum((xi-xj)*(xi-xj))
    nodalgauss(iup)=fac*exp(-halfinvs2*xsq)
  enddo
enddo
write(stdout,'(a)')'complete!'
end subroutine compute_gaussian_kernel
!===============================================================================

subroutine convolve_with_gaussian(ismpi,gnod,sigma,nodalf,nodalfg,errcode,errtag)
use global
use math_library,only:determinant,norm
use math_constants, only : one,half,two,pi,zero
implicit none
logical,intent(in) :: ismpi
integer,intent(in) :: gnod(4)!geometrical nodes (corner nodes) per element
real(kind=kreal),intent(in) :: sigma
real(kind=kreal),dimension(nnode),intent(in) :: nodalf ! nodal data to be convolved
real(kind=kreal),dimension(nnode),intent(out) :: nodalfg ! nodal convolved data
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt,i_gll,i_node
integer :: num(nenod)
real(kind=kreal) :: fac,halfinvs2,x2
real(kind=kreal) :: conv,gconv
real(kind=kreal) :: detjac,gauss
real(kind=kreal) :: coord(ndim,nenod),f(nenod),x(ndim),xp(ndim),xcenter(ndim)
real(kind=kreal) :: jac(ndim,ndim)
real(kind=kreal),dimension(nnode) :: nodalgauss
integer :: ipart ! partition ID

errtag="ERROR: unknown!"
errcode=-1

ipart=myid-1 ! partition ID starts from 0
!if(ismpi)then
!  write(format_str,*)ceiling(log10(real(nproc)+1))
!  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
!  write(ptail, fmt=format_str)'_proc',ipart
!else
!  ptail=""
!endif

halfinvs2=half/(sigma*sigma)
fac=halfinvs2/pi
nodalfg=zero
do i_node=1,nnode
  xp= g_coord(:,i_node)
 
  conv=zero; gconv=zero 
  do i_elmt=1,nelmt
    num=g_num(:,i_elmt)
    coord=g_coord(:,num)
    xcenter=half*(coord(:,1)+coord(:,ngll))
    if(halfinvs2*sum((xp-xcenter)*(xp-xcenter)).gt.30.0_kreal)then
      !print*,'skipped element:',i_elmt
      cycle
    endif
    f=nodalf(num)
    do i_gll=1,ngll
      x=g_coord(:,num(i_gll))
      x2=sum((xp-x)*(xp-x))

      gauss=fac*exp(-halfinvs2*x2)
      !if(abs(gauss).lt.1.0e-40_kreal)cycle
      conv=conv+f(i_gll)*gauss*storejw(i_gll,i_elmt)
      gconv=gconv+gauss
    end do ! i_gll=1,ngll
  enddo ! i_elmt
  nodalfg(i_node)=conv
  nodalgauss(i_node)=gconv
  !print*,'node:',i_node,'/',nnode,' finished!'
enddo
!!normalize
!where(nodalgauss.ne.zero)nodalfg=nodalfg/nodalgauss
print*,'zeros in nodal gauss:',count(nodalgauss==zero),minval(abs(nodalgauss))

errcode=0
return


end subroutine convolve_with_gaussian
!===============================================================================

subroutine convolve_with_gaussian_precomp(ismpi,gnod,nup,nodalgauss,nodalf,nodalfg,errcode,errtag)
use global
use math_library,only:determinant,norm,upinds_i2js
use math_constants, only : gtol,one,half,two,pi,zero
implicit none
logical,intent(in) :: ismpi
integer,intent(in) :: gnod(4)!geometrical nodes (corner nodes) per element
integer,intent(in) :: nup
real(kind=kreal),intent(in) :: nodalgauss(0:) ! nodal data to be convolved
real(kind=kreal),dimension(nnode),intent(in) :: nodalf ! nodal data to be convolved
real(kind=kreal),dimension(nnode),intent(out) :: nodalfg ! nodal convolved data
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt,i_gll,i_node,j_node
integer :: i,j,iup
integer :: num(nenod),iups(ngll)
real(kind=kreal) :: fac,halfinvs2,x2
real(kind=kreal) :: conv,gconv
real(kind=kreal) :: detjac,gauss
real(kind=kreal) :: coord(ndim,nenod),f(nenod),x(ndim),xp(ndim),xcenter(ndim)
real(kind=kreal) :: jac(ndim,ndim)
real(kind=kreal),dimension(nnode) :: nodalgauss_sum
integer :: ipart ! partition ID

errtag="ERROR: unknown!"
errcode=-1

ipart=myid-1 ! partition ID starts from 0
!if(ismpi)then
!  write(format_str,*)ceiling(log10(real(nproc)+1))
!  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
!  write(ptail, fmt=format_str)'_proc',ipart
!else
!  ptail=""
!endif
nodalgauss_sum=zero
nodalfg=zero
do i_node=1,nnode
  xp= g_coord(:,i_node)
 
  conv=zero; gconv=zero 
  do i_elmt=1,nelmt
    num=g_num(:,i_elmt)

    ! right index list for nodalgauss
    iups=upinds_i2js(nup,nnode,i_node,ngll,num)
    if(maxval(abs(nodalgauss(iups))).lt.gtol)cycle
    f=nodalf(num)
    do i_gll=1,ngll
      !x=g_coord(:,num(i_gll))
      !x2=sum((xp-x)*(xp-x))
      !
      !j_node=num(i_gll)

      !! compute index
      !if(i_node.eq.j_node)then
      !  iup=0
      !else
      !  if(i_node.gt.j_node)then
      !    i=j_node
      !    j=i_node
      !  else
      !    i=i_node
      !    j=j_node
      !  endif
      !  !print*,i_node,j_node
      !  iup=nup - (nnode-i)*(nnode-i+1)/2 +j - i
      !endif
      iup=iups(i_gll)
      if(abs(nodalgauss(iup)).lt.gtol)cycle
     
      conv=conv+f(i_gll)*nodalgauss(iup)*storejw(i_gll,i_elmt)
      gconv=gconv+nodalgauss(iup)*storejw(i_gll,i_elmt)
    end do ! i_gll=1,ngll
  enddo ! i_elmt
  ! normalize
  nodalfg(i_node)=conv!/gconv
  nodalgauss_sum(i_node)=gconv
  !print*,'node:',i_node,'/',nnode,' finished!'
enddo
!!normalize
!where(nodalgauss_sum.ne.zero)nodalfg=nodalfg/nodalgauss_sum
print*,'zeros in nodal gauss:',count(nodalgauss==zero),minval(abs(nodalgauss))

errcode=0
return


end subroutine convolve_with_gaussian_precomp
!===============================================================================
end module blur

