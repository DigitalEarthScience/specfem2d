! this modonins math constants
! math parameters
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module math_constants
use set_precision
implicit none
real(kind=kreal),parameter :: zero=0.0_kreal,half=0.5_kreal,one=1.0_kreal,     &
two=2.0_kreal
real(kind=kreal),parameter :: pi=3.141592653589793_kreal
real(kind=kreal),parameter :: deg2rad=pi/180.0_kreal,rad2deg=180.0_kreal/pi

! tolerance value for zero
real(kind=kreal),parameter :: inftol=1.0e32_kreal,zerotol = 1.0e-12_kreal

! tolerance for gaussian convolution
real(kind=kreal),parameter :: gtol=1e-16_kreal
end module math_constants
!=======================================================

! this module contains math routines
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module math_library
use set_precision
use math_constants
contains
function cross_product(v1,v2) result(v)
real(kind=kreal),dimension(3),intent(in) :: v1,v2
real(kind=kreal),dimension(3) :: v

! cross product
v(1)=v1(2)*v2(3)-v2(2)*v1(3)
v(2)=v2(1)*v1(3)-v1(1)*v2(3)
v(3)=v1(1)*v2(2)-v2(1)*v1(2)

return
end function cross_product
!=======================================================
! this function computes normal to the plane formed by three points
!        1-------2
!         \     /
!          \   /
!           \ /
!            0
function get_normal(x0,x1,x2) result(nx)
real(kind=kreal),dimension(3),intent(in) :: x0,x1,x2
real(kind=kreal),dimension(3) :: nx
real(kind=kreal),dimension(3) :: v1,v2
real(kind=kreal) :: norm

! two vectors
v1=x1-x0
v2=x2-x0
! cross product
nx(1)=v1(2)*v2(3)-v2(2)*v1(3)
nx(2)=v2(1)*v1(3)-v1(1)*v2(3)
nx(3)=v1(1)*v2(2)-v2(1)*v1(2)
norm=sqrt(dot_product(nx,nx))
if(norm<=0.0_kreal)then
  write(*,*)'ERROR: undefined normal!'
  stop
endif
! unit normal
nx=nx/norm


return
end function get_normal
!=======================================================

! compute distance between two points in a n-dimensional space                   
function distance(x1,x2,n) result(r)                                             
implicit none                                                                    
integer,intent(in) :: n                                                          
real(kind=kreal),intent(in) :: x1(n),x2(n)                                       
real(kind=kreal) :: dx(n),r                                                      
                                                                                 
dx=x1-x2                                                                         
r=sqrt(sum(dx*dx))                                                               
return                                                                           
end function distance                                                            
!=======================================================

function angle(x,y) result(theta)
!
! this function calculates angle in radian
!
use math_constants,only:pi,zero,half
implicit none
real(kind=kreal),intent(in) :: x,y
real(kind=kreal) :: theta

! origin
if(x.eq.zero.and.y.eq.zero)then
  theta=zero
  return
endif

! x-axis
if(y.eq.zero)then
  if(x.gt.zero)then
    theta=zero
    return
  else
    theta=pi
    return
  endif
endif

! y-axis
if(x.eq.zero)then
  if(y.gt.zero)then
    theta=half*pi
    return
  else
    theta=3.0_kreal*half*pi
    return
  endif
endif

theta=atan(abs(y/x))
! QI
if(x.gt.zero.and.y.gt.zero)then  
  return
! QII
elseif(x.lt.zero.and.y.gt.zero)then
  theta=pi-theta
  return
!QIII
elseif(x.lt.zero.and.y.lt.zero)then
  theta=pi+theta
  return
!QIV
elseif(x.gt.zero.and.y.lt.zero)then
  theta=2.0_kreal*pi-theta
  return
endif
end function angle
!=======================================================

function norm(x) result(l2n)
!
! this function calculates the l2 norm of vector x
!
implicit none
real(kind=kreal),intent(in) :: x(:)
real(kind=kreal)::l2n
l2n=sqrt(sum(x**2))
return
end function norm
!=======================================================

! the function below determines whether the given point (xp) is on or inside
! the hexahedron defined by its corner points (xcorner)
!          8______________7
!         /|             /|
!        / |            / |
!       /  |           /  |
!     5/___|_________6/   |
!     |    |         |    |
!     |    4---------|----3
!     |   /          |   /
!     |  /           |  /
!     | /            | /
!     1/_____________2/
function IsPointInHexahedron(ielmt,xcorner,xp) result(isinside)
implicit none
integer :: ielmt
real(kind=kreal),intent(in) :: xcorner(3,8),xp(3)
logical :: isinside

real(kind=kreal),parameter :: tol=1e-16_kreal
integer :: fnode3(3,6)
integer :: i_face,ndim,ngnode
real(kind=kreal) :: cosine,normp
real(kind=kreal) :: dxp(3),normvec(3),uvecp(3),x0(3),x1(3),x2(3)

! check number of nodes
ngnode=ubound(xcorner,2)
if(ngnode.ne.8)then
  print*,'ERROR: number of nodes must be 8 for IsPointInHexahedron(..)!'
  stop
endif
ndim=ubound(xcorner,1)
if(ndim.ne.3)then
  print*,'ERROR: dimenionsion must be 3 for IsPointInHexahedron(..)!'
  stop
endif
! choose three nodes for each face such that the normal is always outward
fnode3(1,1)=1; fnode3(2,1)=2; fnode3(3,1)=5 ! front
fnode3(1,2)=2; fnode3(2,2)=3; fnode3(3,2)=6 ! right
fnode3(1,3)=3; fnode3(2,3)=4; fnode3(3,3)=7 ! back
fnode3(1,4)=4; fnode3(2,4)=1; fnode3(3,4)=8 ! left
fnode3(1,5)=1; fnode3(2,5)=4; fnode3(3,5)=2 ! bottom
fnode3(1,6)=5; fnode3(2,6)=6; fnode3(3,6)=8 ! top
isinside=.false.
if(xp(1).lt.minval(xcorner(1,:)) .and. xp(1).gt.maxval(xcorner(1,:)) .and.     &
   xp(2).lt.minval(xcorner(2,:)) .and. xp(2).gt.maxval(xcorner(2,:)) .and.     &
   xp(3).lt.minval(xcorner(3,:)) .and. xp(3).gt.maxval(xcorner(3,:)))then
  ! point is outside the element
  return
endif
do i_face=1,6
  x0=xcorner(:,fnode3(1,i_face))
  x1=xcorner(:,fnode3(2,i_face))
  x2=xcorner(:,fnode3(3,i_face))
  normvec=-get_normal(x0,x1,x2) ! negative makes inward normal
  
  ! vector to the point
  dxp=xp-x0
  normp=norm(dxp)
  if(normp.eq.tol)then
    ! xp is x0 hence xp is on the plane
    ! no need to check for all faces
    isinside=.true.
    return
  endif
  uvecp=dxp/normp
  cosine=dot_product(normvec,uvecp)
  if(cosine.gt.one.or.cosine.lt.-one)then
    print*,'ERROR: invalid cosine value!',cosine
    stop
  endif
  if(abs(cosine).lt.tol)then
    ! point is on the surface but due to numerical precision is outside.
    ! consider this as inside. this my happen due to meshing.
    cycle
  endif
  if(cosine.lt.zero)then
   ! point is outside the element
    return
  endif
enddo
isinside=.true.
return
end function IsPointInHexahedron
!===============================================================================

recursive function factorial(n) result(nfact)
implicit none
integer, intent(in) :: n
integer :: nfact
if(n > 0) then
  nfact = n * factorial(n-1)
  return
elseif (n==0)then
  nfact = 1
else
  write(*,*)'ERROR: undefined factorial!'
  stop
end if
end function factorial
!===============================================================================

! this function returns the determinant of a 1x1, 2x2 or 3x3
! jacobian matrix.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
function determinant(jac)result(det)
implicit none
real(kind=kreal),intent(in)::jac(:,:)
real(kind=kreal)::det
integer::it
it=ubound(jac,1)
select case(it)
case(1)
  det=1.0_kreal
case(2)
  det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
case(3)
  det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
  det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
  det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
case default
  write(*,*)'ERROR: wrong dimension for jacobian matrix!'
end select
return
end function determinant
!=======================================================

! this subroutine inverts a small square matrix onto itself.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine invert(matrix)
implicit none
real(kind=kreal),intent(in out)::matrix(:,:)
real(kind=kreal)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
integer::ndim,i,k
ndim=ubound(matrix,1)
if(ndim==2)then
  det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
  j11=matrix(1,1)
  matrix(1,1)=matrix(2,2)
  matrix(2,2)=j11
  matrix(1,2)=-matrix(1,2)
  matrix(2,1)=-matrix(2,1)
  matrix=matrix/det
else if(ndim==3)then
  det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
  det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
  det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
  j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
  j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
  j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
  j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
  j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
  j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
  j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
  j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
  j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
  matrix(1,1)=j11
  matrix(1,2)=j12
  matrix(1,3)=j13
  matrix(2,1)=j21
  matrix(2,2)=j22
  matrix(2,3)=j23
  matrix(3,1)=j31
  matrix(3,2)=j32
  matrix(3,3)=j33
  matrix=matrix/det
else
  do k=1,ndim
    con=matrix(k,k)
    matrix(k,k)=1.0_kreal
    matrix(k,:)=matrix(k,:)/con
    do i=1,ndim
      if(i/=k)then
        con=matrix(i,k)
        matrix(i,k)=0.0_kreal
        matrix(i,:)=matrix(i,:)-matrix(k,:)*con
      end if
    end do
  end do
end if
return
end subroutine invert
!=======================================================

! this subroutine forms the stress invariants in 2- or 3-d.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine stress_invariant(stress,sigm,dsbar,theta)
implicit none
real(kind=kreal),intent(in)::stress(:)
real(kind=kreal),intent(out),optional::sigm,dsbar,theta
real(kind=kreal)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,s1,s2,s3,s4,s5,s6,ds1,ds2,ds3,&
  d2,d3,sq3,zero=0.0_kreal,small=1.e-12_kreal,one=1.0_kreal,two=2.0_kreal,     &
  three=3.0_kreal,six=6.0_kreal,thpt5=13.5_kreal
integer::nst
nst=ubound(stress,1)
select case(nst)
case(4)
  sx=stress(1)
  sy=stress(2)
  txy=stress(3)
  sz=stress(4)
  sigm=(sx+sy+sz)/three
  dsbar=sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)/sqrt(two)
  if(dsbar<small)then
    theta=zero
  else
    dx=(two*sx-sy-sz)/three
    dy=(two*sy-sz-sx)/three
    dz=(two*sz-sx-sy)/three
    xj3=dx*dy*dz-dz*txy**2
    sine=-thpt5*xj3/dsbar**3
    if(sine>=one)sine=one
    if(sine<-one)sine=-one
    theta=asin(sine)/three
  end if
case(6)
  sq3=sqrt(three)
  s1=stress(1)
  s2=stress(2)
  s3=stress(3)
  s4=stress(4)
  s5=stress(5)
  s6=stress(6)
  sigm=(s1+s2+s3)/three
  d2=((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/six+s4*s4+s5*s5+s6*s6

  if(d2<small)d2=small ! special case of hydrostatic pressure or just at the tip

  ds1=s1-sigm
  ds2=s2-sigm
  ds3=s3-sigm
  d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+two*s4*s5*s6
  dsbar=sq3*sqrt(d2)
  if(dsbar<small)then
    theta=zero
  else
    sine=-three*sq3*d3/(two*sqrt(d2)**3)
    if(sine>=one)sine=one
    if(sine<-one)sine=-one
    theta=asin(sine)/three
  end if
case default
  write(*,*)"ERROR: wrong size for nst in invar!"
end select
return
end subroutine stress_invariant
!=======================================================

! quick sort of integer list                                                     
function iquick_sort(x,n) result(xnew)                                           
integer,intent(in) :: n ! size of the vector data x                              
integer, dimension(n) :: x ! data vector to sort                                 
integer :: temp                                                                  
integer :: i,j                                                                   
integer,dimension(n) :: xnew                                                     
                                                                                 
do i = 2, n                                                                      
  j = i - 1                                                                      
  temp = x(i)                                                                    
  do while (j>=1 .and. x(j)>temp)                                                
    x(j+1) = x(j)                                                                
    j = j - 1                                                                    
  end do                                                                         
  x(j+1) = temp                                                                  
end do                                                                           
xnew=x                                                                           
end function iquick_sort                                                         
!======================================================= 

! quick sort of real list
function rquick_sort(x,n) result(xnew)
integer,intent(in) :: n ! size of the vector data x
real(kind=kreal), dimension(n) :: x ! data vector to sort
real(kind=kreal) :: temp
integer :: i,j
real(kind=kreal),dimension(n) :: xnew

do i = 2, n
  j = i - 1
  temp = x(i)
  do while (j>=1 .and. x(j)>temp)
    x(j+1) = x(j)
    j = j - 1
  end do
  x(j+1) = temp
end do
xnew=x
end function rquick_sort
!=======================================================

! insertion sort of integer list
subroutine insertion_sort(x,n)
integer,intent(in) :: n ! size of the vector data x
real, intent(inout), dimension(n) :: x ! data vector to sort
real :: temp
integer :: i, j

do i = 2, n
  j = i - 1
  temp = x(i)
  do while (j>=1 .and. x(j)>temp)
    x(j+1) = x(j)
    j = j - 1
  end do
  x(j+1) = temp
end do
end subroutine insertion_sort
!=======================================================

! Author: Michel Olagnon
! orderpack 2.0
! source: http://www.fortran-2000.com/rank/
Subroutine i8_uniinv (XDONT, IGOEST)
! UNIINV = Merge-sort inverse ranking of an array, with removal of
! duplicate entries.
! this routine is similar to pure merge-sort ranking, but on
! the last pass, it sets indices in IGOEST to the rank
! of the value in the ordered set with duplicates removed.
! for performance reasons, the first 2 passes are taken
! out of the standard loop, and use dedicated coding.
implicit none
integer,parameter :: kint8=selected_int_kind(13)
integer(kind=kint8),intent(in)  :: XDONT(:)
integer(kind=kint8),intent(out) :: IGOEST(:)

integer(kind=kint8) :: XTST, XDONA, XDONB
integer(kind=kint8), dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
integer(kind=kint8) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
integer(kind=kint8) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
select case (NVAL)
case (:0)
  return
case (1)
  IGOEST (1) = 1
  return
case default
  continue
end select

! fill-in the index array, creating ordered couples
do IIND = 2, NVAL, 2
  if(XDONT(IIND-1) < XDONT(IIND)) then
    IRNGT (IIND-1) = IIND - 1
    IRNGT (IIND) = IIND
  else
    IRNGT (IIND-1) = IIND
    IRNGT (IIND) = IIND - 1
  endif
enddo
if(modulo(NVAL,2_kint8) /= 0) then
  IRNGT (NVAL) = NVAL
endif

! we will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into     C   -   C   - ...
LMTNA = 2
LMTNC = 4

! first iteration. The length of the ordered subsets goes from 2 to 4
do
  if (NVAL <= 4) Exit
  ! loop on merges of A and B into C
  do IWRKD = 0, NVAL - 1, 4
    if ((IWRKD+4) > NVAL) then
      if ((IWRKD+2) >= NVAL) Exit
      !   1 2 3
      if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !   1 3 2
      if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
        IRNG2 = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNG2
        !   3 1 2
      else
        IRNG1 = IRNGT (IWRKD+1)
        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNG1
      endif
      exit
    endif
    !   1 2 3 4
    if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !   1 3 x x
    if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
      IRNG2 = IRNGT (IWRKD+2)
      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
      if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
        !   1 3 2 4
        IRNGT (IWRKD+3) = IRNG2
      else
        !   1 3 4 2
        IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
        IRNGT (IWRKD+4) = IRNG2
      endif
      !   3 x x x
    else
      IRNG1 = IRNGT (IWRKD+1)
      IRNG2 = IRNGT (IWRKD+2)
      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
      if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
        IRNGT (IWRKD+2) = IRNG1
        if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
          !   3 1 2 4
          IRNGT (IWRKD+3) = IRNG2
        else
          !   3 1 4 2
          IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
          IRNGT (IWRKD+4) = IRNG2
        endif
      else
        !   3 4 1 2
        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
        IRNGT (IWRKD+3) = IRNG1
        IRNGT (IWRKD+4) = IRNG2
      endif
    endif
  enddo

! the Cs become As and Bs
  LMTNA = 4
  Exit
enddo

! iteration loop. Each time, the length of the ordered subsets
! is doubled.
do
  if (2*LMTNA >= NVAL) Exit
  IWRKF = 0
  LMTNC = 2 * LMTNC

  ! loop on merges of A and B into C
  do
    IWRK = IWRKF
    IWRKD = IWRKF + 1
    JINDA = IWRKF + LMTNA
    IWRKF = IWRKF + LMTNC
    if (IWRKF >= NVAL) then
      if (JINDA >= NVAL) Exit
      IWRKF = NVAL
    endif
    IINDA = 1
    IINDB = JINDA + 1

    ! ONE steps in the C subset, that we create in the final rank array
    ! make a copy of the rank array for the iteration
    JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
    XDONA = XDONT (JWRKT(IINDA))
    XDONB = XDONT (IRNGT(IINDB))
    do
      IWRK = IWRK + 1
      ! we still have unprocessed values in both A and B
      if (XDONA > XDONB) then
        IRNGT (IWRK) = IRNGT (IINDB)
        IINDB = IINDB + 1
        if (IINDB > IWRKF) then
          ! only A still with unprocessed values
          IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
          Exit
        endif
        XDONB = XDONT (IRNGT(IINDB))
      else
        IRNGT (IWRK) = JWRKT (IINDA)
        IINDA = IINDA + 1
        if (IINDA > LMTNA) Exit! Only B still with unprocessed values
        XDONA = XDONT (JWRKT(IINDA))
      endif

    enddo
  enddo

  ! the Cs become As and Bs
  LMTNA = 2 * LMTNA
enddo

! last merge of A and B into C, with removal of duplicates.
IINDA = 1
IINDB = LMTNA + 1
NUNI = 0

! ONE steps in the C subset, that we create in the final rank array
JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
if (IINDB <= NVAL) then
  XTST = i8_nearless(Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
else
  XTST = i8_nearless(XDONT(JWRKT(1)))
endif

do IWRK = 1, NVAL
  ! we still have unprocessed values in both A and B
  if (IINDA <= LMTNA) then
    if (IINDB <= NVAL) then
      if (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) then
        IRNG = IRNGT (IINDB)
        IINDB = IINDB + 1
      else
        IRNG = JWRKT (IINDA)
        IINDA = IINDA + 1
      endif
    else
      ! only A still with unprocessed values
      IRNG = JWRKT (IINDA)
      IINDA = IINDA + 1
    endif
  else
    ! only B still with unprocessed values
    IRNG = IRNGT (IWRK)
  endif
  if (XDONT(IRNG) > XTST) then
    XTST = XDONT (IRNG)
    NUNI = NUNI + 1
  endif
  IGOEST (IRNG) = NUNI
enddo
return
end subroutine i8_uniinv
!=======================================================

function i8_nearless (XVAL) result (I_nl)
! nearest value less than given value
implicit none
integer,parameter :: kint8=selected_int_kind(13)
integer(kind=kint8),intent(in) :: XVAL
integer(kind=kint8) :: I_nl
I_nl = XVAL - 1
return
end function i8_nearless
!=======================================================

Subroutine i_uniinv (XDONT, IGOEST)
! UNIINV = Merge-sort inverse ranking of an array, with removal of
! duplicate entries.
! this routine is similar to pure merge-sort ranking, but on
! the last pass, it sets indices in IGOEST to the rank
! of the value in the ordered set with duplicates removed.
! for performance reasons, the first 2 passes are taken
! out of the standard loop, and use dedicated coding.
implicit none
integer,intent(in)  :: XDONT(:)
integer,intent(out) :: IGOEST(:)

integer :: XTST, XDONA, XDONB
integer, dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
select case (NVAL)
case (:0)
  return
case (1)
  IGOEST (1) = 1
  return
case default
  continue
end select

! fill-in the index array, creating ordered couples
do IIND = 2, NVAL, 2
  if(XDONT(IIND-1) < XDONT(IIND)) then
    IRNGT (IIND-1) = IIND - 1
    IRNGT (IIND) = IIND
  else
    IRNGT (IIND-1) = IIND
    IRNGT (IIND) = IIND - 1
  endif
enddo
if(modulo(NVAL,2) /= 0) then
  IRNGT (NVAL) = NVAL
endif

! we will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into     C   -   C   - ...
LMTNA = 2
LMTNC = 4

! first iteration. The length of the ordered subsets goes from 2 to 4
do
  if (NVAL <= 4) Exit
  ! loop on merges of A and B into C
  do IWRKD = 0, NVAL - 1, 4
    if ((IWRKD+4) > NVAL) then
      if ((IWRKD+2) >= NVAL) Exit
      !   1 2 3
      if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !   1 3 2
      if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
        IRNG2 = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNG2
        !   3 1 2
      else
        IRNG1 = IRNGT (IWRKD+1)
        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNG1
      endif
      exit
    endif
    !   1 2 3 4
    if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !   1 3 x x
    if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
      IRNG2 = IRNGT (IWRKD+2)
      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
      if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
        !   1 3 2 4
        IRNGT (IWRKD+3) = IRNG2
      else
        !   1 3 4 2
        IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
        IRNGT (IWRKD+4) = IRNG2
      endif
      !   3 x x x
    else
      IRNG1 = IRNGT (IWRKD+1)
      IRNG2 = IRNGT (IWRKD+2)
      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
      if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
        IRNGT (IWRKD+2) = IRNG1
        if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
          !   3 1 2 4
          IRNGT (IWRKD+3) = IRNG2
        else
          !   3 1 4 2
          IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
          IRNGT (IWRKD+4) = IRNG2
        endif
      else
        !   3 4 1 2
        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
        IRNGT (IWRKD+3) = IRNG1
        IRNGT (IWRKD+4) = IRNG2
      endif
    endif
  enddo

! the Cs become As and Bs
  LMTNA = 4
  Exit
enddo

! iteration loop. Each time, the length of the ordered subsets
! is doubled.
do
  if (2*LMTNA >= NVAL) Exit
  IWRKF = 0
  LMTNC = 2 * LMTNC

  ! loop on merges of A and B into C
  do
    IWRK = IWRKF
    IWRKD = IWRKF + 1
    JINDA = IWRKF + LMTNA
    IWRKF = IWRKF + LMTNC
    if (IWRKF >= NVAL) then
      if (JINDA >= NVAL) Exit
      IWRKF = NVAL
    endif
    IINDA = 1
    IINDB = JINDA + 1

    ! ONE steps in the C subset, that we create in the final rank array
    ! make a copy of the rank array for the iteration
    JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
    XDONA = XDONT (JWRKT(IINDA))
    XDONB = XDONT (IRNGT(IINDB))
    do
      IWRK = IWRK + 1
      ! we still have unprocessed values in both A and B
      if (XDONA > XDONB) then
        IRNGT (IWRK) = IRNGT (IINDB)
        IINDB = IINDB + 1
        if (IINDB > IWRKF) then
          ! only A still with unprocessed values
          IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
          Exit
        endif
        XDONB = XDONT (IRNGT(IINDB))
      else
        IRNGT (IWRK) = JWRKT (IINDA)
        IINDA = IINDA + 1
        if (IINDA > LMTNA) Exit! Only B still with unprocessed values
        XDONA = XDONT (JWRKT(IINDA))
      endif

    enddo
  enddo

  ! the Cs become As and Bs
  LMTNA = 2 * LMTNA
enddo

! last merge of A and B into C, with removal of duplicates.
IINDA = 1
IINDB = LMTNA + 1
NUNI = 0

! ONE steps in the C subset, that we create in the final rank array
JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
if (IINDB <= NVAL) then
  XTST = i_nearless(Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
else
  XTST = i_nearless(XDONT(JWRKT(1)))
endif

do IWRK = 1, NVAL
  ! we still have unprocessed values in both A and B
  if (IINDA <= LMTNA) then
    if (IINDB <= NVAL) then
      if (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) then
        IRNG = IRNGT (IINDB)
        IINDB = IINDB + 1
      else
        IRNG = JWRKT (IINDA)
        IINDA = IINDA + 1
      endif
    else
      ! only A still with unprocessed values
      IRNG = JWRKT (IINDA)
      IINDA = IINDA + 1
    endif
  else
    ! only B still with unprocessed values
    IRNG = IRNGT (IWRK)
  endif
  if (XDONT(IRNG) > XTST) then
    XTST = XDONT (IRNG)
    NUNI = NUNI + 1
  endif
  IGOEST (IRNG) = NUNI
enddo
return
end subroutine i_uniinv
!=======================================================

function i_nearless (XVAL) result (I_nl)
! nearest value less than given value
implicit none
integer,intent(in) :: XVAL
integer :: I_nl
I_nl = XVAL - 1
return
end function i_nearless
!=======================================================

function upind_i2j(nup,n,inode,jnode) result(iup)
implicit none
integer,intent(in) :: nup ! number of upper trianglular elements
integer,intent(in) :: n ! number of rows or columns of a full square matrix 
integer,intent(in) :: inode,jnode ! row, column 
integer :: iup ! index for upper triangular element

integer :: i,j

iup=-1
if(inode.eq.jnode)then
  iup=0
else
  if(inode.gt.jnode)then
    i=jnode
    j=inode
  else
    i=inode
    j=jnode
  endif
  !print*,i_node,j_node
  iup=nup - (n-i)*(n-i+1)/2 +j - i
endif
return
end function upind_i2j
!=======================================================

function upinds_i2js(nup,n,inode,nj,jnodes) result(iups)
implicit none
integer,intent(in) :: nup ! number of upper trianglular elements
integer,intent(in) :: n ! number of rows or columns of a full square matrix 
integer,intent(in) :: inode,nj,jnodes(nj) ! row, number of columns, columns 
integer :: iups(nj) ! indices for upper triangular elements

integer :: i,j,i_j,jnode

iups=-1
do i_j=1,nj
  jnode=jnodes(i_j)
  iups(i_j)=upind_i2j(nup,n,inode,jnode)
enddo
return
end function upinds_i2js
!=======================================================

end module math_library
