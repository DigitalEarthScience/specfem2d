! this module contains preprocessing library routines
! REVISION:
!   HNG, Jul 07,2011
module preprocess
contains
subroutine allocate_store_arrays
use global
use math_constants
implicit none
allocate(storederiv(ndim,ngll,ngll,nelmt))
allocate(storejw(ngll,nelmt))
storederiv=zero
storejw=zero
end subroutine allocate_store_arrays
!===============================================================================

subroutine precompute_quadrature(gnod,lagrange_gll)
use global
use math_constants
use math_library
use shape_library
use gll_library
implicit none
integer,intent(in) :: gnod(4)
real(kind=kreal),intent(out) :: lagrange_gll(ngll,ngll)

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal

integer :: num(nenod)
integer :: i,i_elmt
real(kind=kreal) :: coord(ngnode,ndim),jac(ndim,ndim)
real(kind=kreal) :: dshape_quad4(ndim,ngnode,ngll)
real(kind=kreal) :: detjac
real(kind=kreal) :: xigll(ngllx),wxgll(ngllx),etagll(nglly),wygll(nglly)
! compute gauss-lobatto-legendre quadrature information
real(kind=kreal) :: gll_weights(ngll),gll_points(ndim,ngll),                   &
dlagrange_gll(ndim,ngll,ngll)

! get gll points and weights
call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)

! get derivatives of shape functions for 4-noded hex
call dshape_function_quad4(ngnode,ngllx,nglly,xigll,etagll,dshape_quad4)
call gll_quadrature2d(ndim,ngllx,nglly,ngll,gll_points,gll_weights,        &
lagrange_gll,dlagrange_gll)

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))   
  do i=1,ngll
    !interpf(:,1)=lagrange_gll(i,:)

    ! compute Jacobian
    jac=matmul(dshape_quad4(:,:,i),coord)  
    detjac=determinant(jac)
    if(detjac.le.ZERO)stop 'negative or zero jacobian!'
    call invert(jac)
    
    storederiv(:,:,i,i_elmt)=matmul(jac,dlagrange_gll(:,i,:))

    storejw(i,i_elmt)=detjac*gll_weights(i)
  end do ! i=1,ngll
end do 

end subroutine precompute_quadrature
!===============================================================================

subroutine stiffness_load_edgefield(gdof_elmt,interpfgll,nodalm0,storekmat,dprecon,bodyload)
use global
use math_constants,only:zero
use math_library,only:norm
implicit none
integer,intent(in) :: gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
real(kind=kreal),intent(in) :: nodalm0(1,nnode)
real(kind=kreal),intent(out) :: storekmat(nedof,nedof,nelmt)
real(kind=kreal),intent(out) :: dprecon(0:neq),bodyload(0:neq)

integer :: egdof(nedof),num(ngll)
real(kind=kreal) :: dinterpf(ndim,ngll),interpf(ngll,1)
real(kind=kreal) :: kmat(nedof,nedof),kmat_gradnu(nedof,nedof),kmat_nu(nedof,nedof)
real(kind=kreal) :: eload(nedof),gradm(ndim,1),mgllmat(nedof,1)
real(kind=kreal) :: gradmsq
integer :: i_elmt,i

! compute stiffness, body load, and precondioner
storekmat=zero
dprecon=zero
bodyload=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  egdof=gdof_elmt(:,i_elmt) 
  mgllmat(:,1)=nodalm0(1,num)
  kmat=zero; eload=zero
  do i=1,ngll
    interpf(:,1)=interpfgll(i,:)
    dinterpf=storederiv(:,:,i,i_elmt)
    gradm=matmul(dinterpf,mgllmat)
    gradmsq=norm(gradm(:,1))
  
    kmat_gradnu = eps_img*matmul(transpose(dinterpf),dinterpf)
    kmat_nu     = (0.25_kreal/eps_img+gam_img*gradmsq/eta_img)*matmul(interpf,transpose(interpf))
    kmat = kmat + (kmat_gradnu+kmat_nu)*storejw(i,i_elmt)
    eload=eload+(0.25_kreal/eps_img)*interpfgll(i,:)*storejw(i,i_elmt)
  end do ! i=1,ngll
  storekmat(:,:,i_elmt)=kmat
  do i=1,nedof
    dprecon(egdof(i))=dprecon(egdof(i))+kmat(i,i)
  end do
  bodyload(egdof)=bodyload(egdof)+eload
end do 
return
end subroutine stiffness_load_edgefield
!===============================================================================

subroutine stiffness_load_model(gdof_elmt,interpfgll,nodalg,nodalnu0,storekmat,dprecon,bodyload)
use global
use math_constants,only:zero
use math_library,only:determinant,norm
use math_constants,only:half,pi
implicit none
integer,intent(in) :: gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
real(kind=kreal),intent(in) :: nodalg(1,nnode),nodalnu0(1,nnode)
real(kind=kreal),intent(out) :: storekmat(nedof,nedof,nelmt)
real(kind=kreal),intent(out) :: dprecon(0:neq),bodyload(0:neq)

integer :: egdof(nedof),num(ngll),iups(ngll)
real(kind=kreal) :: dinterpf(ndim,ngll),interpf(ngll,1)
real(kind=kreal) :: kmat(nedof,nedof),kmat_gradm(nedof,nedof),kmat_m(nedof,nedof)
real(kind=kreal) :: eload(nedof),gradm(ndim,1),ggll(nedof),nugll(nedof),nugllmat(nedof,1)
real(kind=kreal) :: gradnu2(ndim,1)
real(kind=kreal) :: gradmsq
real(kind=kreal) :: cgfac
integer :: i_elmt,i

! compute stiffness, body load, and precondioner
storekmat=zero
dprecon=zero
bodyload=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  
  egdof=gdof_elmt(:,i_elmt) 
  nugll=nodalnu0(1,num)
  nugllmat(:,1)=nodalnu0(1,num)
  ggll=nodalg(1,num)
  kmat=zero; kmat_gradm=zero; kmat_m=zero
  eload=zero
  do i=1,ngll
    interpf(:,1)=interpfgll(i,:)
    dinterpf=storederiv(:,:,i,i_elmt)
     
    kmat_gradm = gam_img*nugll(i)*nugll(i)*matmul(transpose(dinterpf),dinterpf)
    kmat_m = matmul(interpf,transpose(interpf))
    kmat = kmat + (kmat_gradm+kmat_m)*storejw(i,i_elmt)
    eload=eload+ggll(i)*interpfgll(i,:)*storejw(i,i_elmt)
  end do ! i=1,ngll
  storekmat(:,:,i_elmt)=kmat
  do i=1,nedof
    dprecon(egdof(i))=dprecon(egdof(i))+kmat(i,i)
  end do
  bodyload(egdof)=bodyload(egdof)+eload
end do 
return
end subroutine stiffness_load_model
!===============================================================================

subroutine compute_bodyload_model(gdof_elmt,interpfgll,nodalg,bodyload)
use global
use math_constants,only:zero
implicit none
integer,intent(in) :: gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
real(kind=kreal),intent(in) :: nodalg(1,nnode)
real(kind=kreal),intent(out) :: bodyload(0:neq)

integer :: i_elmt,i
integer :: egdof(nedof),num(ngll)
real(kind=kreal) :: eload(nedof),ggll(nedof)
real(kind=kreal) :: interpf(ngll,1)

! compute stiffness, body load, and precondioner
bodyload=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  egdof=gdof_elmt(:,i_elmt) 
  ggll=nodalg(1,num)
  eload=zero
  do i=1,ngll
    interpf(:,1)=interpfgll(i,:)
     
    eload=eload+ggll(i)*interpfgll(i,:)*storejw(i,i_elmt)
  end do ! i=1,ngll
  bodyload(egdof)=bodyload(egdof)+eload
end do 
return
end subroutine compute_bodyload_model
!===============================================================================

subroutine stiffness_load_modelconv(gdof_elmt,interpfgll,nodalg,nodalnu0,gsigma,nup,nodalgauss,storekmat,dprecon,bodyload,xmterm)
use global
use math_constants,only:zero
use math_library,only:determinant,norm,upinds_i2js
use math_constants,only:half,pi
implicit none
integer,intent(in) :: gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
real(kind=kreal),intent(in) :: nodalg(1,nnode),nodalnu0(1,nnode)
real(kind=kreal),intent(in) :: gsigma
integer,intent(in) :: nup
real(kind=kreal),intent(in) :: nodalgauss(0:) ! nodal data to be convolved
real(kind=kreal),intent(out) :: storekmat(nedof,nedof,nelmt)
real(kind=kreal),intent(out) :: dprecon(0:neq),bodyload(0:neq)
logical,intent(in) :: xmterm ! if this is .true. kmat_m term will be discarded 
! because that term goes to the RHS.

integer :: egdof(nedof),num(ngll),iups(ngll)
real(kind=kreal) :: dinterpf(ndim,ngll),interpf(ngll,1),cginterpf(ngll,1)
real(kind=kreal) :: kmat(nedof,nedof),kmat_gradm1(nedof,nedof),kmat_gradm2(nedof,nedof),kmat_m(nedof,nedof)
real(kind=kreal) :: eload(nedof),gradm(ndim,1),ggll(nedof),nugll(nedof),nugllmat(nedof,1)
real(kind=kreal) :: gradnu2(ndim,1)
real(kind=kreal) :: gradmsq
real(kind=kreal) :: cgfac
integer :: i_elmt,i

cgfac=half/(pi*2.0_kreal*gsigma*gsigma)
! compute stiffness, body load, and precondioner
storekmat=zero
dprecon=zero
bodyload=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  
  egdof=gdof_elmt(:,i_elmt) 
  nugll=nodalnu0(1,num)
  nugllmat(:,1)=nodalnu0(1,num)
  ggll=nodalg(1,num)
  kmat=zero; kmat_gradm1=zero; kmat_gradm2=zero; kmat_m=zero
  eload=zero
  do i=1,ngll
    iups=upinds_i2js(nup,nnode,num(i),ngll,num)
    interpf(:,1)=interpfgll(i,:)
    dinterpf=storederiv(:,:,i,i_elmt)
     
    gradnu2=matmul(dinterpf,nugllmat*nugllmat)
    
    cginterpf(:,1)=nodalgauss(iups)*storejw(:,i_elmt)
    !kmat_gradm1 = gam_img*matmul(interpf,matmul(transpose(gradnu2),dinterpf))
    kmat_gradm2 = gam_img*nugll(i)*nugll(i)*matmul(transpose(dinterpf),dinterpf)
    if(.not.xmterm)kmat_m = matmul(interpf,transpose(cginterpf))
    kmat = kmat + (kmat_gradm2+kmat_m)*storejw(i,i_elmt)
    eload=eload+ggll(i)*interpfgll(i,:)*storejw(i,i_elmt)
  end do ! i=1,ngll
  storekmat(:,:,i_elmt)=kmat
  do i=1,nedof
    dprecon(egdof(i))=dprecon(egdof(i))+kmat(i,i)
  end do
  bodyload(egdof)=bodyload(egdof)+eload
end do 
return
end subroutine stiffness_load_modelconv
!===============================================================================

end module preprocess

