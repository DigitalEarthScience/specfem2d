! this module contains post processing library routines
module postprocess
use set_precision
contains

subroutine compute_image_energy(nodalg,nodalm,nodalnu,energy)
use global
use math_constants,only:one,zero
use math_library,only:norm
implicit none
real(kind=kreal),intent(in) :: nodalg(1,nnode),nodalm(1,nnode),nodalnu(1,nnode)
real(kind=kreal),intent(out) :: energy
integer :: i_elmt,i
integer :: num(ngll)
real(kind=kreal) :: mgllmat(nedof,1),nugllmat(nedof,1)
real(kind=kreal) :: dmgsq(ngll),gradm(ndim,1),gradmsq,gradnu(ndim,1),gradnusq
real(kind=kreal) :: dinterpf(ndim,ngll)
! compute energy
energy=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  mgllmat(:,1)=nodalm(1,num)
  nugllmat(:,1)=nodalnu(1,num)
  dmgsq(:)=(nodalm(1,num)+nodalg(1,num))**2
  do i=1,ngll
    dinterpf=storederiv(:,:,i,i_elmt)

    gradm=matmul(dinterpf,mgllmat)
    gradmsq=norm(gradm(:,1))
    
    gradnu=matmul(dinterpf,nugllmat)
    gradnusq=norm(gradnu(:,1))
 
    energy=energy+( dmgsq(i)+                                                &
                    gam_img*nugllmat(i,1)**2*gradnusq+                       &
                    eta_img*(eps_img*gradnusq+0.25_kreal*(nugllmat(i,1)-one)**2/eps_img)) &
                  *storejw(i,i_elmt)
  end do ! i=1,ngll
end do 
end subroutine compute_image_energy
!===============================================================================

! this routine save data to files Ensight Gold format
! TODO: make it optional
subroutine save_nodal_data(ptail,format_str,istep,nnode,nodalv,vtag,vext,isstep)
use global,only:nndof,ngll,nst,out_path,file_head,savedata
use math_constants
use visual
implicit none
character(len=20),intent(in) :: format_str,ptail
integer,intent(in) :: istep,nnode
real(kind=kreal),intent(in) :: nodalv(nndof,nnode)
character(len=*),intent(in) :: vext
character(len=*),intent(in) :: vtag
logical,optional,intent(in) :: isstep

integer :: npart
character(len=250) :: out_fname
character(len=80) :: destag

! write displacement vector
if(present(isstep).and. .not.isstep)then
  write(out_fname,'(a)')trim(out_path)//trim(file_head)//trim(ptail)//'.'//trim(vext)
else
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',istep,trim(ptail)//'.'//trim(vext)
endif
npart=1;
destag=trim(vtag)
!print*,trim(out_fname)
call write_ensight_pernodeVECAS(out_fname,destag,npart,1,nnode,real(nodalv))

end subroutine save_nodal_data

!===============================================================================
! Jeroen's note on Mumford-Shah eq (2)
subroutine save_specfem_output(ismpi,gnod,format_str,gdof_elmt,interpfgll,nodalnu,nodalm,nvalency)
use global
use math_constants,only:ONE,TWO,HALF,zero
use math_library,only:norm
use string_library,only:parse_file
use blur,only:convolve_with_gaussian
implicit none
logical,intent(in) :: ismpi ! .true.=parallel version, .false.=serial version
integer,intent(in) :: gnod(4)!geometrical nodes (corner nodes) per element
character(len=20),intent(in) :: format_str
integer,intent(in) :: gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
real(kind=kreal),intent(in) :: nodalnu(1,nnode)
real(kind=kreal),intent(in) :: nodalm(1,nnode)
integer,intent(in) :: nvalency(nnode)

integer :: num(ngll)
real(kind=kreal) :: TwoGam,TwoEtaEps,HalfEtaOverEps
real(kind=kreal) :: dinterpf(ndim,ngll),interpf(ngll,1)
real(kind=kreal) :: gradmsq(ngll),gradm(ndim,1),nusqgradm(ndim,ngll),mgllmat(ngll,1)
real(kind=kreal) :: divgradnu,gradnu(ndim,1),gradnugll(ndim,ngll),nugllmat(ngll,1)
real(kind=kreal),allocatable :: dmterm(:,:),dnuterm(:,:)
real(kind=kreal),allocatable :: nodal_dmterm(:),nodal_dnuterm(:)
real(kind=kreal),allocatable :: nodal_dmtermc(:),nodal_dnutermc(:)
integer :: i_elmt,i

character(len=250) :: matfile_head
character(len=150) :: path
character(len=20) :: ext

character(len=250) :: errtag ! error message
integer :: errcode

errtag=""; errcode=-1

! factors
TwoGam=two*gam_img
TwoEtaEps=two*eta_img*eps_img
HalfEtaOverEps=half*eta_img/eps_img
! compute and write \nabla.(\nu^2 \nabla m)
allocate(dmterm(ngll,nelmt),dnuterm(ngll,nelmt))
dmterm=zero
dnuterm=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  
  mgllmat(:,1)=nodalm(1,num)
  
  nugllmat(:,1)=nodalnu(1,num)

  ! compute and store \nu^2 \grad m on gll points
  do i=1,ngll
    interpf(:,1)=interpfgll(i,:)
    dinterpf=storederiv(:,:,i,i_elmt)
    
    gradm=matmul(dinterpf,mgllmat)
    gradmsq(i)=dot_product(gradm(:,1),gradm(:,1))
    nusqgradm(:,i)=nugllmat(i,1)*nugllmat(i,1)*gradm(:,1)
    
    gradnu=matmul(dinterpf,nugllmat)
    gradnugll(:,i)=gradnu(:,1)
    
  end do ! i=1,ngll
  ! compute and store \nabla.\nu^2 \grad m on gll points
  do i=1,ngll
    dinterpf=storederiv(:,:,i,i_elmt)
    
    dmterm(i,i_elmt) = -TwoGam*(dot_product(dinterpf(1,:),gradnu(1,:)) +       &
                                dot_product(dinterpf(2,:),gradnu(2,:)))
    
    divgradnu = dot_product(dinterpf(1,:),gradnu(1,:)) +                       &
                dot_product(dinterpf(2,:),gradnu(2,:))
    dnuterm(i,i_elmt) = TwoGam*nugllmat(i,1)*gradmsq(i) -                              & 
                        TwoEtaEps*divgradnu +                                  &
                        HalfEtaOverEps*(nugllmat(i,1)-one)
    
  end do ! i=1,ngll
end do
!print*,maxval(abs(dnuterm)),maxval(abs(dmterm)) 
call parse_file(matfile,path,matfile_head,ext)
open(11,file=trim(out_path)//trim(matfile_head)//'_dnu.'//trim(ext),            &
access='stream',form='unformatted',action='write',status='replace')
write(11)dnuterm
close(11)
allocate(nodal_dnuterm(nnode),nodal_dnutermc(nnode))
nodal_dnuterm=zero

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  nodal_dnuterm(num)=nodal_dnuterm(num)+dnuterm(:,i_elmt)
enddo
nodal_dnuterm=nodal_dnuterm/nvalency
!call convolve_with_gaussian(ismpi,gnod,2.0_kreal,nodal_dnuterm,nodal_dnutermc,errcode,errtag)
call save_nodal_data(ptail,format_str,0,nnode,nodal_dnuterm,'dnu','dnu',.false.)
deallocate(nodal_dnuterm,nodal_dnutermc)

open(11,file=trim(out_path)//trim(matfile_head)//'_dm.'//trim(ext),             &
access='stream',form='unformatted',action='write',status='replace')
write(11)dmterm
close(11)
allocate(nodal_dmterm(nnode),nodal_dmtermc(nnode))
nodal_dmterm=zero

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  nodal_dmterm(num)=nodal_dmterm(num)+dmterm(:,i_elmt)
enddo
nodal_dmterm=nodal_dmterm/nvalency
!call convolve_with_gaussian(ismpi,gnod,2.0_kreal,nodal_dmterm,nodal_dmtermc,errcode,errtag)
call save_nodal_data(ptail,format_str,0,nnode,nodal_dmterm,'dm','dm',.false.)
deallocate(nodal_dmterm,nodal_dmtermc)
return
end subroutine save_specfem_output
!===============================================================================
end module postprocess
