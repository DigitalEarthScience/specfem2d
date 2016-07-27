!AUTHORS:
!Hom Nath Gharti
!Stefano Zhampini
!REFERENCE:
!PETSC documentation
! All routines below are empty. these have to be modified for petsc serial solver
! if available. Serial version curerntly use builtin solver.
! -----------------------------------------------------------------------
module solver_petsc
use ksp_constants
use global
implicit none
contains
!=======================================================
subroutine petsc_initialize()
implicit none

end subroutine petsc_initialize
!=======================================================

subroutine petsc_create_vector()
implicit none

end subroutine petsc_create_vector
!=======================================================

subroutine petsc_matrix_preallocate_size()
implicit none

end subroutine petsc_matrix_preallocate_size
!=======================================================

subroutine petsc_create_matrix()
implicit none

end subroutine petsc_create_matrix
!=======================================================

subroutine petsc_create_solver()
implicit none

end subroutine petsc_create_solver
!=======================================================

subroutine petsc_set_ksp_operator
implicit none

end subroutine petsc_set_ksp_operator
!=======================================================

subroutine petsc_set_ksp_operatorconv
implicit none

end subroutine petsc_set_ksp_operatorconv
!=======================================================

subroutine petsc_set_stiffness_matrix(storekmat)
use ieee_arithmetic
implicit none
real(kind=kreal),intent(in) :: storekmat(:,:,:)

end subroutine petsc_set_stiffness_matrix
!=======================================================

subroutine petsc_set_stiffness_matrixconv(nup,nvalency,interpfgll,gdof_elmt,nodalnu0,nodalgauss)
use math_constants,only:gtol
use math_library,only:upinds_i2js
!use math_library_mpi,only:sumscal
use global,only:storederiv,storejw
use ieee_arithmetic
implicit none
integer,intent(in) :: nup
integer,intent(in) :: nvalency(:)
real(kind=kreal),intent(in) :: interpfgll(ngll,ngll)
integer,intent(in) :: gdof_elmt(:,:)
real(kind=kreal),intent(in) :: nodalnu0(:,:),nodalgauss(0:)

end subroutine petsc_set_stiffness_matrixconv
!=======================================================

subroutine petsc_set_vector(rload)
!use global,only:l2gdof,nelmt,NEDOF
use ieee_arithmetic
implicit none
!PetscScalar,intent(in) :: rload(0:)
real(kind=kreal),intent(in) :: rload(0:)

end subroutine petsc_set_vector
!=======================================================

subroutine petsc_set_initialguess(rload)
!use global,only:l2gdof,NEDOF
implicit none
!PetscScalar,intent(in) :: rload(0:)
real(kind=kreal),intent(in) :: rload(0:)

end subroutine petsc_set_initialguess
!=======================================================

subroutine petsc_solve(sdata,cg_iter)
implicit none
!PetscScalar sdata(:)
!PetscInt    cg_iter
real(kind=kreal),intent(in) :: sdata(:)
integer :: cg_iter

end subroutine petsc_solve
!=======================================================

subroutine petsc_load()
implicit none

end subroutine petsc_load
!=======================================================

subroutine petsc_save()
implicit none

end subroutine petsc_save
!=======================================================

subroutine petsc_finalize()
implicit none

end subroutine petsc_finalize
!=======================================================

end module solver_petsc
!=======================================================
