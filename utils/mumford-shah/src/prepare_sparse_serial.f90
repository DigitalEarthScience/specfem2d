! this module contains an empty routine for equivalent serial version
module sparse_serial
contains
subroutine prepare_sparse(nup)!,nodal_gauss)
use set_precision
implicit none
integer,intent(in) :: nup
!real(kind=kreal),intent(in) :: nodal_gauss(:)
end subroutine prepare_sparse
end module sparse_serial
