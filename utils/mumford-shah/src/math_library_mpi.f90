! MPI math library
module math_library_mpi
use math_constants
use set_precision_mpi

private :: iminscal,fminscal
private :: iminvec,fminvec
private :: isumscal,fsumscal
private :: imaxscal,fmaxscal
private :: imaxvec,fmaxvec

! global sum of a scalar in all processors
interface sumscal
  module procedure isumscal
  module procedure fsumscal
end interface

! global maximum of a scalar in all processors
interface minscal
  module procedure iminscal
  module procedure fminscal
end interface

! global maximum of a scalar in all processors
interface maxscal
  module procedure imaxscal
  module procedure fmaxscal
end interface

! global maximum of a vector in all processors
interface maxvec
  module procedure imaxvec
  module procedure fmaxvec
end interface

! global minimum of a scalar in all processors
interface minvec
  module procedure iminvec
  module procedure fminvec
end interface
contains
!=======================================================
!=======================================================

function iminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmin
integer :: scal_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmin_mpi,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function iminscal
!=======================================================

function fminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmin
real(kind=kreal) :: scal_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmin_mpi,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function fminscal
!=======================================================

function imaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmax
integer :: scal_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr
scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmax_mpi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function imaxscal
!=======================================================

function fmaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmax
real(kind=kreal) :: scal_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmax_mpi,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function fmaxscal
!=======================================================

function imaxvec(vec) result(gmax)
implicit none
integer,intent(in)::vec(:)
integer :: gmax ! global
integer :: lmax_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

lmax_mpi(1)=maxval(vec)
call MPI_ALLREDUCE(lmax_mpi,gmax_mpi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function imaxvec
!=======================================================

function fmaxvec(vec) result(gmax)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: gmax ! global
real(kind=kreal) :: lmax_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

lmax_mpi(1)=maxval(vec)
call MPI_ALLREDUCE(lmax_mpi,gmax_mpi,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function fmaxvec
!=======================================================

function iminvec(vec) result(gmin)
implicit none
integer,intent(in)::vec(:)
integer :: gmin ! global
integer :: lmin_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

lmin_mpi(1)=minval(vec)
call MPI_ALLREDUCE(lmin_mpi,gmin_mpi,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function iminvec
!=======================================================

function fminvec(vec) result(gmin)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: gmin ! global
real(kind=kreal) :: lmin_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

lmin_mpi(1)=minval(vec)
call MPI_ALLREDUCE(lmin_mpi,gmin_mpi,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function fminvec
!=======================================================

function isumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gsum ! global
integer :: scal_mpi(1),gsum_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gsum_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
gsum=gsum_mpi(1)
return
end function isumscal
!=======================================================

function fsumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gsum ! global
real(kind=kreal) :: scal_mpi(1),gsum_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gsum_mpi,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
gsum=gsum_mpi(1)
return
end function fsumscal
!=======================================================

function dot_product_par(vec1,vec2) result(gdot)
!
! this finds global dot product of TWO vectors across the processors
!
implicit none
real(kind=kreal),intent(in)::vec1(:),vec2(:)
real(kind=kreal) :: gdot ! global
real(kind=kreal) :: ldot_mpi(1),gdot_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

! find local dot
ldot_mpi(1)=dot_product(vec1,vec2)
call MPI_ALLREDUCE(ldot_mpi,gdot_mpi,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
gdot=gdot_mpi(1)
return
end function dot_product_par
!=======================================================

end module math_library_mpi
