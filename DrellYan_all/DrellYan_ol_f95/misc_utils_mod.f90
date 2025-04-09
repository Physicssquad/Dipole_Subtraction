module misc_utils_mod
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: dot

contains

  !==========================================================
  real(real64) function dot(p, q) result(d)

    ! Lorentz scalar product: dot(p, q) = p0*q0 - p1*q1 - ...
    real(real64), intent(in) :: p(0:3), q(0:3)
!    real(real64) :: d
    d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot
  !==========================================================

  !==========================================================
!  real(real64) function dot2(p,i,j,l) result(ans)
!
!      real(real64), intent(in) :: p(0:3,3,2)
!      integer,      intent(in) :: i,j,l
!
!      ans = p(0,i,l)*p(0,j,l) - p(1,i,l)*p(1,j,l) - p(2,i,l)*p(2,j,l) - p(3,i,l)*p(3,j,l)
!  end function dot2
  !==========================================================


end module misc_utils_mod

