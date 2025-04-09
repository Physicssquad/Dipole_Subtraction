module matrix_elements_mod
  use openloops
  use globals_mod, only : ge, id_LO, id_NLO_1r
  use misc_utils_mod
  use constants_mod
  implicit none
  private

  public :: evaluate_born,evaluate_1real

contains

  subroutine evaluate_born(p1, p2, p3, p4,answer)
    real(8), intent(out) :: answer(5) 
    real(8), intent(in) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)

    ! Internal variables
    real(8) :: s12, s13, s23
    real(8) ::  e
    real(8) :: p_ex(0:3,4)
!    real(8) :: qu2,gew, qe, qu, cw, sw, cw2, sw2, zprop
    integer :: i

    ! Initialize coupling
    e = dsqrt(ge * 4.d0 * PI)

    ! Mandelstam variables
    s13 = 2.d0 * dot(p1, p3) ! t
    s23 = 2.d0 * dot(p1, p4) ! u
    s12 = 2.d0 * dot(p1, p2) ! s

!    qe = e
!    qu = 1d0 
!    qu2 = qu**2 
!    zprop = 1.d0 / (s - mZ**2)
!      BornAmp = (2*e**4*qu2*(-2*s13*s23 + s12*(s13 + s23)))/(3d0*s12**2)

! Set OpenLoops parameters and evaluate
!    call SET_PARAMETER('alpha_s', AL)
!    call SET_PARAMETER('mu', scale)

    call p1d_to_p2d_4(p1, p2, p3, p4, p_ex)
    do i = 1,5
    call evaluate_tree(id_LO(i), p_ex, answer(i))
    enddo

  end subroutine evaluate_born

  subroutine evaluate_1real(AL, p1, p2, p3, p4, p5,answer)
    real(8), intent(out) :: answer(5) 
    real(8), intent(in) :: AL, p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)

    ! Internal variables
    real(8) ::  e
    real(8) :: p_ex(0:3,5)
!    real(8) :: qu2,gew, qe, qu, cw, sw, cw2, sw2, zprop
    integer :: i

    ! Initialize coupling
    e = dsqrt(ge * 4.d0 * PI)

! Set OpenLoops parameters and evaluate
    call SET_PARAMETER('alpha_s', AL)
!    call SET_PARAMETER('mu', scale)

    call p1d_to_p2d_5(p1, p2, p3, p4 ,p5, p_ex)
    do i = 1,5
    call evaluate_tree(id_NLO_1r(i), p_ex, answer(i))
    enddo
  end subroutine evaluate_1real

end module matrix_elements_mod
