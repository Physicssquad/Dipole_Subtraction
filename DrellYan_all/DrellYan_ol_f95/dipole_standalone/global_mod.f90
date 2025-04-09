!#################################################
module globals_mod
    use iso_fortran_env
    implicit none
    real(real64) :: s
    real(real64) :: ge

! OPENLOOOPS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)

! VEGAS PARAMETERS
    real(real64) :: pt_real 
    integer :: its_real,npt_real

! MACHINE PARAMETERS
    real(real64) :: xq, xq_initial, step_size, ecm
    character(len=50)  :: pdf_name_NLO
    integer :: it_max
  
end module globals_mod
!#################################################
