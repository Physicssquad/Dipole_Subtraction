module parameters
 implicit none
 integer, parameter :: max_size = 50

 double precision :: distr1_LO(max_size),integral1_LO(max_size),error1_LO(max_size)
 double precision :: distr2_LO(max_size),integral2_LO(max_size),error2_LO(max_size)

 double precision :: distr1_NLO(max_size),integral1_NLO(max_size),error1_NLO(max_size)
 double precision :: distr2_NLO(max_size),integral2_NLO(max_size),error2_NLO(max_size)

 double precision :: distr1_real(max_size),integral1_real(max_size),error1_real(max_size)
 double precision :: distr2_real(max_size),integral2_real(max_size),error2_real(max_size)

 double precision :: distr1_virtual(max_size),integral1_virtual(max_size),error1_virtual(max_size)
 double precision :: distr2_virtual(max_size),integral2_virtual(max_size),error2_virtual(max_size)

 double precision :: distr1_PK(max_size),integral1_PK(max_size),error1_PK(max_size)
 double precision :: distr2_PK(max_size),integral2_PK(max_size),error2_PK(max_size)

 double precision :: distr1_Plus(max_size),integral1_Plus(max_size),error1_Plus(max_size)
 double precision :: distr2_Plus(max_size),integral2_Plus(max_size),error2_Plus(max_size)

 double precision :: distr1_regular(max_size),integral1_regular(max_size),error1_regular(max_size)
 double precision :: distr2_regular(max_size),integral2_regular(max_size),error2_regular(max_size)

 double precision :: distr1_delta(max_size),integral1_delta(max_size),error1_delta(max_size)
 double precision :: distr2_delta(max_size),integral2_delta(max_size),error2_delta(max_size)

 character(len=100) :: run_tag
 character(len=100) :: real_dipole,virtual,PK,LO,LO_ref,NLO_ref
 character(len=100) :: regular_dat,plus_dat,delta_dat 

 integer :: test_LO,test_LO_ref
 integer :: test_real,test_real_ref
 integer :: test_PK,test_PK_ref
 integer :: test_virtual,test_virtual_ref
 integer :: test_plus,test_plus_ref
 integer :: test_regular,test_regular_ref
 integer :: test_delta,test_delta_ref

 integer :: it_max,i

end module parameters
!.................................................

