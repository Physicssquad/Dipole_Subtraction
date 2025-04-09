!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_LO_init
        use globals_mod, only : id_LO
        use constants_mod
        use openloops

      
        call set_parameter("order_ew", 2)
        call set_parameter("order_qcd", 0)

        call set_parameter("mass(23)", mass(23))
        call set_parameter("mass(24)", mass(24))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("ew_scheme",2)

        call set_parameter("verbose",-1)
        call set_parameter("model","sm") 
                                               ! OPENLOOPS MANUAL
!Page 16 :http://link.springer.com/10.1140/epjc/s10052-019-7306-2


       id_LO(2) = register_process("2 -2 -> 11 -11", 1)  ! u ū
       id_LO(1) = register_process("1 -1 -> 11 -11", 1)  ! d d̄

       id_LO(4) = register_process("4 -4 -> 11 -11", 1)  ! c c̄
       id_LO(3) = register_process("3 -3 -> 11 -11", 1)  ! s s̄

       id_LO(5) = register_process("5 -5 -> 11 -11", 1)  ! b b̄

      return
      end subroutine ol_LO_init
   !=============================
       subroutine ol_NLO_1real_init
        use globals_mod, only : id_NLO_1r
        use constants_mod
        use openloops

      
        call set_parameter("order_ew", 2)
        call set_parameter("order_qcd", 1)


        call set_parameter("mass(23)", mass(23))
        call set_parameter("mass(24)", mass(24))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("ew_scheme",2)

        call set_parameter("verbose",-1)
        call set_parameter("model","sm") 
                                               ! OPENLOOPS MANUAL
!Page 16 :http://link.springer.com/10.1140/epjc/s10052-019-7306-2


       id_NLO_1r(2) = register_process("2 -2 -> 11 -11 21", 1)  ! u ū
       id_NLO_1r(1) = register_process("1 -1 -> 11 -11 21", 1)  ! d d̄
       id_NLO_1r(4) = register_process("4 -4 -> 11 -11 21", 1)  ! c c̄
       id_NLO_1r(3) = register_process("3 -3 -> 11 -11 21", 1)  ! s s̄
       id_NLO_1r(5) = register_process("5 -5 -> 11 -11 21", 1)  ! b b̄

      return
      end subroutine ol_NLO_1real_init
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!       subroutine ol_NLO_real_init(id_NLO_1r)
!        use openloops
!        implicit none
!        integer :: id_NLO_1r, error, k,icall
!        real(8) :: m2_tree, m2_loop(0:2), acc
!        real(8) :: p_ex(0:3,4),p1(0:3),p2(0:3),p3(0:3),p4(0:3)
!        real(8) :: mH = 125.0
!        real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
!c      common/usedalpha/alpha_s,ge
!      common/energy/energy
!      common/amass/am1,am2,mH,am4,am5
!      
!      mu = mH/2d0
!c	if ( icall .eq. 0d0 ) then
!        call set_parameter("mass(25)", mH)
!        call set_parameter("verbose",-1)
!        call set_parameter("model","heft") 
!c        call set_parameter("alpha_s", mu)
!        call set_parameter("mu", mu)
!
!        id_NLO_1r = register_process("21 21 -> 25 21", 1)
!c	endif
!      
!      end
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!c     subroutine ol_get_NLO_1loop(p_ex,alpha_s,energy,mu,m2_tree)
!      subroutine ol_NLO_1loop_init(id_NLO_1loop)
!      use openloops
!      implicit none
!      integer :: id_NLO_1loop, error, k
!      real(8) :: m2_tree, m2_loop(0:2), acc
!      real(8) :: p_ex(0:3,3)
!      real(8) :: mH = 125.0
!      real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
!      common/amass/am1,am2,mH,am4,am5
!      
!      ! Set QCD and EW orders
!      !  call set_parameter("order_ew", 1)
!      !  call set_parameter("order_qcd", 2)
!
!      mu = mH/2d0 
!      
!        ! Increase verbosity level to list loaded libraries
!      call set_parameter("mass(25)", mH)
!      call set_parameter("verbose",-1)
!      call set_parameter("model","sm") 
!c      call set_parameter("mu", mu)
!
!      id_NLO_1loop = register_process("2 -2 -> 11 -11", 11)
!      return
!      end
