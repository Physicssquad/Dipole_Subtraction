c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_LO_init(id_LO)
        use openloops
        implicit none
        integer :: id_LO, error, k,icall
        real(8) :: m2_tree, m2_loop(0:2), acc
        real(8) :: p_ex(0:3,3),p1(0:3),p2(0:3),p3(0:3)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
c      common/usedalpha/alpha_s,ge
      common/energy/energy
      common/amass/am1,am2,mH,am4,am5
      common/caller/icall
      
         mu = mH/2d0
c	alpha_s = 125d0/2d0
      !  call set_parameter("order_ew", 1)
      !  call set_parameter("order_qcd", 2)
        if ( icall .eq. 0d0 ) then
        call set_parameter("mass(25)", mH)
        call set_parameter("verbose",-1)
        call set_parameter("model","heft") 
c        call set_parameter("alpha_s", alpha_s) !  AL will be set at the integral.f
        call set_parameter("mu", mu)

        id_LO = register_process("21 21 -> 25", 1)
        endif
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_NLO_real_init(id_NLO_1r)
        use openloops
        implicit none
        integer :: id_NLO_1r, error, k,icall
        real(8) :: m2_tree, m2_loop(0:2), acc
        real(8) :: p_ex(0:3,4),p1(0:3),p2(0:3),p3(0:3),p4(0:3)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
c      common/usedalpha/alpha_s,ge
      common/energy/energy
      common/amass/am1,am2,mH,am4,am5
      
      mu = mH/2d0
c	if ( icall .eq. 0d0 ) then
        call set_parameter("mass(25)", mH)
        call set_parameter("verbose",-1)
        call set_parameter("model","heft") 
c        call set_parameter("alpha_s", mu)
        call set_parameter("mu", mu)

        id_NLO_1r = register_process("21 21 -> 25 21", 1)
c	endif
      
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_get_NLO_1loop(p_ex,alpha_s,energy,mu,m2_tree)
        use openloops
        implicit none
        integer :: id, error, k
        real(8) :: m2_tree, m2_loop(0:2), acc
        real(8) :: p_ex(0:3,3)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy
      
      
      ! Set QCD and EW orders
      !  call set_parameter("order_ew", 1)
      !  call set_parameter("order_qcd", 2)
      
        call set_parameter("mass(25)", mH)
      
        ! Increase verbosity level to list loaded libraries
        call set_parameter("verbose",-1)
        call set_parameter("model","heft") 

        ! second argument of register_process:
        ! 1 for tree-like matrix elements (tree, color and spin correlations),
        ! 11 for loop, 12 for loop^2
        !id = register_process("1 -1 -> 23 2 -2", 1)
        id = register_process("21 21 -> 25", 11)
      
        ! start
        call start()
      
        if (id > 0) then
          ! set strong coupling
          call set_parameter("alpha_s", alpha_s)
          ! set renormalisation scale
          call set_parameter("mu", mu)

          ! evaluate tree matrix element
          call evaluate_tree(id, p_ex, m2_tree)
        end if
      
        call finish()
      
      end
