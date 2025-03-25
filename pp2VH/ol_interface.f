c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_LO_init(id_LO)
        use openloops
        implicit none
        integer :: id_LO, error, k,icall
        real(8) :: m2_tree, m2_loop(0:2), acc
        real(8) :: p_ex(0:3,3),p1(0:3),p2(0:3),p3(0:3)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
        integer :: a,b
      common/usedalpha/alpha_s,ge
      common/energy/energy
      common/amass/am1,am2,mH,am4,am5
      common/caller/icall
      
!         mu = mH/2d0
c	alpha_s = 125d0/2d0
!        call set_parameter("order_ew",2)
!        call set_parameter("order_qcd",2)
c        if ( icall .eq. 0d0 ) then
!        call set_parameter("mass(25)", mH)
        call set_parameter("verbose",-5)
        call set_parameter("model","sm_yuksel") 

	call set_parameter("order_ew",2)  
c... what is the order_ew value, it doesnt matter weather 
c... i use 0,1,2 |loop_induced|^2 is same 

	call set_parameter("order_qcd",2)
c        call set_parameter("alpha_s", alpha_s) !  AL will be set at the integral.f
c        call set_parameter("muren",250d0)

        id_LO = register_process("21 21 -> 23 25",12)
        return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ol_NLO_real_init(id_NLO_1r)
      use openloops
       implicit none
       integer :: id_NLO_1r, error, k,icall
       real(8) :: m2_tree, m2_loop(0:2), acc
       real(8) :: p_ex(0:3,4),p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       real(8) :: mH = 125.0
       real(8) :: mu, alpha_s, energy,ge,am1,am2,am3,am4,am5
       common/energy/energy
       common/amass/am1,am2,am3,am4,am5
      
       call set_parameter("mass(25)", am3) !Higgs
       call set_parameter("mass(23)", am4) ! Z boson
       call set_parameter("verbose",-1)
       call set_parameter("model","sm") 
!      call set_parameter("muren", mur)
c... this can be called from integrand itself as it is a dynamic value

        id_NLO_1r = register_process("21 21 -> 25 23 21", 1)
      
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     subroutine ol_get_NLO_1loop(p_ex,alpha_s,energy,mu,m2_tree)
      subroutine ol_NLO_1loop_init(id_NLO_1loop)
      use openloops
      implicit none
      integer :: id_NLO_1loop, error, k
      real(8) :: m2_tree, m2_loop(0:2), acc
      real(8) :: p_ex(0:3,3)
      real(8) :: mH = 125.0
      real(8) :: mu, alpha_s, energy,ge,am1,am2,amH,am4,am5
      common/amass/am1,am2,mH,am4,am5
      
      ! Set QCD and EW orders
      !  call set_parameter("order_ew", 1)
      !  call set_parameter("order_qcd", 2)

      mu = mH/2d0 
      
        ! Increase verbosity level to list loaded libraries
      call set_parameter("mass(25)", mH)
      call set_parameter("verbose",-1)
      call set_parameter("model","heft") 
      call set_parameter("mu", mu)

      id_NLO_1loop = register_process("21 21 -> 25", 11)
      return
      end
