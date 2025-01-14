
! Example  program how to use the native interface of OpenLoops.
! It calculates the Tree and loop matrix element of the process
! d dbar -> Z u ubar for a random phase-space point.

       subroutine main(p_ex,alpha_s,energy,mu,m2_tree)
        use openloops
        implicit none
        integer :: id, error, k
        real(8) :: m2_tree, m2_loop(0:2), acc
        real(8) :: p_ex(0:3,4)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy
      
      
      ! Set QCD and EW orders
      !  call set_parameter("order_ew", 1)
      !  call set_parameter("order_qcd", 2)
      
        call set_parameter("mass(25)", mH)
      
        ! Increase verbosity level to list loaded libraries
        call set_parameter("verbose",-1)
        call set_parameter("model","heft") 
        ! register one-loop amplitude for process d dbar -> Z u ubar
        ! The "ppzjj" process library must be installed before via
        ! $ ./scons auto=ppzjj
        !
        ! second argument of register_process:
        ! 1 for tree-like matrix elements (tree, color and spin correlations),
        ! 11 for loop, 12 for loop^2
        !id = register_process("1 -1 -> 23 2 -2", 1)
        id = register_process("21 21 -> 25 21", 1)
      
        ! start
        call start()
      
        if (id > 0) then
!          ! generate a random phase space point with rambo
!      !    call phase_space_point(id, energy, p_ex)
!             p_ex(           0 ,1)=   66.160389500771672     
!             p_ex(           1 ,1)=   0.0000000000000000     
!             p_ex(           2 ,1)=   0.0000000000000000     
!             p_ex(           3 ,1)=   66.160389500771672     
!        
!             p_ex(           0 ,2)=   62.838603026614308     
!             p_ex(           1 ,2)=   0.0000000000000000     
!             p_ex(           2 ,2)=   0.0000000000000000     
!             p_ex(           3 ,2)=  -62.838603026614308     
!        
!             p_ex(           0 ,3)=   125.20250857416981     
!             p_ex(           1 ,3)=   4.3036554831809612E-002
!             p_ex(           2 ,3)=   0.0000000000000000     
!             p_ex(           3 ,3)=   7.1180264905390089     
!        
!             p_ex(           0 ,4)=   3.7964839532161676     
!             p_ex(           1 ,4)=  -4.3036554831809612E-002
!             p_ex(           2 ,4)=   0.0000000000000000     
!             p_ex(           3 ,4)=  -3.7962400163816454     
!      
          ! set strong coupling
          call set_parameter("alpha_s", alpha_s)
          ! set renormalisation scale
          call set_parameter("mu", mu)
!      
!          print *
!          print *, "Tree and loop matrix element of the process"
!          print *, "g g > h"
!          print *, "for the phase-space point"
!          do k = 1, size(p_ex,2)
!            print *, 'P[', int(k,1), '] =', p_ex(:,k)
!          end do
      
          ! evaluate tree matrix element
          call evaluate_tree(id, p_ex, m2_tree)
          ! print tree result
c          print *
c          print *, "evaluate_tree"
c          print *, "Tree:       ", m2_tree
      
          ! evaluate tree and loop matrix elements
!          call evaluate_loop(id, p_ex, m2_tree, m2_loop(0:2), acc)
!          ! print loop result
!          print *
!          print *, "evaluate_loop"
!          print *, "Tree:       ", m2_tree
!          print *, "Loop ep^0:  ", m2_loop(0)
!          print *, "Loop ep^-1: ", m2_loop(1)
!          print *, "Loop ep^-2: ", m2_loop(2)
!          print *, "accuracy:   ", acc
!          print *
        end if
      
        ! finish
        call finish()
      
      end subroutine main
