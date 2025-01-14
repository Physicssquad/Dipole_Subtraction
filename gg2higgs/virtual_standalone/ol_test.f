
! Example  program how to use the native interface of OpenLoops.
! It calculates the Tree and loop matrix element of the process
! d dbar -> Z u ubar for a random phase-space point.

       subroutine main(p_ex,alpha_s,energy,mu,m2_tree,m2_loop)
        use openloops
        implicit none
        integer :: id, error, k
        real(8) :: m2_tree, m2_loop(0:2),m2_ir(0:2),acc,m2ir(0:2)
        real(8) :: p_ex(0:3,3)
        real(8) :: mH = 125.0
        real(8) :: mu, alpha_s, energy,m212,m210
      
      
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
        id = register_process("21 21 -> 25", 11)
      
        ! start
        call start()
      
        if (id > 0) then
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
!          call evaluate_tree(id, p_ex, m2_tree)
          ! print tree result
c          print *
c          print *, "evaluate_tree"
c          print *, "Tree:       ", m2_tree
      
          ! evaluate tree and loop matrix elements
        call evaluate_loop(id, p_ex, m2_tree, m2_loop,acc)
          ! print loop result
          print *
          print *, "evaluate_loop"
          print *, "Tree:       ", m2_tree
          print *, "Loop ep^0:  ", m2_loop(0)
          print *, "Loop ep^-1: ", m2_loop(1)
          print *, "Loop ep^-2: ", m2_loop(2)
          print *, "accuracy:   ", acc
          print *
c        call evaluate_iop(id, p_ex, m2_tree, m2_loop(0:2))
c          print *
c          print *, "evaluate_iop"
c          print *, "Tree:       ", m2_tree
c          print *, "Loop ep^0:  ", m2_loop(0)
c          print *, "Loop ep^-1: ", m2_loop(1)
c          print *, "Loop ep^-2: ", m2_loop(2)
c          print *, "accuracy:   ", acc
c          print *
!       call evaluate_iop(id, p_ex , m210,m2ir)
!           print*,"m212:",m212
!          print *, "Iopr ep^0:  ", m2_ir(0)
!          print *, "Iopr ep^-1: ", m2_ir(1)
!          print *, "Iopr ep^-2: ", m2_ir(2)
!          print*," "
!	stop

        end if
      
        ! finish
        call finish()
      
      end subroutine main
