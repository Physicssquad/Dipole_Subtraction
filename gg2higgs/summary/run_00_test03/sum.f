      program sum_columns
       implicit none
       integer, parameter :: n = 7  ! Number of rows in the data
       real(8) :: xq(n), plus(n), regular(n), delta(n), total(n)
       integer :: i
       character(len=50) :: input_file, output_file
   
       ! Input and output file names
       input_file = "PK_mod.dat"
       output_file = "output.txt"
   
       ! Open input file and read the data
       open(unit=10, file=input_file, status="old", action="read")
       do i = 1, n
           read(10, *) xq(i), plus(i), regular(i), delta(i)
       end do
       close(10)
   
       ! Compute the total for each row
       do i = 1, n
           total(i) = plus(i) + regular(i) + delta(i)
       end do
   
       ! Open output file and write the results
       open(unit=20, file=output_file, status="unknown", action="write")
       write(20, '(A)') "   xq         Total"
       do i = 1, n
           write(20, '(F10.0, F15.5)') xq(i), total(i)
       end do
       close(20)
   
       print *, "Data has been processed and written to", output_file
   
      end program sum_columns

