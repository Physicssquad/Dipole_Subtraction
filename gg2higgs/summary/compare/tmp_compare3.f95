module parameters
  implicit none
  integer, parameter :: max_size = 50
  double precision :: distr1(max_size), integral1(max_size), xqvir(max_size), xintvir(max_size)
  double precision :: xqreal(max_size), xintreal(max_size), xqPK(max_size), xintPK(max_size)
  double precision :: xqPKterm2(max_size), xintPKterm2(max_size), distr2(max_size), integral2(max_size)
  double precision :: xqch(max_size), xintch(max_size)
  integer :: mid, it_max, ierr1, ierr2
  double precision :: ecm, xq_initial, step_size
  character(len=100) :: run_tag
  character(len=50) :: name, firstfile, secondfile
end module parameters

program compare
  use parameters

  call read_machine_data('../../run.machine.dat')
  call check_directory()
  call read_machine_data('../' // trim(run_tag) // '/run.machine.dat')
  call print_run_info()

  ! File comparison setup
  firstfile = 'PK_all.dat'
  secondfile = 'PK_ref.dat'

  call check_and_read_file(firstfile, distr1, integral1, ierr1)
  call check_and_read_file(secondfile, distr2, integral2, ierr2)

  if (ierr1 + ierr2 == 2) then
    call print_ratio(firstfile, secondfile, distr2, integral1, integral2)
    call print_ratio(secondfile, firstfile, distr2, integral2, integral1)
  endif
end program compare

subroutine read_machine_data(filepath)
  use parameters
  character(len=*) :: filepath
  open(unit=15, file=filepath, status='unknown')
  read (15,*) mid
  read (15,*) ecm
  read (15,*) name
  read (15,*) it_max
  read (15,*) xq_initial
  read (15,*) step_size
  read (15,*) run_tag
  close(15)
end subroutine read_machine_data

subroutine check_directory()
  use parameters
  call execute_command("test -d ../" // trim(run_tag) // " && echo 1 > command.txt || echo 0  > command.txt", ierr1)
  if (ierr1 == 0) then
    print *, 'Directory not found: ', run_tag
    print *, 'Currently stored data are in:'
    call system('cd ../ && ls')
    stop
  endif
end subroutine check_directory

subroutine execute_command(cmd, result)
  character(len=*) :: cmd
  integer :: result
  call system(cmd)
  open(unit=13, file='command.txt', status='unknown')
  read(13,*) result
  close(13)
  call system('rm command.txt')
end subroutine execute_command

subroutine print_run_info()
  use parameters
  print *, 'Reading Data from directory: /summary/' // trim(run_tag)
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '            ecm:', int(ecm), '[GeV]'
  print *, '     LHApdfname: ', name
  print *, '         it_max:', int(it_max)
  print *, 'initial Q value:', int(xq_initial), '[GeV]'
  print *, '     step size :', int(step_size)
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
end subroutine print_run_info

subroutine check_and_read_file(file_name, xq, xint, ierr)
  use parameters
  character(len=*) :: file_name
   double precision, dimension(50) :: xq, xint
  integer :: i,ierr

  call execute_command("test -f ../" // trim(run_tag) // "/" // trim(file_name) // & 
  "&& echo 1 > command.txt || echo 0 > command.txt", ierr)

  if (ierr == 1) then
    open(unit=17, file='../' // trim(run_tag) // '/' // trim(file_name), status='unknown')
    do i = 1, it_max
      read(17, *) xq(i), xint(i)
    end do
    close(17)
    print *, '/' // trim(file_name)
    print *, '   xq        Integral ', file_name
    do i = 1, it_max
      write(*,'(i7,3e27.15)') int(xq(i)), xint(i)
    end do
  endif
end subroutine check_and_read_file

subroutine print_ratio(file1, file2, xq, xint1, xint2)
  use parameters
  integer :: i
  character(len=*) file1, file2
  double precision xq(50)
  double precision xint1(1:50),xint2(1:50)
  print *, '/' // trim(file1) // ' and  /' // trim(file2)
  print *, 'Ecm/mH        first / second'
  do i = 1, it_max
    write(*,'(i7,3e27.15)') int(xq(i)), xint1(i) / xint2(i)
  end do
end subroutine print_ratio

