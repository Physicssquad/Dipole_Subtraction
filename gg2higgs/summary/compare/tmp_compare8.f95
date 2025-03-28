program compare_combine_and_display
 use parameters

! READING SUMMARY IMFORMATION AND CHECKING THE AVAILABILITY
 call read_machine_data('../../run.machine.dat')
 call read_filenames('../../output_files.dat')
 call check_directory()
 call read_machine_data_main('../' // trim(run_tag) // '/run.machine.dat')
 call print_run_info()

! PRELEMINARY CHECKS ARE COMPLETED GOING TO MAIN DATA HANDELING PART

!............................................................................................
 ! Leading Order
 call check_and_read_file(LO, distr1_LO, integral1_LO, error1_LO, test_LO)
 call check_and_read_file(LO_ref, distr2_LO, integral2_LO, error2_LO, test_LO_ref )
!............................................................................................

!............................................................................................
 ! Real- Dipole
 call check_and_read_file(real_dipole, distr1_real, integral1_real, error1_real, test_real)
 call check_and_read_file('real_ref.dat', distr2_real, integral2_real, error2_real, test_real_ref)
!............................................................................................

!............................................................................................
 ! Virtual 
 call check_and_read_file(virtual, distr1_virtual, integral1_virtual, error1_virtual,test_virtual)
 call check_and_read_file('virtual_ref.dat', distr2_virtual, integral2_virtual, error2_virtual,test_virtual_ref)
!............................................................................................

!............................................................................................
 !PK terms
100 continue
 call check_and_read_file(PK, distr1_PK, integral1_PK, error1_PK, test_PK)
 call loading
 if (test_PK == 0 ) goto 100
call check_and_read_file('PK_ref.dat', distr2_PK, integral2_PK, error2_PK, test_PK_ref)
!............................................................................................

! We have stored all the available data, now time to play around with them



 call play_around_with_data(LO, distr1_LO, integral1_LO, distr2_LO, integral2_LO)
! call play_around_with_data(virtual, distr1_virtual, integral1_virtual, distr2_virtual, integral2_virtual)
! call play_around_with_data(real_dipole, distr1_real, integral1_real, distr2_real, integral2_real)
! call play_around_with_data(PK, distr1_PK, integral1_PK, distr2_PK, integral2_PK)

end program compare_combine_and_display
!..........................................................................................................

subroutine play_around_with_data(identifier, distr_tmp1, integral_tmp1,distr_tmp2, integral_tmp2)
use parameters
character(len=*),INTENT(IN) :: identifier
double precision, dimension(max_size) :: distr_tmp1, integral_tmp1
double precision, dimension(max_size) :: distr_tmp2, integral_tmp2
integer :: j,tester1,tester2

    if ( trim(identifier) == trim(real_dipole)) then
    tester1  = test_real
    tester2  = test_real_ref
    elseif ( trim(identifier) == trim(virtual)) then
    tester1  = test_virtual
    tester2  = test_virtual_ref
    elseif ( trim(identifier) == trim(PK)) then
    tester1  = test_PK
    tester2  = test_PK_ref
    elseif ( trim(identifier) == trim(LO)) then
    tester1  = test_LO
    tester2  = test_LO_ref
    else
      PRINT *, "Unknown identifier: ", identifier
    endif

IF (tester1 + tester2 == 2.0) THEN
   PRINT *, " "
   PRINT *, "Finding ratio of ", TRIM(identifier), " with reference"
   PRINT *, " "

   ! Display Our Data
   PRINT *, "Our Data:"
   DO j = 1, it_max
      PRINT *, INT(distr_tmp1(j)), integral_tmp1(j)
   ENDDO
   PRINT *, " "

   ! Display Reference Data
   PRINT *, "Reference Data:"
   DO j = 1, it_max
      PRINT *, INT(distr_tmp2(j)), integral_tmp2(j)
   ENDDO
   PRINT *, " "

   ! Display Header for Comparison
   WRITE(*,*) ACHAR(27) // '[1;32m' // '     Distribution       Ratio' // ACHAR(27) // '[0m'

   ! Perform Comparison and Calculate Ratio
   DO j = 1, it_max
      IF (distr_tmp1(j) == distr_tmp2(j)) THEN
         PRINT *, INT(distr_tmp1(j)), "        ", integral_tmp1(j) / integral_tmp2(j) 
!         PRINT *, INT(distr_tmp1(j)), "        ", (integral_tmp1(j) - integral_tmp2(j)) / integral_tmp2(j) * 10000.0d0
      ENDIF
   ENDDO

ELSE
   PRINT *, "Cannot compare - File not present or mismatch in data."
ENDIF
end subroutine play_around_with_data











!..........................................................................................................
subroutine loading
use parameters, only : test_PK,PK
if (test_PK == 0 ) then
print *, "PK data does not exist."
print *, "Checking for individual data files..."
do i = 1, 3
  write(*,'(A1)', advance='no') '.'
 call sleep(1) ! You can adjust the duration or remove it if unavailable
end do
 call read_combine_PK_data(PK)
print *, " Done!"
  call sleep(1)
 end if
end subroutine loading


!..........................................................................................................
subroutine read_machine_data(filepath)
  use parameters, only : run_tag
  character(len=*) :: filepath
  open(unit=15, file=filepath, status='unknown')
  read (15,*) 
  read (15,*)
  read (15,*) 
  read (15,*) 
  read (15,*) 
  read (15,*)
  read (15,*) run_tag
  close(15)
end subroutine read_machine_data
!..........................................................................................................

!..........................................................................................................
subroutine read_machine_data_main(filepath)
  use parameters, only : run_tag,it_max
  character(len=*) :: filepath
  open(unit=15, file=filepath, status='unknown')
  read (15,*) 
  read (15,*)
  read (15,*) 
  read (15,*) it_max 
  read (15,*) 
  read (15,*)
  read (15,*) run_tag
  close(15)
 call check_machine_file(filepath)
end subroutine read_machine_data_main
!..........................................................................................................

!..........................................................................................................
subroutine read_filenames(filepath)
  use parameters, only : real_dipole,virtual,PK,LO,LO_ref,NLO_ref,regular_dat,plus_dat,delta_dat 
  character(len=*) :: filepath
  open(unit=16, file=filepath, status='unknown')
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*)
      read (16,*) real_dipole
      read (16,*) virtual
      read (16,*) PK
      read (16,*) LO
      read (16,*) LO_ref
      read (16,*) plus_dat
      read (16,*) regular_dat
      read (16,*) delta_dat
      read (16,*) NLO_ref
      close(16)
end subroutine read_filenames

subroutine check_directory()
  use parameters, only : run_tag
  call execute_command("test -d ../" // trim(run_tag) // " && echo 1 > command.txt || echo 0  > command.txt", ierr1)
  if (ierr1 == 0) then
    print *, 'Directory not found: ', run_tag
    print *, 'Currently stored data are in:'
    call system('cd ../ && ls')
    stop
  endif
end subroutine check_directory

!..........................................................................................................
subroutine check_machine_file(filename)
  use parameters, only : run_tag
  character(len=100) :: message
  character(len=*) :: filename
  integer :: ierr1
  call execute_command("test -f " // trim(filename) // " && echo 1 > command.txt || echo 0  > command.txt", ierr1)
  if (ierr1 == 0) then
      call system("cp ../../run.machine.dat ../"// trim(run_tag))
       print*,"Enter short text about this run for future references."
       print*,"Run Discription:"
       read(*, '(A)', advance='yes')message
       call system("echo " //message// " >> " //trim(filename))
  endif
end subroutine check_machine_file
!..........................................................................................................
!..........................................................................................................
subroutine read_combine_PK_data(PK_data_file)
 use parameters
 character(len=100) :: filename1,filename2,filename3,filename4
 character(len=*) :: PK_data_file

 filename1 = '../../PK_Isolated/summary/'//trim(run_tag)//'/'//trim(regular_dat)
 filename2 = '../../PK_Isolated/summary/'//trim(run_tag)//'/'//trim(plus_dat)
 filename3 = '../../PK_Isolated/summary/'//trim(run_tag)//'/'//trim(delta_dat)
 filename4 = '../'//trim(run_tag)//'/'//trim(PK_data_file)
 print*,"I hope filenames are correct"
 print*,filename1
 print*,filename2
 print*,filename3
 print*," Data from these files will be combined to"
 print*,filename4

 call check_and_read_file(filename1, distr1_Plus, integral1_Plus, error1_Plus,test_Plus)
 call check_and_read_file(filename2, distr1_Regular, integral1_Regular, error1_Regular,test_regular)
 call check_and_read_file(filename3, distr1_Delta, integral1_Delta, error1_Delta,test_delta)
 if ( test_plus + test_regular + test_delta == 3 ) then 

 print*," Combining Data into output "
do i = 1, 3
  write(*,'(A1)', advance='no') '.'
  call sleep(1) ! You can adjust the duration or remove it if unavailable
end do
print *, " Done! Combined to "
print*,PK_data_file
  call sleep(1)



 do i=1,it_max
    integral1_PK(i) = integral1_Plus(i) + integral1_regular(i) + integral1_Delta(i)
    distr1_PK(i) = distr1_Plus(i)
    error1_PK(i) = error1_Plus(i) + error1_regular(i) + error1_Delta(i)
 enddo
 open(unit=17,file=trim(filename4),status='unknown')
 do i=1,it_max
 write(17,'(i7,3e27.15)')int(distr1_PK(i)),integral1_PK(i),error1_PK(i)
 enddo
 close(17)
 else
 print*,"                  ⚠️"
           Print*,"Opps one or more filenames are not correct"
           print*,"or maybe some of the files doesn't exist "
           print*,"Take a closer look 👁️ ^👁️ "
           print*,"                  "
           stop
 endif
end subroutine read_combine_PK_data


subroutine print_run_info()
  use parameters, only :  run_tag
  print *, 'Reading Data from directory: /summary/' // trim(run_tag)
end subroutine print_run_info

!..........................................................................................................
subroutine execute_command(cmd, result)
  character(len=*) :: cmd
  integer :: result
  call system(cmd)
  open(unit=13, file='command.txt', status='unknown')
  read(13,*) result
  close(13)
  call system('rm command.txt')
end subroutine execute_command
!..........................................................................................................

!..........................................................................................................
subroutine check_and_read_file(file_name, distr_tmp, integral_tmp, error_tmp,checked)
use parameters, only : max_size,run_tag,it_max
character(len=*) :: file_name
double precision, dimension(max_size) :: distr_tmp, integral_tmp, error_tmp
integer :: i,ierr,ios,checked

call execute_command("test -f ../" // trim(run_tag) // "/" // trim(file_name) // &
  "&& echo 1 > command.txt || echo 0 > command.txt", ierr)
checked  = 0
if (ierr == 1) then
   open(unit=17, file='../' // trim(run_tag) // '/' // trim(file_name), status='unknown')
   do i = 1, it_max
      read(17, *,IOSTAT=ios) distr_tmp(i), integral_tmp(i), error_tmp(i)
! If error data is unavailable, set error to zero
      if (ios /= 0) then
        error_tmp(i) = 0.0d0
      endif
   end do
   close(17)

   print *, '/' // trim(file_name)
   write(*,*)achar(27) //'[1;32m'// 'distribution      Integral                     error(sd)'// achar(27) // '[0m'
   do i = 1, it_max
     write(*,'(i7,3e27.15)') int(distr_tmp(i)), integral_tmp(i),error_tmp(i)
!     print*,int(distr_tmp(i)), integral_tmp(i),error_tmp(i)
   end do
checked = 1
endif
end subroutine check_and_read_file
!..........................................................................................................
