program Drell_yan_dipoleSubtraction 
  use iso_fortran_env
  use globals_mod
  use printframes
  use constants_mod

! LOCAL PARAMETERS
    integer            :: iprint
    real(real64)       :: ans, sd, chi2
    character(len=100) :: run_tag, filename, mode
    real(real64), allocatable :: ai_nlo3(:), err(:)

! EXTERNAL FUNCTIONS
    external fnlo3
    external dipole_uU_g

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      read (10,*) pt_real         ! vegas points     LO 2 body
      read (10,*) its_real         ! vegas iterations LO 2 body
      npt_real= pt_real
      close(10)

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge      ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) pdf_name_NLO                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      read (15,*) iprint                ! to save data in output file         
      close(15)
      s=ecm*ecm
      

! ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        
      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*) filename
      close(20)

      if(iprint == 1) call output(run_tag,filename)
! ~~~~~~~~~~~~~~~~[--------------------------]~~~~~~~~~~~~~~~~~~~c        

        call init_mass
        call initpdfsetbyname(pdf_name_NLO)
        Call initPDF(0)
        call ol_LO_init
        call ol_NLO_1real_init

          print*,"  ----------------------------------"
          print*,"  |Initializing Dipole Subtraction  |"
          print*,"  ----------------------------------"
          print*," "
          print*," "
          print*," "
          call printframe1(pt_real,its_real)
          allocate(ai_nlo3(it_max), err(it_max)) 


          xq = xq_initial
          do j = 1,it_max
            
          call printframe2(xq)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call brm48i(40,0,0) 
            call vsup(6,npt_real,its_real,fnlo3,ans,sd,chi1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ai_nlo3(j) = ans
               err(j) = sd

            mode  = "[real-dipole]"
            call printframe3(mode,ans,sd,chi2)
            xq = xq + step_size 
          enddo

 ! Print results to stdout
  xq = xq_initial
  call printframe4(mode)
  do j = 1, it_max
     write(*,'(i7,3e27.15)') int(xq), ai_nlo3(j), err(j)
     xq = xq + step_size
  end do

  ! Write to file
  if (iprint == 1) then
     xq = xq_initial
     open(unit=20, file='../summary/'//trim(run_tag)//'/'//trim(filename), status='unknown')
     do j = 1, it_max
        write(20,*) xq, ai_nlo3(j), err(j)
        xq = xq + step_size
     end do
     close(20)
  end if
end program Drell_yan_dipoleSubtraction
