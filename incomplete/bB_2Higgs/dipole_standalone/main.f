      program gg2aa_dipoleSubtraction 
      implicit double precision (a-h,o-z)
      integer leg
      character*50 green
      parameter (pi=3.14159265358979d0)
      dimension ai_nlo3(0:50),xerr(0:50)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/leg_choice/leg
      common/usedalpha/AL,ge
      common/set/set1
      common/countc/n4
      common/distribution/xq      
      common/t_cuts/e_cut,t_cut


      character*50 name,mode
      character*100 run_tag,filename
      external fnlo3

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

c      open(unit=10,file='../param_card.dat',status='unknown')    
c      read (10,*) ge      ! [ 1/Alpha_ew ]
c      close(10)
       ge = 7.8125000000000000E-003


      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      read (15,*) iprint                ! to save data in output file         
      read (15,*) eps                   ! epsilon for bin width
      read (15,*) e_cut                 ! cut used for small e5 Technical_cut 
      read (15,*) t_cut                 ! for cos(t)ij close to 1 Technical_cut

      close(15)

      aem=1.0D0/128.0D0

      

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        
      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*) filename
      close(20)
      if(iprint .eq. 1) call output(run_tag,filename)
c ~~~~~~~~~~~~~~~~[--------------------------]~~~~~~~~~~~~~~~~~~~c        

        call initpdfsetbyname(name)
        Call initPDF(0)
      am1=0d0
      am2=0d0
      am3=125d0
      am4=0d0
      am5=0d0
      leg=0

      ! energy
      s=ecm*ecm

c      print*,'  '
c      print*,"Press 1 to initialise VEGAS:"
c      print*,"Press 2 to initialise CUBA-VEGAS:"
c      read*,I
       I=1
      if (I .eq. 1 .or. I .eq. 2 ) then
       continue
      else
             print*,"Invalid Input"
             stop
      endif

c        i=2
          print*,"  ----------------------------------"
          print*,"  |Initializing Dipole Subtraction  |"
          print*,"  ----------------------------------"
          print*," "
          print*," "
          print*," "

          call printframe1(pt1,its1)
c         call printframe5(e_cut,t_cut)



          xq = xq_initial
          do j = 1,it_max
            
c          call printframe2(xq)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (I .eq. 1 ) then
            call brm48i(40,0,0) 
            call vsup(4,npt1,its1,fnlo3,ans,sd,chi2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elseif(I .eq. 2) THEN
                CALL cubacheck(ans,sd)
        endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ai_nlo3(j) = ans
               xerr(j) = sd

            mode  = "[real-dipole]"
            call printframe3(mode,ans,sd,chi2)
            xq = xq + step_size 
          enddo
          xq = xq_initial
         
          call printframe4(mode)
          do j=1,it_max
          write(*,'(i7,3e27.15)')int(xq),ai_nlo3(j),xerr(j)
          
          xq = xq + step_size 
          enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .          //trim(filename),status='unknown')
c     .          //trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_nlo3(i),xerr(i)
          xq = xq + step_size 
         enddo
         close(20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
123         continue
       end

