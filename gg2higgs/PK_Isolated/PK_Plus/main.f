      program intPK_Plus
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,headline,headline1,headline2
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg,xint_PlusA,xint_Plus_h
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/amass/am1,am2,amH,am4,am5
      common/prc_id/id_LO,id_NLO_1r
      common/scales/xmuf,xmur
      common/cuts_delta/delta
      external flo2_PlusA,flo2_PlusB

      !input data card
      open(unit=10,file='../../run.vegas.dat',status='unknown')
      do i=1,12     
        read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      close(10)
      
      open(unit=10,file='../../param_card.dat',status='unknown')    
      read (10,*) ge          ! [ 1/Alpha_ew ]
      read (10,*) xmuf
      read (10,*) xmur
      close(10)

      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      read (15,*) iprint        ! to print data in file
      close(15)

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../../output_files.dat',status='unknown')
      do i = 1,5
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
      aem=1.0D0/128.0D0
      am1=0.0d0
      am2=0.0d0
      amH=125d0
      am4=0d0
      am5=0d0
      delta = 1d-5

      call ol_LO_init(id_LO)

      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

      headline = "     [+] distribution     "
      call printframe0(headline)
c... ...... Plus implementation type
      iselect_integrand = 0
      iselect_modified  = 1
      iselect_integral  = 0

!      call printframe9(ecm)
      call printframe6(ecm,xmur,xmuf,name,amH)





      if (iselect_integrand .ne. 1 ) goto 150
      headline1 = "[+] with integrand level subtraction"
      call printframe1(pt1,its1)   ! Prints Vegas points
      call printframe0(headline1)
c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(2,npt1,its1,flo2_Plus,ai_lo2,sd,chi2)
c      -------------------------------------------------
         PKPlus   = ai_lo2
         err_plus = sd




150    continue
      if (iselect_modified .ne. 1 ) goto 151
      call printframe1(pt1,its1)   ! Prints Vegas points
      headline1 = "[+] with intXh(x) "

         call printframe0(headline1)
         call brm48i(40,0,0) 
         call vsup(2,npt1,its1,xint_PlusA,ai_lo2_a,sd,chi2)

         call printframe1(pt1,its1)   ! Prints Vegas points
         call brm48i(40,0,0) 
         call vsup(1,npt1,its1,xint_Plus_h,ai_lo2_h,sd,chi2)

         PKPlus   = ai_lo2_a + ai_lo2_h 
         err_plus = sd






151    continue
       if (iselect_integral .ne. 1 ) goto 152
      headline1 = "PlusA distribution"
      headline2 = "PlusB distribution"
         call printframe1(pt1*10d0,its1)   ! Prints Vegas points

         call printframe0(headline1)
         call brm48i(40,0,0)
         call vsup(2,npt1*10,its1,flo2_PlusA,ai_lo2A,sdA,chi2)

         call printframe0(headline2)
         call brm48i(40,0,0)
         call vsup(2,npt1*10,its1,flo2_PlusB,ai_lo2B,sdB,chi2)

         PKPlus = (ai_lo2A - ai_lo2B)
         err_plus = sdA + sdB







152     continue
      headline = "Plus"
      call printframe3(headline,PKPlus,err_plus,chi2)   
      call printframe4(headline)
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(ecm),PKPlus,err_plus

c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .ne. 1 ) goto 123
c       open(unit=21,file='../../summary/'//trim(run_tag)//
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
          write(21,*)ecm,PKPlus,err_Plus
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
