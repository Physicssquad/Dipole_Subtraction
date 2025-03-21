      program intPK
      use openloops
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,headline,headline1,headline2,headline3,headline4
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg
      external flo2_PlusA,flo2_PlusB,flo2_PlusC
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      common/prc_id/id_LO,id_NLO_1r
      common/cuts_delta/delta

      dimension PKPlus(1:50),err_Plus(1:50)
      dimension PKReg(1:50),err_Reg(1:50)
      dimension PKDel(1:50),err_Del(1:50)
      dimension PK(1:50),err(1:50)
      

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      do i=1,12
      read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      close(10)
      

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge          ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      read (15,*) iprint        ! to print data in file
      read (15,*) eps                   ! epsilon for bin width
      close(15)

      aem=1.0D0/128.0D0
      am1=0.0d0
      am2=0.0d0
      amH=125d0
      am4=0d0
      am5=0d0
 
c Parameters:
       xmuf = amH
       xmur = xmuf   !...being used in all common blocks.
      delta = 1d-5

      iselect_plus=0
      iselect_Regu=1
      iselect_Delt=1
c~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../output_files.dat',status='unknown')
      do i=1,2   
      read (20,*) 
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[P and K terms from Here ] 
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)

c ... Openloops initialization
        call ol_LO_init(id_LO)
c        call ol_NLO_real_init(id_NLO_1r)

       s=ecm*ecm

      headline = "P and K terms"
      call printframe0(headline)
      call printframe6(ecm,xmur,xmuf,name,amH)

        ! HERE RESET ALL THE VALUES TO INITIALISE
        do l = 1,it_max
          PKReg(l)    = 0d0
          err_Reg(l)  = 0d0 
          PKDel(l)    = 0d0
          err_Del(l)  = 0d0
          PKPlus(l)   = 0d0
          err_plus(l) = 0d0 
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus  functions ]
      if (iselect_Plus .eq. 1) then

      headline1 = "[+] distribution"
      headline2 = "PlusA distribution"
      headline3 = "PlusB distribution"
      pt2 = pt1
      npt2 = pt2
      its2 = its1

      call printframe0(headline1)
      print*,"PlusA - PlusB"
      call printframe1(pt2,its2)   ! Prints Vegas points

        do j=1,it_max

c      -------------------------------------------------
         call printframe0(headline2)
         call brm48i(40,0,0) 
         call vsup(2,npt2,its2,flo2_PlusA,ai_lo2A,sdA,chi2)

         call printframe0(headline3)
         call brm48i(40,0,0) 
         call vsup(2,npt2,its2,flo2_PlusB,ai_lo2B,sdB,chi2)
c      -------------------------------------------------
         ai_lo2 = (ai_lo2A - ai_lo2B)

c...     ai_lo2 =-(ai_lo2A - ai_lo2B)
c... this overall negative sign is a problem with the code I
c... am runing, for the comparison I am doing it remove soon.
c... surprisingly with negative sign values are matching
c... this is wrong by the way.

         PKPlus(j)   = ai_lo2
         err_plus(j) = sdA + sdB

      headline = "Plus"
      call printframe3(headline,ai_lo2,sd,chi2)   

        xq=xq + xincr
        enddo

        xq = xq_initial
      call printframe4a(headline)

        do j=1,it_max
          write(*,'(i7,3e27.15)')int(ecm),PKPlus(j),err_plus(j)
c         xq = xq + xincr
        enddo

	endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ regular functions ]
	if (iselect_Regu .eq. 1) then

      !input data card for vegas on regular and delta functions
      open(unit=10,file='../run.vegas.dat',status='unknown')
      do i=1,6
      read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      close(10)
 
        xq = xq_initial

      headline = "Regular Terms "
      call printframe0(headline)
      call printframe1(pt1,its1)

        do j=1,it_max

c     -------------------------------------------------
      call brm48i(40,0,0) 
      call vsup(2,npt1,its1,flo2_PKReg,ai_lo2,sd,chi2)
c     -------------------------------------------------

          PKReg(j) = ai_lo2
          err_Reg(j)=sd

      headline = "Regular"
         call printframe3(headline,ai_lo2,sd,chi2)

           xq=xq + xincr
        enddo
        xq = xq_initial
      call printframe4a(headline)

        do j=1,it_max
          write(*,'(i7,3e27.15)')int(ecm),PKReg(j),err_Reg(j)
        enddo

	endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ delta functions ]
	if (iselect_Delt .eq. 1 ) then

        xq = xq_initial

      headline = "Delta Functions"
      call printframe0(headline)
      call printframe1(pt1,its1)

        do j=1,it_max

      call printframe2(xq)

c     -------------------------------------------------
        call brm48i(40,0,0) 
        call vsup(1,npt1,its1,flo2_PKDel,ai_lo2,sd,chi2)
c     -------------------------------------------------

          PKDel(j) = ai_lo2
          err_Del(j)=sd

      headline = "Delta"
      call printframe3(headline,ai_lo2,sd,chi2)

           xq=xq + xincr
        enddo
        xq = xq_initial

       call printframe4a(headline)
        do j=1,it_max
          write(*,'(i7,3e27.15)')int(ecm),PKDel(j),err_Del(j)
        enddo

	endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Combining All ]
        xq = xq_initial
      print*,"  "
      write(*,*)achar(27)//'[1;32m'//
     . "   xq              Plus                      regular     
     .            delta                    combined PK     
     .    error",achar(27) //'[0m'

        do j=1,it_max
        PK(j) = PKPlus(j) + PKReg(j) + PKDel(j)
        err(j) = err_Plus(j) + err_Reg(j) + err_Del(j)
          write(*,'(i7,3e27.15,3e27.15,3e27.15,3e27.15,3e27.15)')
     .    int(xq),PKPlus(j),PKReg(j),PKDel(j),PK(j),err(j)
          xq = xq + xincr
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      


c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .ne. 1 ) goto 123
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PK(i),err(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
