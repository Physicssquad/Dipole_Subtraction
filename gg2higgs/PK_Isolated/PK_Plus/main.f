      program intPK_Plus
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,mode,mode1
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg,xint_PlusA,xint_Plus_h
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/amass/am1,am2,amH,am4,am5
      common/prc_id/id_LO,id_NLO_1r
      common/scales/xmuf,xmur


      dimension PKPlus(1:50),err_Plus(1:50)
      dimension PKReg(1:50),err_Reg(1:50)
      dimension PKDel(1:50),err_Del(1:50)
      dimension PK(1:50),err(1:50)
      

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

      aem=1.0D0/128.0D0
      am1=0.0d0
      am2=0.0d0
      amH=125d0
      am4=0d0
      am5=0d0
!      xmuf = amH
!      xmur = xmuf

      call ol_LO_init(id_LO)

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../../output_files.dat',status='unknown')
      do i = 1,5
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        


        xq_initial = ecm
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

      mode = "P and K terms"
      call printframe0(mode)
      iselect_integrand = 0
      iselect_modified  = 1
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus  functions ]
      if (iselect_integrand .ne. 1 ) goto 150
      mode1 = "[+] distribution performed at integrand level"

      call printframe0(mode1)
      call printframe1(pt1,its1)   ! Prints Vegas points

        do j=1,it_max

      call printframe2(xq)

c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(2,npt1,its1,flo2_Plus,ai_lo2,sd,chi2)
c      -------------------------------------------------
         PKPlus(j)   = ai_lo2
         err_plus(j) = sd

      mode = "Plus"
      call printframe3(mode,ai_lo2,sd,chi2)   

        xq=xq + xincr
        enddo

        xq = xq_initial
      call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKPlus(j),err_plus(j)
          xq = xq + xincr
        enddo

        xq = xq_initial

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      
c... From here we will be implementing the [+] with 
c... different priscription.
150    continue
      if (iselect_modified .ne. 1 ) goto 151
      mode1 = "[+] distribution performed with intXh(x) "


      call printframe0(mode1)
      call printframe1(pt1,its1)   ! Prints Vegas points

!        do j=1,it_max

      call printframe2(xq)

c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(2,npt1*5,its1,xint_PlusA,ai_lo2,sd,chi2)
c      -------------------------------------------------
         Answer = ai_lo2
c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(1,npt1/10,its1,xint_Plus_h,ai_lo2,sd,chi2)
c      -------------------------------------------------
         xinth = ai_lo2

         PKPlus(1)   = Answer + xinth 
         err_plus(1) = sd

      mode = "Plus modified"
      call printframe3(mode,PKPlus(1),sd,chi2)   
151     continue
c... Output is not yet feeding into the file, I just want to 
c... check the numerical result from the vegas output.



c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .ne. 1 ) goto 123
c       open(unit=21,file='../../summary/'//trim(run_tag)//
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PKPlus(i),err_Plus(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
