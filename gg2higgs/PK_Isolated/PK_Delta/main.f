      program intPK_Delta
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,mode,mode1
      character*100 command,run_tag,dir_path,filename,filename1
      dimension PKPlus(1:50),err_Plus(1:50)
      dimension PKReg(1:50),err_Reg(1:50)
      dimension PKDel(1:50),err_Del(1:50)
      dimension PK(1:50),err(1:50)
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/amass/am1,am2,amH,am4,am5
      common/prc_id/id_LO,id_NLO_1r
      common/scales/xmuf,xmur
      

      !input data card
      open(unit=10,file='../../run.vegas.dat',status='unknown')
      do i=1,6
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
      aem=1.0D0/128.0D0
      am1=0.0d0
      am2=0.0d0
      amH=125d0
      am4=0d0
      am5=0d0
c Parameters:
!       xmuf = amH
!       xmur = xmuf   !...being used in all common blocks.
      delta = 1d-5


c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../../output_files.dat',status='unknown')
      do i =1,7
      read (20,*) 
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[P and K terms from Here ] 
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
      call ol_LO_init(id_LO)
       s=ecm*ecm
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ delta functions ]
        xq = xq_initial

      mode = "Delta Functions"
      call printframe0(mode)
      call printframe1(pt1,its1)

        do j=1,it_max

      call printframe2(ecm)

c     -------------------------------------------------
        call brm48i(40,0,0) 
        call vsup(1,npt1/5,its1,flo2_PKDel,ai_lo2,sd,chi2)
c     -------------------------------------------------

          PKDel(j) = ai_lo2
          err_Del(j)=sd

      mode = "Delta"
      call printframe3(mode,ai_lo2,sd,chi2)

           xq=xq + xincr
        enddo
        xq = xq_initial

       call printframe4(mode)
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKDel(j),err_Del(j)
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
          write(21,*)xq,PKDel(i),err_Del(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
