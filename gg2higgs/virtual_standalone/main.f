      program uU2eE_Virtual 
      use openloops
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50),err(0:50)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,amH,am4,am5
      common/usedalpha/AL,ge
      common/distribution/xq
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      common/scales/xmuf,xmur


      character*50 name,mode
      character*100 run_tag,filename
      external flo2_Vir

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      do i=1,3
      read (10,*)
      enddo
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)


      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge      ! [ 1/Alpha_ew ]
      read (10,*) xmuf      ! [ 1/Alpha_ew ]
      read (10,*) xmur      ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          !lhapdf set
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) xstep         !step
      read (15,*) run_tag       ! dir name to save data
      read (15,*) iprint            !save data in output file ../summary
      close(15)


      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*)
      read (20,*) filename
      close(20)

!      xmuf = xmuf*0.5d0
!      xmur = xmur*0.5d0



      aem=1.0D0/128.0D0

        call initpdfsetbyname(name)
        Call initPDF(0)

c ... Openloops initialization
        call ol_LO_init(id_LO)
        call ol_NLO_real_init(id_NLO_1r)
        call ol_NLO_1loop_init(id_NLO_1loop)
      
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      amH=125d0
      am4=0d0
      am5=0d0
      leg=0
      ! energy
      s=ecm*ecm


c       writes data in output file
        if(iprint .eq. 1)  call output(run_tag,filename)

        mode = "virtual contribrtion"
        call printframe0(mode)
        xq = xq_initial



        call printframe1(pt1,its1)
        call printframe6(ecm,xmur,xmuf,name,amH)
        
        do j=1,it_max

         call brm48i(40,0,0) ! initialize random number generator
         call vsup(1,npt1,its1,flo2_Vir,ans,sd,chi2)
           ai_lo2(j) = ans
              err(j) = sd

         call printframe3(mode,ans,sd,chi2)
        enddo

        call printframe4a(mode)
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(ECM),ai_lo2(j),err(j)
        enddo


        if(iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
c     .  //trim(filename),status='unknown', access='append')
         do i=1,it_max
          write(20,*)ECM,ai_lo2(i),err(i)
         enddo
         close(20)

123     continue 
       end
