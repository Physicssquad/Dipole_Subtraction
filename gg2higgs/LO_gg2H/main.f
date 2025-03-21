      program gg2H_LO 
      implicit double precision (a-h,o-z)
      double precision lambda
      dimension x(10),y(10),ai_lo2(1:50),err(0:50)
      parameter (pi=3.14159265358979d0)

      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/mass/amh
c      common/param2/xmur
      common/caller/icall
      common/scales/xmuf,xmur

c--------------------------------------------
c     common blocks used in couplings.f  
      common/add_par/xms,nd
      common/add_par1/acut
      common/rs_par/aam1,c0,aamh
      common/unpar/xl3,xdu,xlamu
      common/xmcoeff/xc1,xc2
      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      common/chfile/fname8
      common/isub/io,is
      common/max_order/iorder
      common/prc_id/id_LO,id_NLO_1r


      character*50 pdf_name,mode
      character*100 run_tag,filename
      external flo1_LO

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      do i=1,9
      read (10,*)
      enddo
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge       ! [ 1/Alpha_ew ]
      read (10,*) xmuf     ! Higgs mass  
      read (10,*) xmur     ! Higgs mass  
      read (10,*) amh      ! Higgs mass  
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) pdf_name          !lhapdf set
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) step_size         !step
      read (15,*) run_tag           !run_tag saves the output 
      read (15,*) iprint            !save data in output file ../summary
      close(15)

cc ~~~~~~~~~~~~~~~~[files needed by couplings.f]~~~~~~~~~~~~~~~~~~~c        
c
c      open(unit=20,file='../slicing_files/run.param.dat',
c     .    status='unknown')
c      read (20,*) nf            ! No. of flavours
c      read (20,*) ipdfs1        ! LO pdf set
c      read (20,*) xlqcd1        ! LO L_QCD5
c      read (20,*) ipdfs2        ! NLO pdf set
c      read (20,*) xlqcd2        ! NLO L_QCD5
c      close(20)
c
c      open(unit=30,file='../slicing_files/run.add.dat',status='unknown')
c      read (30,*) xms            ! M_s Fundamental Planck scale
c      read (30,*) nd             ! number of extra dimensions, 2<d<6
c      read (30,*) acut           ! \Lambda = acut*M_s
c      close (30)

      aem=1.0D0/128.0D0
      lambda = xlqcd1

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        
      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*) 
      read (20,*) 
      read (20,*) 
      read (20,*) filename
      close(20)
      if(iprint .eq. 1) call output(run_tag,filename)
c ~~~~~~~~~~~~~~~~[--------------------------]~~~~~~~~~~~~~~~~~~~c        

        call initpdfsetbyname(pdf_name)
        Call initPDF(0)
        call ol_LO_init(id_LO)
c        icall = 0
      
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=125d0
      am4=0d0
      am5=0d0
c      amh = 125.0d0
      leg=0
      ! energy
      s=ecm*ecm

        print*," "
        print*," "
        print*,"____________________________________"
        Print*," Calculating LO_gg2H"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        print*," "

      call printframe1(pt1,its1)   ! vegas points print 
      call printframe6(ecm,xmur,xmuf,pdf_name,am3)
      call system("sleep 1")

c      print*,'  '
c      print*,"Press 1 to initialise VEGAS:"
c      print*,"Press 2 to initialise CUBA-VEGAS:"
c      read*,int_choice
	int_choice = 1
        xq = xq_initial
        do j=1,it_max

c        call printframe2(xq)

	if (int_choice .eq. 1 ) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call brm48i(40,0,0) ! initialize random number generator
         call vsup(1,npt1,its1,flo1_LO,ans,sd,chi2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	elseif( int_choice .eq. 2 ) then
	call cubacheck(ans,sd)
	endif
            ai_lo2(j) = ans
              err(j)  = sd

            mode  = "Leading Order"
            call printframe3(mode,ans,sd,chi2)
         xq = xq + step_size
        enddo

         xq = xq_initial

         call printframe7(mode)

        do j=1,it_max
c          write(*,'(i7,3e27.15)')int(xq),ai_lo2(j),err(j)
          write(*,'(i7,3e27.15)')int(dsqrt(s)),ai_lo2(j),err(j)
          xq = xq + step_size
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i),err(i)
          xq = xq + step_size
         enddo
         close(20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

123         continue
c
       end
