      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50),err(0:50)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL,ge
      common/distribution/xq

c--------------------------------------------
c     common blocks used in couplings.f  
      common/add_par/xms,nd
      common/add_par1/acut
      common/rs_par/aam1,c0,aamh
      common/unpar/xl3,xdu,xlamu
      common/xmcoeff/xc1,xc2
c      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      common/chfile/fname8
      common/isub/io,is
      common/max_order/iorder
      common/param/aem,xmur,lambda

      common/cone/ET_iso,r0,rgg  ! this is for the cone part.
      common/counter/ifilter,itot_ev,iselect_scale
      common/counter_diff/diff,eps

c--------------------------------------------




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


c ~~~~~~~~~~~~~~~~[files needed by couplings.f]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../slicing_files/run.param.dat',
     .    status='unknown')
      read (20,*) nf            ! No. of flavours
      read (20,*) ipdfs1        ! LO pdf set
      read (20,*) xlqcd1        ! LO L_QCD5
      read (20,*) ipdfs2        ! NLO pdf set
      read (20,*) xlqcd2        ! NLO L_QCD5
      close(20)

      open(unit=30,file='../slicing_files/run.add.dat',status='unknown')
      read (30,*) xms            ! M_s Fundamental Planck scale
      read (30,*) nd             ! number of extra dimensions, 2<d<6
      read (30,*) acut           ! \Lambda = acut*M_s
      close (30)

      open(unit=50,file='../slicing_files/run.cone.dat',
     . status='unknown')
      read (50,*) ET_iso       ! ET_iso in GeV
      read (50,*) r0           ! r0
      read (50,*) rgg          ! r_gamma_gamma
      read (50,*) niso         ! n value in Frixione's algorithm
      close (50)



      aem=1.0D0/128.0D0
      lambda = xlqcd1





        call initpdfsetbyname(name)
        Call initPDF(0)
      
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=0d0
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
        do j=1,it_max

        call printframe2(xq)

         call brm48i(40,0,0) ! initialize random number generator
         call vsup(3,npt1,its1,flo2_Vir,ans,sd,chi2)
           ai_lo2(j) = ans
              err(j) = sd

         call printframe3(mode,ans,sd,chi2)

        xq = xq + xstep
        enddo

        xq = xq_initial

        call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),ai_lo2(j),err(j)
          xq = xq + xstep
        enddo


        if(iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
c     .  //trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i),err(i)
          xq = xq + xstep
         enddo
         close(20)

123     continue 
       end
