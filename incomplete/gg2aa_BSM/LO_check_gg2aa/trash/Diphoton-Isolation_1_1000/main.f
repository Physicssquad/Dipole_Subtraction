C----------------------------------------------------------------------C
C-                                                                    -C
C-  This code calculates d\sigma/dQ for digamma production.           -C
C-  One can choose to calculate for LHC or TEVATRON.                  -C
C-  Calculation is done using phase-space slicing method, and         -C
C-  uses Frixione's algorithm to completely remove fragmentation      -C
C-  contribution.                                                     -C
C     -                     CTEQ6M                                    -C
C-              NLO pdf used throughout.                              -C
C-    Renormalization and factorization scales are identifield        -C
C-    and set equal to Q.                                             -C
C-                                                                    -C
C----------------------------------------------------------------------C 
C----------------------------------------------------------------------C

      program diphoton
      implicit double precision (a-h,o-z)
      double precision lambda

      parameter (pi=3.14159265358979d0)
      character * 25 machine,pdf_name
      character * 25 pref,pref1,pref2,pref3,pref4
      character * 25 prefa,prefb,prefc,prefd 
      character * 25 fname,fname1,fname2,fname3,fname4,fname5 
      character * 25 fname6,fname7,fname8
      character * 6 ch_model, sub_process
      character * 8 ch_order, ch_pdfs

      common/machine/mid

      common/pdflag/iset
      common/energy/s
      
      common/totordiff/crosssec
      common/slice/deltas,deltac
      common/corr/icorr

      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/angle/cststar

      common/model/model
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/add_par/xms,nd
      common/add_par1/acut
      common/rs_par/aam1,c0,aamh
      common/nflavour/nf
      common/unpar/xl3,xdu,xlamu
      common/iparallel/ilabel
      common/xmcoeff/xc1,xc2
      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      common/chfile/fname8
      common/isub/io,is
      common/max_order/iorder

      integer iw(2)

      dimension sqrs(25)
      dimension smvr(3)
      character chvar(85)*3

      dimension yfdn(10)
      dimension cfdn(20)

      data sqrs/1.0D0,3.0D0,5.0D0,7.0D0,10.0D0,15.0D0,17.0D0,
     &        18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0,24.0D0,
     &        28.0D0,32.0D0,36.0D0,38.0D0,39.0D0,40.0D0,41.0D0,
     &        42.0D0,43.0D0,44.0D0,45.0D0/

       data smvr/0.5D0,1.D0,1.5D0/

        data chvar/ 'a1',  'a2',   'a3', 'a4',  'a5',  'a6',  'a7',
     &      'a8',  'a9',  'a10', 'a11', 'a12', 'a13', 'a14', 'a15',
     &      'a16', 'a17', 'a18', 'a19',
     &      'a20', 'a21', 'a22', 'a23', 'a24', 'a25', 'a26',
     &      'a27', 'a28', 'a29', 'a30', 'a31', 'a32', 'a33', 'a34',
     &      'a35', 'a36', 'a37', 'a38', 'a39', 'a40', 'a41', 'a42',
     &      'a43', 'a44', 'a45', 'a46', 'a47', 'a48', 'a49',
     &      'a50', 'a51', 'a52', 'a53', 'a54', 'a55', 'a56',
     &      'a57', 'a58', 'a59', 'a60', 'a61', 'a62', 'a63', 'a64',
     &      'a65', 'a66', 'a67', 'a68', 'a69', 'a70', 'a71', 'a72',
     &      'a73', 'a74', 'a75', 'a76', 'a77', 'a78', 'a79', 'a80',
     &      'a81', 'a82', 'a83', 'a84', 'a85'/

      data cfdn/-0.95D0,-0.85D0,-0.75D0,-0.65D0,-0.55D0,
     &          -0.45D0,-0.35D0,-0.25D0,-0.15D0,-0.05D0,
     &           0.05D0, 0.15D0, 0.25D0, 0.35D0, 0.45D0,
     &           0.55D0, 0.65D0, 0.75D0, 0.85D0, 0.95D0/

C----------------------------------------------------------------------C
C Arrays for various subprocesses.  These five dimensional arrays are
C written lo,lo2,nlo2,nlo3 wise.  They will be of the form
C sub(mq,mf,id,io,is)
C mq : for distribution like Q
C mf : for scale variation like \mu_F
C id : for variation of slicing parameters, \delta_s, \delta_c
C io : subprocess origin, 1: SM; 2: BSM direct; 3: SM*BSM interference
C is : subprocess, 1 for qqb; 2 for qg; 3 for gg
C----------------------------------------------------------------------C

      dimension sub_lo(25,3,15,3,3)
      dimension sub_lo2(25,3,15,3,3)
      dimension sub_nlo2(25,3,15,3,3)
      dimension sub_nlo3(25,3,15,3,3)

      external flo2
      external fnlo2
      external fnlo3

       tag=50
       dest=0

      aem=1.0D0/128.0D0

      !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
	npt1 = pt1
      read (10,*) its1          ! vegas iterations LO 2 body
      read (10,*) npt2          ! vegas points     NLO 2 body
      read (10,*) its2          ! vegas iterations NLO 2 body
      read (10,*) npt3          ! vegas points     NLO 3 body
      read (10,*) its3          ! vegas iterations NLO 3 body          
      close(10)

      open(unit=15,file='run.machine.dat',status='unknown')
      read (15,*) machine       ! LHC or TEVATRON
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      close(15)

      open(unit=20,file='run.param.dat',status='unknown')
      read (20,*) nf            ! No. of flavours
      read (20,*) ipdfs1        ! LO pdf set
      read (20,*) xlqcd1        ! LO L_QCD5
      read (20,*) ipdfs2        ! NLO pdf set
      read (20,*) xlqcd2        ! NLO L_QCD5
      close(20)

      open(unit=25,file='run.order.dat',status='unknown')
      read (25,*) iorder        ! 0:lo   1:nlo
      read (25,*) model         ! 0: SM 1: ADD  2: RS  3: Unparticle
      read (25,*) isw           ! Background SM Q variation for model 1,2or3
      read (25,*) id            ! id deltas:deltac
      read (25,*) imin          ! imin
      read (25,*) imax          ! imax
      read (25,*) icorr         ! type of correction
      read (25,*) afac          ! afac
      read (25,*) mqmin         ! mqmin
      read (25,*) mqmax         ! mqmax
      read (25,*) aqstep        ! aqstep
      read (25,*) mfmin         ! mfmin
      read (25,*) mfmax         ! mfmax
      read (25,*) mr            ! \mu_R index 1,2,3 for (0.5,1,1.5)*Q 
      read (25,*) xeps          ! xeps
      close(25)

       open(unit=50,file='run.cone.dat',status='unknown')
       read (50,*) ET_iso       ! ET_iso in GeV
       read (50,*) r0           ! r0
       read (50,*) rgg          ! r_gamma_gamma
       read (50,*) niso         ! n value in Frixione's algorithm
       close (50)


      if (model.eq.0) then
      prefa='sm'
      elseif (model.eq.1) then
      prefa='add'
      elseif (model.eq.2) then
      prefa='rs'
      elseif (model.eq.3) then
      prefa='up'
      endif

      if (iorder.eq.0) then
      prefb='qlo_'
      elseif (iorder.eq.1) then
      prefb='qnlo_'
      endif

        if (iorder.eq.1) then
                if (model.eq.0) then
                nprocess=2              ! For NLO \sigma with NLO pdfs
                elseif (model.eq.1.or.model.eq.2.or.model.eq.3) then
                nprocess=5
                endif
        elseif (iorder.eq.0) then
                nprocess=1
        else
                write(*,*)'Unable to assign the processors'
        stop
        endif

       ichavar = mqmin
       prefc=chvar(ichavar)
       prefd='.out'
       call strcat(prefa,prefb,fname6)
       call strcat(fname6,prefc,fname7)
       call strcat(fname7,prefd,fname8)

       do norder=0,iorder

      write (*,*) '____________________'
      write (*,*) 'machine      = ',machine
      write (*,*) 'machine id   = ',mid
      write (*,*) 'ecm          = ',ecm
      write (*,*) '____________________'
      write (*,*) ' ' 

      write (*,*) '______VEGAS_______' 
      write (*,*) 'pts1      = ',npt1
      write (*,*) 'its1      = ',its1
      write (*,*) 'pts2      = ',npt2
      write (*,*) 'its2      = ',its2
      write (*,*) 'pts3      = ',npt3
      write (*,*) 'its3      = ',its3
      write (*,*) '_________________ '
      write (*,*) ' '        
            
      write (*,*) '______MODEL_________'
      write (*,*) 'crosssec             =',crosssec
      write (*,*) 'model                = ',model
      write (*,*) 'SM B/G for model     =',isw
      write (*,*) 'nf                   = ',nf
      write (*,*) 'ipdfs1               = ',ipdfs1
      write (*,*) 'LO L_QCD5            = ',xlqcd1
      write (*,*) 'ipdfs2               = ',ipdfs2   
      write (*,*) 'NLO L_QCD5           = ',xlqcd2
      write (*,*) '____________________'
      write (*,*) ' '

      write (*,*) '______LO/NLO________'
      write (*,*) 'iorder     = ',iorder
      write (*,*) 'id        = ',id
      write (*,*) 'imin      = ',imin
      write (*,*) 'imax      = ',imax
      write (*,*) 'icorr     = ',icorr
      write (*,*) 'afac      = ',afac
      write (*,*) 'mqmin     = ',mqmin
      write (*,*) 'mqmax     = ',mqmax
      write (*,*) 'mfmin     = ',mfmin
      write (*,*) 'mfmax     = ',mfmax
      write (*,*) 'mr        = ',mr
      write (*,*) 'aqstep    = ',aqstep
      write (*,*) 'xeps      = ', xeps
      write (*,*) '___________________'
      write (*,*) ' '

      write (*,*) 'E_T iso      =', ET_iso
      write (*,*) 'cone radius  =', r0
      write (*,*) 'R_gg         =', rgg
      write (*,*) 'n in H(R)    =', niso
      write (*,*) '___________________'
      write (*,*) ' '

      if (model.eq.1) then
      open(unit=30,file='run.add.dat',status='unknown')
      read (30,*) xms            ! M_s Fundamental Planck scale
      read (30,*) nd             ! number of extra dimensions, 2<d<6
      read (30,*) acut           ! \Lambda = acut*M_s
      close (30)
      write (*,*) 'ADD model'
      write (*,*) 'M_s = ',xms,'GeV'
      write (*,*) 'ND=',nd
      write (*,*) 'acut=',acut

      elseif (model.eq.2 .or. isw.eq.2) then
      open(unit=35,file='run.rs.dat',status='unknown')
      read (35,*) aam1           ! M1 mass of the 1st excited RS mode
      read (35,*) c0_a           ! c0 effective RS coupling
      read (35,*) aamh           ! Higgs mass
      close (35)

      c0 = c0_a/dsqrt(8.0d0*pi)
c      c0 = c0_a

      write (*,*) 'RS model'
c      write (*,*) 'M1 = ',aam1,'GeV'
      write (*,*) 'C0 bar= ',c0_a
      write (*,*) 'M_H = ',aamh,'GeV'
      elseif (model.eq.3) then
      open(unit=40,file='run.up.dat',status='unknown')
      read (40,*) xlamu         ! Lambda_U 
      read (40,*) xdu           ! du
      read (40,*) xl3           ! lam_T
      read (40,*) xlams         ! lam_S 
      close(40)
      write (*,*) 'Lambda_U  = ',xlamu
      write (*,*) 'du        = ',xdu
      write (*,*) 'lam_T     = ',xl3
      endif

      pdf_name = 'MMHT2014nlo68cl'
      call initpdfsetbyname(pdf_name)
      Call initPDF(0)

        do mq = mqmin, mqmax

      ! energy
      s=ecm*ecm
      write(*,*) ' '

        if (norder.eq.0) then
        iset=ipdfs1
        lambda=xlqcd1           ! L_QCD5 for MRST 2001 LO pdfs
        ch_pdfs = 'LO PDFs'
        write(*,*)'iset=',iset
        write(*,*) 'LO L_QCD5=',lambda

                if (iset.eq.9)  write(*,*)'MRST 2001 LO PDFs'
                if (iset.eq.51) write(*,*)'CTEQ6L LO PDFs'
                if (iset.eq.71) write(*,*)'ALEKHIN LO PDFs (NO L_QCD5)'

        elseif (norder.eq.1) then
        iset=ipdfs2
        lambda=xlqcd2           ! L_QCD for MRST 2001 NLO pdfs
        ch_pdfs = 'NLO PDFs'
        write(*,*)'iset=',iset
        write(*,*) 'NLO L_QCD5=',lambda

                if (iset.eq.5)  write(*,*)'MRST 2001 NLO PDFs'
                if (iset.eq.52) write(*,*)'CTEQ6M NLO PDFs'
                if (iset.eq.72) write(*,*)'ALEKHIN NLO PDFs (NO L_QCD5)'

c  For ALEKHIN pdfs, \alphs_s will be supplied with the grid, hence 
c  L_QCD is not required for both LO and NLO pdfs

      endif

c     show what alphas this lambda gives at the Z mass

      zm=91.19d0
      xnf =nf

      alphas1 = ALFAS1(zm,lambda,xnf)
      alphas2 = ALFAS2(zm,lambda,xnf)

      write (*,*) '______ALPHAS________'
      write (*,*)'Alphas at Z mass'
      write (*,*)'Z mass   :',zm
      write (*,*)'one loop :',alphas1
      write (*,*)'two loop :',alphas2
      write (*,*) '___________________'

        do i = imin, imax

        if(norder.eq.0)then
c        i=imin
        write(*,*)'LO'
        elseif (norder.eq.1) then
c        i=imin

                if(id.eq.0)then
                deltas=(1.0d-6)*(10**i)
                deltac=1d-05
                elseif(id.eq.1)then
                deltas=1d-03
                deltac=(1d-6)*(5**i)
                elseif(id.eq.2)then
                deltas = (1.0d-5)*(10**i)
                deltac = deltas/afac
                else
                write(*,*) 'unknown id: ',id
                stop
                write(*,*)'NLO'
                endif

         write(*,*) '| deltas = ',deltas,' |'
         write(*,*) '| deltac = ',deltac,' |'
         write(*,*) ' '
         endif

            xc2 = smvr(mr)
            mf=mfmin
            xc1 = smvr(mf)


c             aam1 = mq*aqstep
c             xq =  mq*aqstep
                aam1 = mq + aqstep
                  xq = mq + aqstep

             xlow = xq - xeps
             xhigh = xq + xeps

            write(*,*)'M1 =', aam1
            write(*,*)' QF =',xq
            write(*,*)'Q_low ',xlow,' Q_high ',xhigh

c    ----------------------SM------------------------
        if (model.eq.0) then

                write (*,*) 'SM model distributions'

        elseif (model.eq.1 .or. model.eq.2 .or. model.eq.3) then

                if(model.eq.1) then
                write (*,*) 'ADD model distributions'
                elseif (model.eq.2) then
                write(*,*) 'RS model distributions'
                elseif (model.eq.3) then
                write(*,*) 'Unparticle distributions'
                endif

        endif         ! endif model

c        call class(iirank, io_min, io_max, is_min, is_max)

        io_min = 2
        io_max = 2

        is_min = 3
        is_max = 3
c        is_max = 2

        do io = io_min,io_max

        if (io .eq.1) then
        ch_model = 'SM'
        elseif (io.eq.2) then
        ch_model = 'BSM'
        elseif (io.eq.3) then
        ch_model = 'INTF'
        endif

        do is = is_min, is_max          ! 1:qqb, 2:qg, 3: gg

        if (is.eq.1) then
        sub_process = 'qqb'
        elseif  (is.eq.2) then
        sub_process = 'qg'
        elseif (is.eq.3) then
        sub_process = 'gg'
        endif

c--------------- LEADING ORDER ----------------------------------
        if(norder.eq.0) then            
c----------------------- Leading order with LO pdfs--------------
        write(*,*) ch_model, sub_process, 'LO--', ch_pdfs
c-------No qg contribution at LO --------------------------------
        if (is.eq.2) then
        sub_lo(mq,mf,i,io,is) = 0d0
        else
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(3,npt1,its1,flo2,ai_lo,sd,chi2)
        sub_lo(mq,mf,i,io,is) = ai_lo
        endif
c----------------------------------------------------------------

c================================================================
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c================================================================

c----------------NEXT-TO-LEADING ORDER --------------------------
        elseif(norder.eq.1)then
c----------------------- Leading order with NLO pdfs-------------
        write(*,*) ch_model, sub_process, 'LO2--', ch_pdfs
c-------No qg contribution at LO --------------------------------
        if (is.eq.2) then
        sub_lo2(mq,mf,i,io,is) = 0d0
        else
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(3,npt1,its1,flo2,ai_lo2,sd,chi2)
        sub_lo2(mq,mf,i,io,is) = ai_lo2
        endif
c----------------------------------------------------------------
c-------Next-to-leading order 2-body contribution ---------------
        write(*,*) ch_model, sub_process, 'NLO2--', ch_pdfs
c-------No SM gg NLO 2-body contribution-------------------------
        if ((io.eq.1 .or.  io.eq.3) .and. is.eq.3) then
        sub_nlo2(mq,mf,i,io,is) = 0d0
        else
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(3,npt2,its2,fnlo2,ai_nlo2,sd,chi2)
        sub_nlo2(mq,mf,i,io,is) = ai_nlo2
        endif
c----------------------------------------------------------------
c-------Next-to-leading order 3-body contribution ---------------
        write(*,*) ch_model, sub_process, 'NLO3--', ch_pdfs
c-------No SM gg 3-body contribution-----------------------------
        if ((io.eq.1 .or. io.eq.3) .and. (is.eq.3)) then
        sub_nlo3(mq,mf,i,io,is) = 0d0
        else
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(6,npt3,its3,fnlo3,ai_nlo3,sd,chi2)
        sub_nlo3(mq,mf,i,io,is) = ai_nlo3
        endif
c----------------------------------------------------------------
       else
       write(*,*) 'ERROR: unknown norder: ',norder
       stop
       endif
c================================================================

 111    continue


        enddo           ! i variation
        enddo           ! is variation
        enddo           ! io variation


       enddo            ! mq variation
       enddo            ! norder enddo 

       if (model.eq.0) then
       io_min = 1                      ! SM
       io_max = 1
       elseif (model.eq.1 .or. model.eq.2 .or. model.eq.3) then
       io_min = 2                      ! BSM
       io_max = 3                      ! SM*BSM
       endif

       do norder =0,iorder

        do io = io_min, io_max

       if(model.eq.0) then 
       pref='sm.lhc.'

       elseif (model.eq.1) then

       if (io.eq.2) then
       pref='add.lhc.'
       elseif (io.eq.3) then
       pref='add.int.lhc.'
       endif

       elseif (model.eq.2) then

       if (io.eq.2) then
       pref='rs.lhc.'
       elseif (io.eq.3) then
       pref='rs.int.lhc.'
       endif

       elseif (model.eq.3) then

       if (io.eq.2) then
       pref='up.lhc.'
       elseif (io.eq.3) then
       pref='up.int.lhc.'
       endif

       endif

       if (norder.eq.0) then
       pref1='lo.dat'
       elseif (norder.eq.1) then
       pref1='lo2.dat'
       pref2='nlo2.dat'
       pref3='nlo3.dat'
       endif

       if (norder.eq.0) then
       call strcat(pref,pref1,fname)
       open(unit=56,file=fname,status='unknown')

       elseif (norder.eq.1) then
       call strcat(pref,pref1,fname1)
       call strcat(pref,pref2,fname2)
       call strcat(pref,pref3,fname3)
       call strcat(pref,'tot.nlo.dat',fname4)
       call strcat(pref,'sub.nlo.dat',fname5)

       open(unit=57,file=fname1,status='unknown')
       open(unit=58,file=fname2,status='unknown')
       open(unit=59,file=fname3,status='unknown')
       open(unit=80,file=fname4,status='unknown')
       open(unit=85,file=fname5,status='unknown')
       endif


       mf = mfmin
       do i=imin,imax

        if(id.eq.0)then
                deltas=(1.0d-6)*(10**i)
                deltac=1d-05
        elseif(id.eq.1)then
                deltas=1d-03
                deltac=(1d-6)*(5**i)
        elseif(id.eq.2)then
                deltas = (1.0d-5)*(10**i)
                deltac = deltas/afac
        else
        write(*,*) 'unknown id: ',id
        stop
        endif

       do mq=mqmin,mqmax

 
c       aam1 = mq*aqstep
c       xq   = aam1

                aam1 = mq + aqstep
                  xq = mq + aqstep
c       xq = sqrs(mq)*aqstep + 300d0
c           cf = cfdn(mq)

       if (norder.eq.0) then

       tot_lo = sub_lo(mq,mf,i,io,1) 
     .        + sub_lo(mq,mf,i,io,2)
     .        + sub_lo(mq,mf,i,io,3)  

       write(56, 106) xq,model,norder,
     .  (sub_lo(mq,mf,i,io,is),is=1,3),tot_lo

       elseif (norder.eq.1) then
C Total LO2, NLO2 and NLO3 contributions.
       tot_lo2 = sub_lo2(mq,mf,i,io,1)
     .         + sub_lo2(mq,mf,i,io,2)
     .         + sub_lo2(mq,mf,i,io,3)

       tot_nlo2 = sub_nlo2(mq,mf,i,io,1)
     .          + sub_nlo2(mq,mf,i,io,2)
     .          + sub_nlo2(mq,mf,i,io,3)

       tot_nlo3 = sub_nlo3(mq,mf,i,io,1)
     .          + sub_nlo3(mq,mf,i,io,2)
     .          + sub_nlo3(mq,mf,i,io,3)

C Subprocess contribution by summing the difference PSSM pieces.
       tot_qqb = sub_lo2(mq,mf,i,io,1)
     .         + sub_nlo2(mq,mf,i,io,1) 
     .         + sub_nlo3(mq,mf,i,io,1)

       tot_qg = sub_lo2(mq,mf,i,io,2) 
     .        + sub_nlo2(mq,mf,i,io,2) 
     .        + sub_nlo3(mq,mf,i,io,2)

       tot_gg = sub_lo2(mq,mf,i,io,3) 
     .        + sub_nlo2(mq,mf,i,io,3)
     .        + sub_nlo3(mq,mf,i,io,3)

C Total contribution at NLO
       tot_nlo = tot_qqb + tot_qg + tot_gg

       write(57,106) xq, model, norder, 
     &  (sub_lo2(mq,mf,i,io,is),is=1,3), tot_lo2

       write(58,108) xq, model, norder, deltas, deltac,
     &  (sub_nlo2(mq,mf,i,io,is),is=1,3), tot_nlo2

       write(59, 108) xq, model, norder, deltas, deltac,
     &  (sub_nlo3(mq,mf,i,io,is),is=1,3), tot_nlo3

       write(80,108) xq, model, norder, deltas, deltac,
     &  tot_lo2, tot_nlo2, tot_nlo3, tot_nlo

       write(85,108)xq, model, norder, deltas, deltac,
     &  tot_qqb, tot_qg, tot_gg, tot_nlo

        endif           ! 'norder' endif
        enddo           !  'mq' enddo 
        enddo           !  'i' enddo 
        enddo           !  io variation for model

 106   format((f8.2),2(i2),4(g15.4))
 108   format((f8.2),2(i2),2(g15.2),4(g15.4))

       enddo       ! norder enddo

       close(56)
       close(57)
       close(58)
       close(59)
       close(80)


       stop
       end
