C      Program test
C      IMPLICIT REAL*8(A-H,O-Z)
C      REAL*8 LAMBDA,NF
C      common/couplings/aem,alps,ak2D 
C      MODEL=3
C      RS=91.10D0 
C      LAMBDA=0.2D0 
C      NF=5.0D0
C      XMH=200.0D0
C      p34=0.9D0
C
C      call coupfact(model,p34,rs,lambda,xnf)
C      write(*,*) aem,alps,AK2D 
C
C      RETURN
C      END
*----------------------------------------------------------------------
C      call the following before you call subroutine SMqqb etc        C
C      call coupfact(model,p34,rs,lambda,xnf)                          C
*----------------------------------------------------------------------
*               BLACK BOX
*----------------------------------------------------------------------
      subroutine coupfact(model,rp34,AK2D,AK2DINTF)
      implicit double precision (a-h,o-z)

      call ak2dmodel(model,rp34,bk2d,bk2dintf)
      
      AK2D=bk2d
      AK2DINTF=bk2dintf
      return
      end


C THE ONE-LOOP RUNNING COUPLING CONSTANT

      REAL*8 FUNCTION ALFAS1(RS,LAMBDA,XNF)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA,LAMBD2
      data PI/3.141592653589793238462643D0/

      LAMBD2=LAMBDA*LAMBDA
      BF=(33.0D0-2.0D0*XNF)/12.0D0/PI
      RS2=RS*RS
      DL1=DLOG(RS2/LAMBD2)
      ALFAS1=1.0D0/DL1/BF
      RETURN
      END
C THE TWO-LOOP RUNNING COUPLING CONSTANT

      REAL*8 FUNCTION ALFAS2(RS,LAMBDA,XNF)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA,LAMBD2
      data PI/3.141592653589793238462643D0/
 
      LAMBD2=LAMBDA*LAMBDA
      BF=(33.0D0-2.0D0*XNF)/12.0D0/PI
      BFP=(153.0D0-19.0D0*XNF)/(66.0D0-4.0D0*XNF)/PI
      RS2=RS*RS
      DL1=DLOG(RS2/LAMBD2)
      DL2=DLOG(DL1)
      ALFAS2=(1.0D0-(BFP*DL2)/(BF*DL1))/DL1/BF
      RETURN
      END

C THE THREE-LOOP RUNNING COUPLING CONSTANT

      REAL*8 FUNCTION ALFAS3(RS,LAMBDA,NFL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA,NFL
      data PI/3.141592653589793238462643D0/
 
      AL=2.0D0*DLOG(RS/LAMBDA)
      AL2=AL*AL
      AL3=AL2*AL
      BETA0=(11-2.0D0*NFL/3.0D0)/4.0D0
      BETA02=BETA0*BETA0
      BETA03=BETA02*BETA0
      B1=(102.0D0-38.0D0*NFL/3.0D0)/BETA0/16.0D0
      B2=(2857.0D0/2.0D0-5033.0D0*NFL/18.0D0
     1+325.0D0*NFL*NFL/54.0D0)/BETA0/64.0D0
      ALFAS3=PI*(1.0D0/BETA0/AL-B1*DLOG(AL)/BETA02/AL2
     2+(B1*B1*(DLOG(AL)**2.0D0-DLOG(AL)-1.0D0)+B2)
     3/BETA03/AL3)
      RETURN
      END

       SUBROUTINE AK2DMODEL(MODEL,XMH,AK2D,AK2DINTF)
*      -----------------------------------
*      AK2D= kappa^2*D = 16 Pi 1/Ms^4 lambda

       IMPLICIT REAL*8 (A-H,O-Z)
       complex*8 lambdaf
       complex*8 sprop
       common/xpar/c00,ams,amh
       common/add_par/xms,nd
       common/add_par1/acut
       common/rs_par/aam1,c0,aamh
       common/unpar/xl3,xdu,xlamu
       common/param/aem,xmur,lambda
       common/sparam/smass,xlam
       common/scalar_couplings/cph,cgl,cfe
       common/sm_couplings/c1,c2,c3,ct,g1,g2,g3
       common/weak/amz,amw,sw2,cw2
       common/amfermion/amtau,amc,amb,amt

       data PI/3.141592653589793238462643D0/


c       write(*,*)'Ph, Gl, Fe =',cph, cgl, cfe
c       write(*,*)'c1,c2,c3,ct =',c1,c2,c3,ct
c       write(*,*)'g1,g2,g3 =',g1,g2,g3
c       write(*,*)'Mz, Mw, Mt, Sw2, Cw2 =', amz, amw, amt, sw2,cw2
c       write(*,*)'Scalar mass, Lambda =', smass, xlam


       XMZ=91.1876D0
       GMZ=2.4952D0
*
       XMH2=XMH*XMH
       XMZ2=XMZ*XMZ
       GMZ2=GMZ*GMZ
*
       ALEM=aem
       ALEM2=ALEM*ALEM
*
       SW2=sw2
       CW2=cw2
*
       QE=-1.0D0
       QU=2.0D0/3.0D0
       QD=-1.0D0/3.0D0
*
       GEV=-1.0D0/4.0D0-SW2*QE
       GUV=1.0D0/4.0D0-SW2*QU
       GDV=-1.0D0/4.0D0-SW2*QD
*
       GEA=1.0D0/4.0D0
       GUA=-1.0D0/4.0D0
       GDA=1.0D0/4.0D0
*
       SU=QE*QU+XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2) 
     &       /CW2/SW2*GUV*GEV
       SD=QE*QD+XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2) 
     &       /CW2/SW2*GDV*GEV
       SUA=XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2) 
     &       /CW2/SW2*GUA*GEA
       SDA=XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2) 
     &       /CW2/SW2*GDA*GEA


        XMS=xms
        ND=nd

        c00=c0
        amh=aamh
        am1=aam1
        am0=am1/3.83D0

       XDU=xdu
       MODL=MODEL

       IF (MODL.EQ.3) GO TO 300
       IF (MODL.EQ.2) GO TO 200
       IF (MODL.EQ.1) GO TO 100

 100   ams=XMS
       xlambda=acut*ams
       rms=xmh/xlambda
       xlambms=xlambda/ams

c  k^2D(Q^2)=[pi/2+i*\Lambda(M_s/Q)].  Hence in the direct contribution
c  there will be a factor of pi^2/4 on taking the complex absolute.
c  Modified on 04-04-2008.

c       xc = cabs(lambdaf(0,0,ND,rms,0))           ! (Direct)
c       xd=rmsa**(ND-2)
c       xpi2=PI/2*xd
c       abslam = 32.0D0*PI/ams/ams/ams/ams/2.D0
c       abslam=abslam*dsqrt(xpi2**2+xc**2)

       xc = cabs(lambdaf(0,0,ND,rms))           ! (Direct)
       xd = xlambms**(ND-2)
       abslam = 16.0D0*PI/ams/ams/ams/ams*xd
       abslam=abslam*xc

        FGRU=abslam

c "-ive" sign will be there for the integral to be consistent with the formula.
c Intereference: sign introduced which comes from the k^2D(Q^2)
c In the direct contribution sqaure will be taken and sign can be neglected.

c       abslam_intf = -32.0D0*PI*cabs(lambdaf(0,0,ND,rms,0))
c     &  /ams/ams/ams/ams
c       abslam_intf=abslam_intf/2.d0

       xc_intf = (-1.0d0)*real(lambdaf(0,0,ND,rms))
       xd = xlambms**(ND-2)
       abslam_intf=16.0D0*PI/ams/ams/ams/ams*xd
       abslam_intf=abslam_intf*xc_intf

       FGRU1=4.0D0*PI*ALEM*abslam_intf
       GO TO 110

  200  ams=am0
       rms=xmh/ams
c       abslam = 32.0D0*PI*C00**2*cabs(lambdaf(1,0,ND,rms,0))
c     &  /ams/ams/ams/ams
c       abslam=abslam/2.d0
       xc=cabs(lambdaf(1,0,ND,rms))
       abslam = 16.0D0*PI*C00**2/ams/ams/ams/ams
       abslam=abslam*xc
c       write(*,*)'am1,c0,amh,nd,rms',am1,c0,amh,ND,rms
       FGRU=abslam
c       abslam_intf = 32.0D0*PI*C00**2*cabs(lambdaf(1,0,ND,rms,1))
c     &  /ams/ams/ams/ams
c       abslam_intf=abslam_intf/2.d0
       xc_intf = real(lambdaf(1,0,ND,rms))
       abslam_intf = 16.0D0*PI*C00**2/ams/ams/ams/ams
       abslam_intf = abslam_intf*xc_intf
       FGRU1=4.0D0*PI*ALEM*abslam_intf
       GO TO 110

C  AK2D = kap^2*D   (D is propagator/i)

  300  sdw = SWidth(smass,xlam)

       SBWRe = xmh**2 - smass**2
       SBWIm = smass*sdw
       Den   = SBWRe**2 + SBWIm**2

       SRe =  SBWRe/Den
       SIm = -SBWIm/Den
        
       sprop = cmplx(SRe, SIm)

       FGRU = cabs(sprop)
       FGRU1 = 4.0d0*PI*ALEM*real(sprop)

c       write(*,*)'Scalar k^2 D(Q^2), Re(k^2 D(Q^2)) =',FGRU, FGRU1

       GO TO 110

c FGRU1 is the interference term containing e**2.  AK2D contains the propagator
c term only.  Its square will be taken in the matrix elements.

  110  AK2D=FGRU
       AK2DINTF=FGRU1

c       write(*,*)'alem,pi',ALEM,PI
c      write(*,*)'M_s,nd,rms,AK2D, AK2DINTF',ams,ND,rms,AK2D,
c     &          AK2DINTF/4/PI/ALEM

       RETURN
       END

C       |----------------------------------------------------|
C       | This function calculates the sum over Kaluza-Klein |
C       | graviton propagators in                            |
C       |   (a) the ADD model for both s and t channels      |
C       |   (b) the  RS model for both s and t channels      |
C       |----------------------------------------------------|
C       | Inputs: Switch  IModel = 0 for ADD; 1 for RS       |
C       |                 IChanl = 0 for s;   1 for t        |
C       |     Parameters  Dim for ADD Model                  |
C       |           read in through the Common block /Param/ |
C       | Variable:  X = Sqrt(s) / M_S for ADD Model, s-chanl|
C       |               (range 0 < X < 1)                    |
C       |              = Sqrt(-t)/ M_S for ADD Model, t-chanl|
C       |               (range 0 < X < 1)                    |
C       |              = Sqrt(s) / m_0 for  RS Model, s-chanl|
C       |               (range 0 < X < 1000)                 |
C       |              = Sqrt(-t)/ m_0 for  RS Model, t-chanl|
C       |               (range 0 < X < 100)                  |
C       | Index : Isep =0 for Direct BSM; 1 for SM*BSM       |
C       | Index "Isep" included on 05-04-2008.               |
C       |----------------------------------------------------|
c          Function Lambdaf (IModel, IChanl,Dim, X,Isep)
          Function Lambdaf (IModel, IChanl,Dim, X)
          Implicit None
          Integer  IModel, IChanl, Nfit, NRes, NRestart
          Integer  NUppr, NLowr, KReal, KImag
c          Integer  Isep 
          Complex*8  Lambdaf
          Real*8  X, FitY, Pi, EvenLog, OddLog
          Real*8  Xfit(200), Yfit(200) 
          Real*8  zJ1, BWigRe, BWigIm, BWReTm, BWImTm, GW, Arg1, Arg2
          Real*8     GWidth, DiGamma
          Integer        Dim
          Real*8             ams, c00, amh
          Real*8              m0, c0, mH
c         Common /Param/Dim,m0,c0,mH
          common/xpar/c00,ams,amh
C          common/dpar/dim
          External GWidth, DiGamma
          data Pi/3.141592653589793238462643D0/
          Data Xfit/ 
     >  0., 1.0D-02, 2.0D-02, 3.0D-02, 4.0D-02, 5.0D-02, 6.0D-02, 
     >  7.0D-02, 8.0D-02, 9.0D-02, 1.0D-01, 0.20, 0.30, 0.40, 0.50, 
     >  0.60, 0.70, 0.80, 0.90, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
     >  1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
     >  2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
     >  4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2,
     >  5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4,
     >  6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
     >  7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8,
     >  8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0,
     >  11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22.,
     >  23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34.,
     >  35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46.,
     >  47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58.,
     >  59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
     >  70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80.,
     >  81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91.,
     >  92., 93., 94., 95., 96., 97., 98., 99., 100./ 
          Data Yfit/ 
     > 12.5768, 12.5768, 12.5766, 12.5763, 12.5760, 12.5755, 12.5749,
     > 12.5742, 12.5735, 12.5726, 12.5716, 12.5559, 12.5299, 12.4938,
     > 12.4479, 12.3924, 12.3278, 12.2545, 12.1730, 12.0838, 11.9875,
     > 11.8846, 11.7758, 11.6616, 11.5426, 11.4194, 11.2926, 11.1627,
     > 11.0302, 10.8957, 10.7596, 10.6224, 10.4844, 10.34598, 10.2077,
     > 10.06949, 9.93197, 9.79527, 9.65963, 9.52523, 9.39227, 9.26088,
     > 9.13119, 9.00333, 8.87736, 8.75338, 8.63144, 8.51159, 8.39386,
     > 8.27828, 8.16485, 8.05359, 7.94449, 7.83754, 7.73274, 7.63005,
     > 7.52946, 7.43094, 7.33446, 7.23998, 7.14748, 7.05692, 6.96825,
     > 6.88145, 6.79647, 6.71327, 6.63182, 6.55207, 6.47399, 6.39753,
     > 6.32265, 6.24933, 6.17751, 6.10716, 6.03825, 5.97074, 5.90459,
     > 5.83977, 5.77625, 5.71398, 5.65295, 5.59311, 5.53444, 5.47690,
     > 5.42047, 5.36512, 5.31083, 5.25755, 5.20528, 5.15398, 5.10362,
     > 5.05419, 5.00567, 4.95802, 4.91122, 4.86527, 4.82012, 4.77578,
     > 4.73220, 4.68938, 4.64731, 4.60595, 4.56529, 4.52532, 4.48601,
     > 4.44736, 4.40935, 4.37197, 4.33519, 4.29900, 3.96703, 3.68164,
     > 3.43388, 3.21691, 3.02541, 2.85520, 2.70296, 2.56600, 2.44216,
     > 2.32964, 2.22699, 2.13295, 2.04649, 1.96674, 1.89295, 1.82448,
     > 1.76077, 1.70135, 1.64580, 1.59375, 1.54488, 1.49892, 1.45560,
     > 1.41472, 1.37606, 1.33945, 1.30474, 1.27178, 1.24045, 1.21061,
     > 1.18218, 1.15505, 1.12913, 1.10436, 1.08064, 1.05792, 1.03613,
     > 1.01523, 0.995147, 0.975845, 0.957277, 0.939401, 0.922180, 
     > 0.905579, 0.889565, 0.874106, 0.859176, 0.844746, 0.830793, 
     > 0.817293, 0.804225, 0.791567, 0.779302, 0.767411, 0.755876, 
     > 0.744684, 0.733818, 0.723264, 0.713009, 0.703041, 0.693348, 
     > 0.683918, 0.674741, 0.665808, 0.657107, 0.648631, 0.640371, 
     > 0.632318, 0.624465, 0.616805, 0.609330, 0.602035, 0.594912, 
     > 0.587955, 0.581160, 0.574519, 0.568029, 0.561683, 0.555478, 
     > 0.549408, 0.543469, 0.537658, 0.531969, 0.526399, 0.520945, 
     > 0.515603, 0.510369, 0.505240, 0.500213, 0.495286/
*--------------------------------------------------------------------*



*  Branch between models:          
*  ~~~~~~~~~~~~~~~~~~~~~                                              
          If (IModel .eq. 0) Then                        ! (ADD Model)
*                                                                     
*            Restrict to cases of 2 - 6 extra dimensions:
*            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             If (      Dim .lt. 2
     >            .or. Dim .gt. 6) Then
                       Lambdaf = Cmplx(0.d0,0.d0)
                       Return
             Endif
*                                                                     
*            Check that the argument 0 < X < 1:
*            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             If (      X  .lt. 0.d0
     >            .or. X  .gt. 1.d0 ) Then
                       Lambdaf = Cmplx(0.d0,0.d0)
                       Return
             Endif

*     ---------------------------------------------
*            Branch between channels:         
*            ~~~~~~~~~~~~~~~~~~~~~~~~                                        
              If (IChanl .eq. 0) Then                  ! (s - channel)
*                                                                             
                    EvenLog = Log ( 1.d0/X**2 - 1.d0 ) 
                    OddLog  = Log ( (1.d0 - X)/(1.d0 + X) )
*                                                       
                    If (Dim .eq. 2) Then                ! (case d = 2)
*                   ~~~~~~~~~~~~~~~~~~~~                                   
                        Lambdaf = 0.5d0 
     >                               * Cmplx 
     >                                (
     >                                  EvenLog,
     >                                  Pi
C     >                                  0.0D0 
     >                                )

                    Endif
*                                                         
                    If (Dim .eq. 3) Then                ! (case d = 3)
*                   ~~~~~~~~~~~~~~~~~~~~                             
                        Lambdaf = 0.5d0
     >                               * Cmplx 
     >                                (
     >                                  X * OddLog + 2.d0, 
     >                                  X * Pi
C     >                                  0.0D0
     >                                )
                    Endif
*                                                
                    If (Dim .eq. 4) Then                ! (case d = 4)
*                   ~~~~~~~~~~~~~~~~~~~~                             
                        Lambdaf = 0.5
     >                               * Cmplx 
     >                                (
     >                                  X*X * EvenLog + 1.d0, 
     >                                  X*X * Pi
C     >                                0.0D0
     >                                )
                    Endif
*                                                     
                    If (Dim .eq. 5) Then                ! (case d = 5)
*                   ~~~~~~~~~~~~~~~~~~~~                                    
                        Lambdaf = 0.5
     >                               * Cmplx
     >                                (
     >                           X**3 * OddLog + 2.d0/3.d0 + 2.d0*X*X,
     >                           X**3 * Pi
C     >                                0.0D0
     >                                )
                    Endif
* 
                    If (Dim .eq. 6) Then                ! (case d = 6)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf =  0.5 
     >                               * Cmplx     
     >                                (
     >                                  X**4 * EvenLog + 0.5d0 + X*X, 
     >                                  X**4 * Pi
C     >                                0.0D0
     >                                )
                    Endif
*     ---------------------------------------------
*             ____                            !(s-channel options end)

              Else
*             ~~~~
*             ---------
*             t-channel
*             ---------
                    EvenLog = Log ( 1.d0/X**2 + 1.d0)
                     OddLog = ATan (1.d0/X)
*
                    If (Dim .eq. 2) Then                ! (case d = 2)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf  = 0.5          
     >                               * Cmplx
     >                                (
     >                                 EvenLog,
     >                                 0.d0 
     >                                )
                    Endif
*                                                 
                    If (Dim .eq. 3) Then                ! (case d = 3)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf  =   1.d0 
     >                               * Cmplx
     >                                (
     >                                 1. - X * OddLog,
     >                                 0.d0
     >                                )
                    Endif
*                                                 
                    If (Dim .eq. 4) Then                ! (case d = 4)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf  =   0.5d0
     >                               * Cmplx
     >                                (
     >                                 1.d0 - X*X * EvenLog,
     >                                 0.d0
     >                                )
                    Endif
*                                                 
                    If (Dim .eq. 5) Then                ! (case d = 5)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf  = 1.d0
     >                               * Cmplx
     >                                (
     >                           1.d0/3.d0 - X*X + X**3 * OddLog,
     >                           0.d0
     >                                )
                    Endif
*                                                 
                    If (Dim .eq. 6) Then                ! (case d = 6)
*                   ~~~~~~~~~~~~~~~~~~~~                              
                        Lambdaf  = 0.5d0
     >                               * Cmplx
     >                                (
     >                                 0.5d0 - X*X + X**4*EvenLog,
     >                                 0.d0
     >                                )
                    Endif
*                                                                     
              Endif                          ! (t-channel options end)
*             ____                                                    
              Elseif (IModel.eq.1) THEN      !(ADD Model Options End)
*             ~~~~                             (RS Model)                       

*         -----------
*         RS Model 
*         -----------
*             Branch between channels:
*             ~~~~~~~~~~~~~~~~~~~~~
               m0=ams
               mh=amh
               c0=c00
              If (IChanl .eq. 0) Then
*             ---------
*             s-channel
*             ---------

* Initialize:
* ~~~~~~~~~~
                       BWigRe = 0.d0
                       BWigIm = 0.d0
                             NRes   = 0
 99                          NRes = NRes + 1   ! (increment res. no)
              IF(NRes. Eq. 1) Then
                              zJ1 = 3.83170605d0 ! (first zero of J1(x))
              Else
                              zJ1 = 7.01558685d0 + (NRes - 2) * Pi
              EndIf
              If (zJ1 .lt. X) Go To 99
                              NRestart = NRes

 88                        GW = GWidth (NRes, m0, c0, mH) ! (G_n width)
                       BWReTm = (X**2 - zJ1**2)
*                              ----------------------------------
     >                        /((X**2 - zJ1**2)**2 + (zJ1*GW)**2)
                       If (DAbs(BWReTm) .Gt. 1.D-2*DAbs(BWigRe)) Then
                               BWigRe = BWigRe + BWReTM !(real part)
                               KReal  = 1
                       Else
                               KReal  =  0
                       Endif
                       BWImTm = - zJ1 * GW
*                              ----------------------------------
     >                        /((X**2 - zJ1**2)**2 + (zJ1*GW)**2)
                       If (DAbs(BWImTm) .Gt. 1.D-2*DAbs(BWigIm)) Then
                               BWigIm = BWigIm + BWImTm !(real part)
                               KImag  = 1
                       Else
                               KImag  =  0
                       Endif
                       If ((KReal+KImag) .ne. 0) Then
                             NRes = NRes + 1   ! (increment res. no)
                             zJ1  = zJ1 + Pi   ! (next resonance)
                             Go To 88
                       Endif
                              NUppr = NRes
                              NRes = NRestart - 1
              If (Nres .Ge. 1) Then
              IF(NRes. Eq. 1) zJ1 = 3.83170605d0 ! (first zero of J1(x))
                              zJ1 = 7.01558685d0 + (NRes - 2) * Pi
 77                        GW = GWidth (NRes, m0, c0, mH) ! (G_n width)
                       BWReTm = (X**2 - zJ1**2)
*                              ----------------------------------
     >                        /((X**2 - zJ1**2)**2 + (zJ1*GW)**2)
                       If (DAbs(BWReTm) .Gt. 1.D-2*DAbs(BWigRe)) Then
                               BWigRe = BWigRe + BWReTM !(real part)
                               KReal  = 1
                       Else
                               KReal  =  0
                       Endif
                       BWImTm = - zJ1 * GW
*                              ----------------------------------
     >                        /((X**2 - zJ1**2)**2 + (zJ1*GW)**2)
                       If (DAbs(BWImTm) .Gt. 1.D-2*DAbs(BWigIm)) Then
                               BWigIm = BWigIm + BWImTM !(real part)
                               KImag  = 1
                       Else
                               KImag  =  0
                       Endif
                       If ((KReal+KImag) .ne. 0) Then
                             NRes = NRes - 1   ! (increment res. no)
                             If (NRes .Gt. 2) zJ1  = zJ1 - Pi
                             If (NRes. Eq. 1) zJ1 = 3.83170605d0
                             If (NRes. Lt. 1) Go To 66
                             Go To 77
                       Endif
              EndIf
 66                           NLowr = NRes

                         Arg1 = (7.01558685d0 + X) /Pi + NUppr - 2
                         Arg2 = (7.01558685d0 - X) /Pi + NUppr - 2
                       BWigRe = BWigRe
     >                        + ( DiGamma(Arg1) - DiGamma(Arg2))
*                              ---------------------------------
     >                          /(2.d0 * Pi * X)

              If (NLowr .Ne. 0) Then
                         Arg1 = (7.01558685d0 + X) /Pi - 1
                         Arg2 = (7.01558685d0 - X) /Pi - 1
                       BWigRe = BWigRe
     >                        + ( DiGamma(Arg1) - DiGamma(Arg2))
*                              ---------------------------------
     >                          /(2.d0 * Pi * X)

                         Arg1 = (7.01558685d0 - X) /Pi + NLowr - 1
                         Arg2 = (7.01558685d0 + X) /Pi + NLowr - 1
                       BWigRe = BWigRe
     >                        + ( DiGamma(Arg1) - DiGamma(Arg2))
*                              ---------------------------------
     >                          /(2.d0 * Pi * X)
              EndIf

c          if (Isep.eq.0) then
c          Lambdaf = Cmplx(BWigRe,BWigIm)
c          ElseIf (Isep .eq. 1) then
c          Lambdaf = Cmplx(BWigRe,0d0)
c          endif

          Lambdaf = Cmplx(BWigRe,BWigIm)

******************************************************************************
* This Isep is included in the RS case only on 05-04-2008.
******************************************************************************

*             ____
              Else
*             ~~~~     
*             ---------
*             t-channel
*             ---------
*
* (The function is calculated in Mathematica and here we just fit it)
*
* range x > 100 :
* ~~~~~~~~~~~~~~~
          If ( x .Gt. 100.d0) Then            
                            Lambdaf = Cmplx(0.d0,0.d0)
                            Return
          EndIf
* range x = 10 to 100:
* ~~~~~~~~~~~~~~~~~~~
          If ( x .Ge. 10.d0) Then    
                           Nfit = x + 100
                           FitY = ((Yfit(Nfit+1) - Yfit(Nfit) )*x
     >                           +( Yfit(Nfit)   * Xfit(Nfit+1)
     >                             -Yfit(Nfit+1) * Xfit(Nfit) )
     >                            )/ (Xfit(Nfit+1) - Xfit(Nfit) )
                           Lambdaf = Cmplx(-FitY,0.d0)/ (32.*Pi)
          EndIf 
* range x = 0.1 to 10:
* ~~~~~~~~~~~~~~~~~~~
          If ( x .Lt. 10.d0 .And. x .Ge. 0.1d0) Then
                           Nfit = 10*x + 10
                           FitY = ((Yfit(Nfit+1) - Yfit(Nfit) )*x
     >                           +( Yfit(Nfit)   * Xfit(Nfit+1)
     >                             -Yfit(Nfit+1) * Xfit(Nfit) )
     >                            )/ (Xfit(Nfit+1) - Xfit(Nfit) )
                           Lambdaf = Cmplx(-FitY,0.d0)/ (32.*Pi)
          EndIf
* range x = 0.01 to 0.1:
* ~~~~~~~~~~~~~~~~~~~~~~
          If ( x .Lt. 0.1d0 .And. x .Ge. 0.01d0) Then 
                           Nfit = 100*x + 1
                           FitY = ((Yfit(Nfit+1) - Yfit(Nfit) )*x
     >                           +( Yfit(Nfit)   * Xfit(Nfit+1)
     >                             -Yfit(Nfit+1) * Xfit(Nfit) )
     >                            )/ (Xfit(Nfit+1) - Xfit(Nfit) )
                           Lambdaf = Cmplx(-FitY,0.d0)/ (32.*Pi)
          EndIf
* range x < 0.01:
* ~~~~~~~~~~~~~~~
          If ( x. Lt. 0.01d0) Then    
                           Lambdaf = Cmplx(-12.5768d0,0.d0)/ (32.*Pi)
          EndIf
* range x < 0:
* ~~~~~~~~~~~~
          If (x. Lt. 0.d0) Then
                            Lambdaf = Cmplx(0.d0,0.d0)/ (32.*Pi)
                            Return
          EndIf
*             _____
              Endif
*             ~~~~~
*             _____                                                   
              Endif                           ! (RS Model Options End)
*             ~~~~~                                                     
*--------------------------------------------------------------------*
          Return
          End
**********************************************************************
************************|-------------------|*************************
************************| Program: GWidth.f |*************************
************************|-------------------|*************************
**********************************************************************
*****************|-------------------------------|********************
*****************| Last modified:  June 1, 2000  |********************
*****************|-------------------------------|********************
**********************************************************************
********|----------------------------------------------------|********
********| This function calculates the width of Kaluza-Klein |********
********| graviton resonances in the RS model. It actually   |********
********| returns the width scaled to the parameter m0.      |********
********|----------------------------------------------------|********
**********************************************************************
*--------------------------------------------------------------------*
          Function GWidth (n, m0, c0, mH)
*--------------------------------------------------------------------*
*  Declarations:
*  ~~~~~~~~~~~~
          Implicit None
          Integer  n, Il, Iq
          Real*8   m0, c0, MW, MZ, mH, ml(3), mq(6)
          Real*8   zJ1,Pi, rW, rZ, rH, rl(3), rq(6)
          Real*8   GWidth, Delta, Mn
          Data Pi/3.141592653589793238462643D0/
          Data mW/80.4d0/, mZ/91.12d0/ 
          Data ml/0.000511d0, 0.105d0, 1.777d0/
          Data mq/0.005d0, 0.009d0, 1.3d0, 0.14d0, 173.8d0, 4.3d0/
*--------------------------------------------------------------------*
* Construct Graviton Mass (using zeros of Bessel function J1(x)):
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          If (n .eq. 1) zJ1 = 3.83170605d0                 ! (first zero)
          If (n .Gt. 1) zJ1 = 7.01558685d0 + (n - 2) * Pi  ! (higher zeros)
                        Mn = zJ1 * m0                      ! (graviton mass)
* Construct Mass Ratios:
* ~~~~~~~~~~~~~~~~~~~~~
             rW = (mW /Mn)**2                     ! (W boson)      
             rZ = (mZ /Mn)**2                     ! (Z boson)      
             rH = (mH /Mn)**2                     ! (Higgs boson)  
          Do Il = 1,3
             rl(Il) = (ml(Il) /Mn)**2             ! (charged leptons)
          EndDo
          Do Iq = 1,6
             rq(Iq) = (mq(Iq) /Mn)**2             ! (quarks)
          EndDo
*--------------------------------------------------------------------*
* Construct Width:
* ~~~~~~~~~~~~~~~~
*  contribution from massless final states (always present):
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             Delta =    0.2d0           !   (photon contribution)
     >                + 1.6d0           !    (gluon contribution)
     >                + 0.3d0           ! (neutrino contribution)
*  W pair contribution:
*  ~~~~~~~~~~~~~~~~~~~
             If (Mn .Gt. 2*mW) Then     ! (kinematic check)
             Delta = Delta 
c     >                + 0.4d0 * Sqrt( 1.d0 - 4d0*rW**2) 
     >                + 0.4d0 * Sqrt( 1.d0 - 4d0*rW) 
     >                  *(13.d0/12 + 14.d0*rW/39 + 4.d0*rW**2/13)
             EndIf
*  Z pair contribution:
*  ~~~~~~~~~~~~~~~~~~~
             If (Mn .Gt. 2*mZ) Then     ! (kinematic check)
             Delta = Delta 
c     >                + 0.2d0 * Sqrt( 1.d0 - 4d0*rZ**2) 
     >                + 0.2d0 * Sqrt( 1.d0 - 4d0*rZ) 
     >                  *(13.d0/12 + 14.d0*rZ/39 + 4.d0*rZ**2/13)
             EndIf
*  Higgs pair contribution:
*  ~~~~~~~~~~~~~~~~~~~~~~~
             If (Mn .Gt. 2*mH) Then     ! (kinematic check)
             Delta = Delta 
     >                + (1.d0/30)* Sqrt( (1.d0 - 4d0*rH)**5 )
             EndIf
*  charged lepton pair contribution:
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Do Il = 1,3
             If (Mn .Gt. 2*ml(Il)) Then     ! (kinematic check)
             Delta = Delta 
     >                + 0.1d0 * Sqrt( (1. - 4d0*rl(Il))**3 ) 
     >                      *( 1.d0 + 8.d0*rl(Il) /3d0 )
             EndIf
          EndDo
*  quark pair contribution:
*  ~~~~~~~~~~~~~~~~~~~~~~~
          Do Iq = 1,6
             If (Mn .Gt. 2*mq(Iq)) Then     ! (kinematic check)
             Delta = Delta 
     >                + 0.3d0 * Sqrt( (1.d0 - 4d0*rq(Iq))**3 ) 
     >                      *( 1.d0 + 8.d0*rq(Iq) /3d0 )
             EndIf
          EndDo
*--------------------------------------------------------------------*
* Include Overall Factor:
* ~~~~~~~~~~~~~~~~~~~~~~
          GWidth = Delta * c0**2 * zJ1**3 
          GWidth = GWidth/2.d0 ! New Convention
*--------------------------------------------------------------------*
* Function SubProgram Ends:
* ~~~~~~~~~~~~~~~~~~~~~~~~
          Return
          End
**********************************************************************
************************|--------------------|************************
************************| Program: DiGamma.f |************************
************************|--------------------|************************
**********************************************************************
*****************|-------------------------------|********************
*****************| Last modified:  June 1, 2000  |********************
*****************|-------------------------------|********************
**********************************************************************
********|----------------------------------------------------|********
********| This function calculates the Euler Digamma function|********
********| for a real argument using the series expansion:    |********
********|                                                    |********
********|  \psi(x) = - \gamma                                |********
********|             + \sum_{n=1}^{\infty} (x-1)/( n(n+x-1) |********
********|                                                    |********
********| given by Abramowitz and Stegun.                    |********
********| \gamma is Euler's constant = 0.57721566            |********
********|----------------------------------------------------|********
********| Warning: This does not converge for x = 0,-1,-2,.. |********
********|----------------------------------------------------|********
**********************************************************************
*--------------------------------------------------------------------*
          Function DiGamma (X)
*--------------------------------------------------------------------*
*  Declarations:
*  ~~~~~~~~~~~~
          Implicit None
          Integer n
          Real*8  X, DiGamma, EulConst, Term, Eps
          Data EulConst/0.57721566d0/
*--------------------------------------------------------------------*
             n = 0
                             DiGamma = - EulConst
 99          n = n +1
          Term = (X - 1.d0) /n/(n - 1.d0 + X)
          Eps  = DAbs(Term/DiGamma)
          If( Eps. Gt. 1.D-4) Then
                             DiGamma = DiGamma + Term
                             Go To 99
          EndIf
*--------------------------------------------------------------------*
          Return
          End
**********************************************************************


c >> >> >> >> >> >> > m1=x1 m0
c >> >> >> >> >> >> > x1=3.83
c >> >> >> >> >> >> >
c >> >> >> >> >> >> > so I had chosen a value of m1 about 1.5 TeV for LHC and
c >> >> about
c >> >> >> a
c >> >> >> >> >> 1TeV
c >> >> >> >> >> >> for TEV.
c >> >> >> >> >> >> >
c >> >> >> >> >> >> > 0.01 < c0 =k/Mp < .1  or so so that there is no further
c >> >> >> >> hierarcy.
*---------------------------------------------------------------------------
*           gamma function
*---------------------------------------------------------------------------

        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute the gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú )
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

       double precision function  SWidth(Ms, Lam)
       implicit double precision (a-h, o-z)
       double precision Ms, Lam
       double precision lambda
       data pi/3.141592653589793238462643D0/
       common/param/aem,xmur,lambda
       common/sparam/smass,xlam
       common/scalar_couplings/cph,cgl,cfe
       common/sm_couplings/c1,c2,c3,ct,g1,g2,g3
       common/weak/amz,amw,sw2,cw2
       common/amfermion/amtau,amc,amb,amt


c       write(*,*)'c1,c2,c3,ct =',c1,c2,c3,ct
c       write(*,*)'g1,g2,g3 =',g1,g2,g3
c       write(*,*)'Mz, Mw, Mt, Sw2, Cw2 =', amz, amw, amt, sw2,cw2
c       write(*,*)'Scalar mass, Lambda =', Ms, Lam

       DWidth = 0.0d0
c -----------------------------------------------
c  Scalar decays to a pair of photons
c -----------------------------------------------
       e = dsqrt(4.0d0*pi*aem)
       Tphph1 = (c1+c2)**2 * e**4 * Ms**3
       Tphph2 = 4.0d0*pi*Lam**2
       DWphph   = Tphph1/Tphph2

       DWidth = DWidth + DWphph

c -----------------------------------------------
c  Scalar decays to a pair of gluons
c -----------------------------------------------
       Tgg1 = 2.0d0 * c3**2 * g3**4 *Ms**3
       Tgg2 = pi*Lam**2
       DWgg = Tgg1/Tgg2 

       DWidth = DWidth + DWgg

c -----------------------------------------------
c  Scalar decays to photon and a Z boson
c -----------------------------------------------
       if (Ms .gt. amz) then
       d1 = c1*g1**2 - c2*g2**2
       d2 = Ms**3 - amz**3
       Tzph1 = sw2*cw2 * d1**2 * d2**2
       Tzph2 = 2.0d0*pi* Ms**3 * Lam**2
       DWzph = Tzph1/Tzph2
       endif

       DWidth = DWidth + DWzph

c -----------------------------------------------
c  Scalar decays to a pair of W bosons
c -----------------------------------------------
       if (Ms .gt. 2.0d0*amw) then
       a1 = Ms**2 - 4.0d0*amw**2
       a2 = Ms**4 - 4.0d0*Ms**2*amw**2 + 6.0d0*amw**4
       Tww1 = c2**2 * g2**4*dsqrt(a1)*a2
       Tww2 = 2.0d0*pi*Ms**2*Lam**2
       DWww = Tww1/Tww2
       endif

       DWidth = DWidth + DWww
c -----------------------------------------------
c  Scalar decays to a pair of Z bosons
c -----------------------------------------------
       if (Ms. gt. 2.0d0*amz) then
       b1 = c2*g2**2*cw2 + c1*g1**2*sw2
       b2 = Ms**2 - 4.0d0*amz**2
       b3 = Ms**4 - 4.0d0*Ms**2*amz**2 + 6.0d0*amz**4
       Tzz1 = b1**2*dsqrt(b2)*b3
       Tzz2 = 4.0d0*pi*Ms**2*Lam**2
       DWzz = Tzz1/Tzz2
       endif

       DWidth = DWidth + DWzz
c -----------------------------------------------
c  Scalar decays to a pair of top quarks (t tbar)
c -----------------------------------------------


       DWchm = fbr(ct, amc, Lam, Ms)
       DWbot = fbr(ct, amb, Lam, Ms)
       DWttb = fbr(ct, amt, Lam, Ms)

       DWidth = DWidth + DWchm
       DWidth = DWidth + DWbot
       DWidth = DWidth + DWttb

c       write(*,*)'Scalar Decay Width, Ms, Lambda =', DWidth, Ms, Lam
       SWidth = DWidth

       return
       end

      real *8 function fbr(cmf, amf, Lam, Ms)
      implicit double precision (a-h, o-z)
      double precision Ms, Lam
      data pi/3.141592653589793238462643D0/

       if (Ms .gt. 2.0d0*amf) then
       aa = Ms**2 - 4.0d0*amf**2
       Tttb1 = 3.0d0*cmf**2 * amf**2 * aa*dsqrt(aa)
       Tttb2 = 8.0d0*pi*Ms**2*Lam**2
       fbr = Tttb1/Tttb2
       else
       fbr = 0.0d0
       endif

       return
       end


