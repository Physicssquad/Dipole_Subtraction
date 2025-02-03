      program eE2uUg_dipoleSubtraction
c       Catani formalism using examples
      implicit double precision(a-h,o-z)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.389379d9)
      common/amass/am1,am2,am3,am4,am5
      common/cmenergy/s
      external fnlo3
      
            !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')
      read (10,*) npt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      close(10)

      open(unit=15,file='run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      close(15)
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=0d0
      am4=0d0
      am5=0d0
      ! energy
      s=ecm*ecm
      
      ge=0.007547169811320755d0
c      print*,"LO RD field :",4d0*pi*ge**2/3d0/10**6
      e= DSQRT(ge*4.d0*PI)

        qu2 = 4d0/9d0
        aM2_norm =  4d0*PI*ge**2*qu2/s        
        Cf = 4d0/3d0
        AL = 0.118d0
        sig_nlo_th = 4d0*CF*AL*ge**2/3d0/s 
        print*,"sig_nlo_th: ",sig_nlo_th

        call brm48i(40,0,0) ! initialize random number generator
        call vsup(3,npt1,its1,fnlo3,ai_nlo3,sd,chi2)

c        write(*,*)'The answer is =', ai_nlo3
c        write(*,*)"Integral      =",ai_nlo3,"+-",sd
c        write(*,*)"with chisq    =",chi2
c        ai_nlo3=1.4227985845555067E-002
c        aM2= (-3.0030847159413110E-004)

        sig1_r = ai_nlo3
        sig2_v = Cf*AL*aM2_norm/pi*hbarc2

        sig_nlo= sig1_r + sig2_v
        
        sig_nlo_th = 4d0*CF*AL*ge**2/3d0/s*hbarc2 

        print*,"                   Sig(NLO) : ",sig_nlo,"  pb"
        print*,"analytical value sigma(NLO) : ",sig_nlo_th," pb"
        print*," "
        print*,"Sigma_r + Sigma_v:",sig1_r," + ",sig2_v   ," pb"
c        print*,"0.389379 pb =  1×10^(−9) GeV−2"
        print*," "
c        sig_nlo=sig_nlo*0.389379d0*10d9
        print*,"                   Sig(NLO) : ",sig_nlo,"  pb "

        
        end

