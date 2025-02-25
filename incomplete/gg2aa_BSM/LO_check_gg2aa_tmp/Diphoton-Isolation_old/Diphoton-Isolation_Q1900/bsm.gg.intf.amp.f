* Function programmes containing the helicity amplitudes of gg box
* diagram.

      function bsmsigma_gg_intf(s,t,u)
      implicit double precision (a-h,o-z)
      complex *8 cmplx_ak2d
      complex *8 smhamp2p2m,smfstu,smfsut
      common/param/aem,xmur,lambda
      common/model/model
      common/nflavour/nf
      complex*8 lambdaf
      common/xpar/c00,ams,amh
      common/add_par/xms,nd
      common/rs_par/aam1,c0,aamh
      common/unpar/xl3,xdu,xlamu
      data Pi/3.141592653589793238462643D0/

        xnf=nf
        alphs=alfas2(xmur,lambda,xnf)
        as   =alphs/4.d0/pi

        xmh=dsqrt(s)
        xmh2=xmh*xmh
        XMS=xms
        ND=nd

        c00=c0
        amh=aamh
        am1=aam1
        am0=am1/3.83D0

       XDU=xdu
       MODL=MODEL

c       write(*,*)'model=',MODEL
       IF (MODL.EQ.3) GO TO 300
       IF (MODL.EQ.2) GO TO 200
       IF (MODL.EQ.1) GO TO 100
       if (model.eq.0) then
       bsmsigma_gg_intf=0d0
       go to 199
       endif

  100  ams=XMS
       rms=xmh/ams
       cmplx_ak2d = -lambdaf(0,0,ND,rms)
       cmplx_ak2d = 16.0D0*PI/ams/ams/ams/ams*cmplx_ak2d
c       write(*,*)xmh,ams,nd
c       write(*,*)'aem,pi',aem,PI
c       write(*,*)'M_s,nd,rms,Cmplx(ak2d)',ams,ND,rms,bsm_cf
       GO TO 110

  200  ams=am0
       rms=xmh/ams
       cmplx_ak2d=lambdaf(1,0,ND,rms)
       cmplx_ak2d = 16.0D0*PI*C00**2/ams/ams/ams/ams*cmplx_ak2d
c       write(*,*)'am1,c0,amh,nd,rms',am1,c0,amh,ND,rms
       GO TO 110

  300  xdu1=2.D0-xdu
       xdu2=xdu+2.0D0
       xdum1=xdu-1.d0
       CALL gamma(xdu1,gamma1) 
       CALL gamma(xdu2,gamma2) 
c      ct=C_T, the tensor coupling,  given in Grinstein et. al.
       ct=1.d0
       delu=ct*gamma1/gamma2/4.d0**xdum1*xdu*xdum1
       XLAMU2=XLAMU*XLAMU
       FGRU=-(2.0D0*XL3)**2*delu*(1.0D0/2.0D0)*
     1 (XMH2/XLAMU2)**(xdu)/XMH2/XMH2
c       write(*,*)'du,ct,Lambda_U,fgru=',xdu,XL3,xlamu,FGRU
c       write(*,*)'G(2-d),G(d+2),d-1=',gamma1,gamma2,xdum1
c       write(*,*)'Q,Q^2=',XMH,XMH2

C  (-1)^du in the interference part gives cos(du\pi) - i sin(du\pi)

       FGRU1=FGRU*(DCOS(XDU*PI))
       FGRU2=FGRU*(-DSIN(XDU*PI))
       cmplx_ak2d=cmplx(FGRU1,FGRU2)

       GO TO 110

  110   continue 

      upq=2d0/3d0
      dnq=-1d0/3d0

      Qqnf5=2*upq**2+3*dnq**2
      Qq2nf5=Qqnf5*Qqnf5

c     smfstu=M+-+-(s,t,u) without the SM couplings put-in
c     smfsut=M+--+(s,t,u) without the SM couplings put-in

      smfstu=smhamp2p2m(s,t,u)
      smfsut=smhamp2p2m(s,u,t)

      Re_cf1=2.d0*real((conjg(smfstu))*cmplx_ak2d)
      Re_cf2=2.d0*real((conjg(smfsut))*cmplx_ak2d)

c     bsm1=M+-+-(s,t,u) for the BSM gg --> \gamma \gamm without the
c     complex k^2D(s) factor.
c     bsm2=M+--+(s,t,u) for the BSM gg --> \gamma \gamm without the
c     complex k^2D(s) factor.

      bsm1=0.5d0*u**2
      bsm2=0.5d0*t**2

      sm_cf=4d0*aem*alphs*Qqnf5
c     \delta^ab*\delta_ab=8d0
      sm_cf=sm_cf*8.d0
      am12=sm_cf*Re_cf1*bsm1
      am22=sm_cf*Re_cf2*bsm2

c     The multiplicity of each kind of the helicity amplitude is 2
      bsmsigma_intf=2d0*am12+2d0*am22

c  Helcity, color average for each gluon and 
c  statistical average= 1/4.1/8.1/8.1/2
      bsmsigma_gg_intf=bsmsigma_intf/4d0/8d0/8d0/2d0
 199    return
        end


