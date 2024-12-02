* Function programmes containing the helicity amplitudes of gg box
* diagram.

      function smsigma_gg(s,t,u)
      implicit double precision (a-h,o-z)
      complex *8 smhamp2p2m
      common/param/aem,xmur,lambda
      common/nflavour/nf
      data Pi/3.141592653589793238462643D0/

      xnf=nf
      alphs=alfas2(xmur,lambda,xnf)
      as   =alphs/4.d0/pi

      upq=2d0/3d0
      dnq=-1d0/3d0

      Qqnf5=2*upq**2+3*dnq**2
      Qq2nf5=Qqnf5*Qqnf5

      capx=smhamp4p(s,t,u)
      capy=cabs(smhamp2p2m(s,t,u))
      capz=cabs(smhamp2p2m(s,u,t))

      ccasq=4d0**2*8d0*aem**2*alphs**2*Qq2nf5
      am12=ccasq*capx**2
      am22=ccasq
      am32=ccasq
      am42=ccasq*capy**2
      am52=ccasq*capz**2

      smsigma=2d0*am12+2d0*am22+8d0*am32+2d0*am42+2d0*am52

c  Helcity, color average for each gluon and 
c statistical average= 1/4.1/8.1/8.1/2
      smsigma_gg=smsigma/4d0/8d0/8d0/2d0
      return
      end

c M^++++_SM (s,t,u) (SM helicity amplitude for ++++)

      function smhamp4p(s,t,u)
      implicit double precision (a-h,o-z)
      data Pi/3.141592653589793238462643D0/

      ubt = dabs(u/t)
      dlnut=dlog(ubt)
      x_stu=1.d0+(u-t)/s*dlnut+0.5d0*(t**2+u**2)/s/s*(dlnut**2+Pi**2)
      smhamp4p=x_stu 
      return
      end


c M^+-++-_SM (s,t,u) (SM helicity amplitude for +-+-)

      function smhamp2p2m(s,t,u)
      implicit double precision (a-h,o-z)
      complex *8 smhamp2p2m
      data Pi/3.141592653589793238462643D0/

      tbs = dabs(t/s)
      dlnts=dlog(tbs)

      x_stu=1.d0+(t-s)/u*dlnts+0.5d0*(s**2+t**2)/u/u*(dlnts**2)
      y_stu=  (t-s)/u+(s**2+t**2)/u/u*dlnts

      RePart=x_stu*(1.d0)
      AImpart=y_stu*Pi*(1.d0)

      smhamp2p2m=cmplx(Repart,AImpart)

      return
      end


