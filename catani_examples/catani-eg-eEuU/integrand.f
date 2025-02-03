c      integer function integrand(ndim,x,ncomp,f,userdata,nvec ,core)
      function fnlo3(xx,weight)
      implicit double precision (a-h,o-z)
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      integer i35,i45,is5,itest
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
c      parameter (hbarc2 =0.389379d0*10**9)
      common/cmenergy/s
      ge=0.007547169811320755d0
      e = DSQRT(ge*4.d0*PI)
      AL=0.118d0
      Cf=4d0/3d0
c      xM2=(4d0*e**4)/(16d0*pi*3d0)             ! Leading Order msq
      PF= 1d0/(16d0*pi**2)                      ! phase space factor
      x1=xx(1)
      x2=xx(2)                                  ! x2 runs from 0 to 1-x1 used theta(x1+x2-1)

       qu2 = 4d0/9d0                            !
c      xM2 = (3d0*ge**2*qu2)/(4d0*s)            ! normalized msq 
       xM2 = 4d0*PI*ge**2*qu2/s                 ! From the RD Field expression
       PS_norm  = 1d0
       xM2_norm = xM2*PS_norm                   !

      if ( x1+x2 .ge. 1d0 ) then
c
          sig = (x1**2+x2**2)/((1-x1)*(1-x2))
      dipole1 = (1d0/(1-x2)*(2d0/(2d0-x1-x2)-(1+x1))+(1-x1)/x2)
      dipole2 = (1d0/(1-x1)*(2d0/(2d0-x1-x2)-(1+x2))+(1-x2)/x1)
        xnorm = hbarc2*xM2_norm*Cf*AL/2d0/PI

        fnlo3 = xnorm*(sig-dipole1-dipole2)
c        fnlo3 = xnorm*(dipole1+dipole2)

      endif
      return
      end
