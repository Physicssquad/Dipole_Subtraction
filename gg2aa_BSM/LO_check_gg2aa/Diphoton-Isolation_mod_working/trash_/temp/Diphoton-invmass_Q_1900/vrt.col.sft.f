c     initial state soft collinear and final state hard collinear terms
      subroutine setcol(s12,scale,deltas,deltac,col)
      implicit double precision (a-h,o-z)
      dimension col(10)
      parameter (cf=4d0/3d0)
      parameter (ca=3.0d0)
      parameter (zeta2=1.64493406684823d0)
      common/corr/icorr
      common/nflavour/nf

      xlds=dlog(deltas)
      xldc=dlog(deltac)
      
      xls=dlog(s12/scale**2)

!qqb sm
      c_i_q=cf*(4.0d0*xlds+3.0d0)*xls

      if(icorr.eq.2)then
         dr=dilog(deltac/deltas)
         c_i_q=c_i_q+4.0d0*cf*dr
      endif
      col(1)=2.d0*c_i_q

!gg sm
      col(2)=0.d0

! qqb bsm
      c_i_q=cf*(4*xlds+3d0)*xls
      
      if(icorr.eq.2)then
         dr=dilog(deltac/deltas)
         c_i_q=c_i_q+4*cf*dr
      endif
      col(3)=2.d0*c_i_q

! gg bsm
      c_i_g=ca*(4*xlds+11.d0/3d0)*xls
      c_i_g=c_i_g - 2.d0/3.d0*nf*xls

      if(icorr.eq.2)then
         dr=dilog(deltac/deltas)
         c_i_g=c_i_g+4*ca*dr  ! check this
      endif
      col(4)=2.d0*c_i_g
!qqb int
      c_i_q=cf*(4*xlds+3d0)*xls
      
      if(icorr.eq.2)then
         dr=dilog(deltac/deltas)
         c_i_q=c_i_q+4*cf*dr
      endif
      col(5)=2.d0*c_i_q
!gg int
      col(6)=0.d0


      return
      end

c     virtual calls
      subroutine setvrt(vrt)
      implicit double precision (a-h,o-z)
      dimension vrt(10)
      parameter (cf=4d0/3d0)
      parameter (zeta2=1.64493406684823d0)
      parameter (pi=3.14159265358979d0)
      parameter (N=3)
      parameter (Tf=0.5D0)
      common/nflavour/nf
      common/xmandelstam/s,t,u
      common/model/model

      ALEM=1.0D0/128.0D0
      e=DSQRT(4.0D0*PI*ALEM)

      xlus=dlog(-u/s)
      xlts=dlog(-t/s)
      xlus2=xlus**2
      xlts2=xlts**2

!born is not factored out unlike Drell-Yan hence 
!watch sig.nlo2.sm.qqb.f 

      vrt(1)=-14d0*(u*u + t*t)/(u*t)
     &     +zeta2*8d0*(u*u + t*t)/(u*t)
     &     +xlus*(4d0 + 6d0*t/u)
     &     +xlts*(4d0 + 6d0*u/t)
     &     +xlus2*(4d0+ (4d0*u*u + 2d0*t*t)/(u*t))
     &     +xlts2*(4d0+ (4d0*t*t + 2d0*u*u)/(u*t))

      vrt(1)=e**4*vrt(1)*cf/n


      vrt(2)=0.0d0
      vrt(3)=cf*(-20.0d0 +8.0d0*zeta2)
      vrt(4) = ( -N*203.d0/9.d0
     &            +nf*Tf*70.d0/9.d0
     &            +N*zeta2*8.0d0 )



      rp34 = dsqrt(s)

      call coupfact(model,rp34,AK2D,AK2DINTF)

!born is not factored out unlike Drell-Yan hence 
!watch sig.nlo2.bsm.intf.qqb.f 
      vrt(5)= (2.d0*u**2/s +u -t)*17.d0/4.d0
     &      +  xlus*(-u**2/s -u +3.0d0*t)/4.0d0
     &      + xlus2*(-u**2/s +u +t)/4.0d0
     &      +  xlts*(-u**2/s +2.0d0*u)/4.0d0
     &      + xlts2*(-u**2/s+2.0d0*t)/4.0d0
     &      + zeta2*(-4.0d0*u**2/s-2.0d0*u+2.0d0*t)

      vrt(5)=vrt(5)*s*AK2DINTF*cf/N     !ak2d =kap^2/s
      vrt(5)=2.0d0*vrt(5)           !interfence has factor 2
      vrt(6)=0.0d0

      return
      end

c. soft calls
      subroutine setsft(deltas,sft)
      implicit double precision (a-h,o-z)
      dimension sft(10)
      parameter (cf=4d0/3d0)
      parameter (ca=3.0d0)

      xlds=dlog(deltas)
      sft(1)=8.0d0*cf*xlds*xlds
      sft(2)=0.0d0
      sft(3)=8.d0*cf*xlds*xlds
      sft(4)=8.d0*ca*xlds*xlds
      sft(5)=8.d0*cf*xlds*xlds
      sft(6)=0.d0

      return
      end









