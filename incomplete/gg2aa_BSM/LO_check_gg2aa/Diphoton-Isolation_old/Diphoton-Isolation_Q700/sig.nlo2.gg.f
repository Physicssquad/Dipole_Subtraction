c.
c. 2 -> 2 nlo cross section
c.
      subroutine sig_nlo2_gg(scale,f1,f2,f1tgg,f2tgg,
     .                           p1,p2,p3,p4,tot)
      implicit double precision (a-h,o-z)
      real*8 lambda
      dimension f1(-6:6),f2(-6:6)
      dimension f1tgg(-6:6),f2tgg(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension col(10),sft(10),vrt(10),sig(10)
      dimension xl(15),xl1gg(15),xl2gg(15)
      parameter (pi=3.14159265358979d0)

      common/slice/deltas,deltac
      common/xmandelstam/s,t,u
      common/param/aem,xmur,lambda
      common/factscale/xmuf
      common/nflavour/nf
      common/isub/io,is
     
      s=2d0*dot(p1,p2)
      t=-2d0*dot(p1,p3)
      u=-(s+t)

      xnf=nf
          
      alphs=alfas2(xmur,lambda,xnf)
      as   =alphs/4.d0/pi
      
      call setlum(f1,f2,xl)
      call setlum(f1tgg,f2,xl1gg)
      call setlum(f1,f2tgg,xl2gg)

      call setcol(s,xmuf,deltas,deltac,col)
      call setsft(deltas,sft)
      call setvrt(vrt)

      call setsigma(s,t,u,sig)

      if (is.eq.3) then                 ! gg

        if (io.eq.2) then               ! BSM direct
        tot= xl(4)*sig(4)*(
     &          vrt(4)+
     &          sft(4)+
     &          col(4) )
        tot=tot+2*(xl1gg(4)+xl2gg(4))*sig(4)
      endif

      tot=tot*as
      endif

c     ##### 2*(xl1(4)+xl2(4))*sig(4): 2 because alphas/2/Pi = 2*as
c      tot=(xl1(4)+xl2(4))*sig(4) ! for testing compton diagrams

      return
      end
