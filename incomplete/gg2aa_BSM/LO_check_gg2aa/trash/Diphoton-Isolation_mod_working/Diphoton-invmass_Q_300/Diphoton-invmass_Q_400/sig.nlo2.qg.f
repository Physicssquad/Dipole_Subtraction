c.
c. 2 -> 2 nlo cross section
c.
      subroutine sig_nlo2_qg(scale,f1,f2,f1tqg,f2tqg,f1tgq,f2tgq,
     &     p1,p2,p3,p4,tot)
      implicit double precision (a-h,o-z)
      real*8 lambda
      dimension f1(-6:6),f2(-6:6)
      dimension f1tqg(-6:6),f2tqg(-6:6)
      dimension f1tgq(-6:6),f2tgq(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension col(10),sft(10),vrt(10),sig(10)
      dimension xl1qg(15),xl2qg(15),xl1gq(15),xl2gq(15)
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

      call setlum(f1tqg,f2,xl1qg)
      call setlum(f1,f2tqg,xl2qg)
      call setlum(f1tgq,f2,xl1gq)
      call setlum(f1,f2tgq,xl2gq)

      call setsigma(s,t,u,sig)
      
      if (is.eq.2) then                 ! qg

        if (io.eq.1) then               ! SM
        tot=0.d0
        tot=tot+2.d0*(xl1qg(1)+xl2qg(1))*sig(1)

        elseif (io.eq.2) then           ! BSM direct
        tot=0.d0
        tot=tot+2.d0*(xl1qg(3)+xl2qg(3))*sig(3)
        tot=tot+2.d0*(xl1gq(4)+xl2gq(4))*sig(4)
      
        elseif (io.eq.3) then           ! SM*BSM interference
        tot =0.d0
        tot=tot+2.d0*(xl1qg(5)+xl2qg(5))*sig(5)  
        endif

      tot=tot*as
      endif

      return
      end


