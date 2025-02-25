c.
c. 2 -> 2 nlo cross section
c.
      subroutine sig_nlo2_qqb(scale,f1,f2,f1tqq,f2tqq,
     .       p1,p2,p3,p4,tot)
      implicit double precision (a-h,o-z)
      real*8 lambda
      dimension f1(-6:6),f2(-6:6)
      dimension f1tqq(-6:6),f2tqq(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension col(10),sft(10),vrt(10),sig(10)
      dimension xl(15),xl1qq(15),xl2qq(15)
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
      call setlum(f1tqq,f2,xl1qq)
      call setlum(f1,f2tqq,xl2qq)

      call setcol(s,xmuf,deltas,deltac,col)
      call setsft(deltas,sft)
      call setvrt(vrt)

      call setsigma(s,t,u,sig)

        if (is.eq.1) then               ! qqb 

        if (io.eq.1) then               ! SM 
        tot = 0d0
        tot= xl(1)*sig(1)*(
     &          sft(1)+
     &          col(1) )
        tot=tot+xl(1)*vrt(1)  
        tot=tot+2.d0*(xl1qq(1)+xl2qq(1))*sig(1)  

        elseif (io.eq.2) then           ! BSM direct
        tot = 0d0
        tot= xl(3)*sig(3)*(
     &          vrt(3)+
     &          sft(3)+
     &          col(3) )
        tot=tot+2.0d0*(xl1qq(3)+xl2qq(3))*sig(3)  

        elseif (io.eq.3) then           ! SM*BSM interference
        tot = 0d0
        tot= xl(5)*sig(5)*(
     &          sft(5)+
     &          col(5) )
        tot=tot+xl(5)*vrt(5)
        tot=tot+2.d0*(xl1qq(5)+xl2qq(5))*sig(5)
        endif

      tot=tot*as
      endif

      return
      end


