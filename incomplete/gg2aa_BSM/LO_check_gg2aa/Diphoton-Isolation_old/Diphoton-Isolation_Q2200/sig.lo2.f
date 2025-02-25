c----------------------------------------C
c     2 -> 2 lo partonic cross section  -C
c----------------------------------------C

      subroutine sig_lo2(f1,f2,p1,p2,p3,p4,tot)
      implicit double precision (a-h,o-z)
      dimension f1(-6:6),f2(-6:6)
      dimension xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension sig(10)
      common/isub/io,is

      s=2d0*dot(p1,p2)
      t=-2d0*dot(p1,p3)
      u=-(s+t)

      call setlum(f1,f2,xl)
      call setsigma(s,t,u,sig)

      if (io.eq.1) then                 ! SM

        if (is.eq.1) then                 ! qqb
        tot=xl(1)*sig(1)
        elseif (is.eq.3) then             ! SM gg
        tot=xl(2)*sig(2)
        endif

      elseif (io.eq.2) then             ! BSM direct

        if (is.eq.1) then                 ! qqb
        tot=xl(3)*sig(3)
        elseif (is.eq.3) then             ! gg
        tot=xl(4)*sig(4)
        endif

      elseif (io.eq.3) then             ! SM*BSM interference

        if (is.eq.1) then                ! qqb
        tot=xl(5)*sig(5)
        elseif (is.eq.3) then             ! gg
        tot=xl(6)*sig(6)
        endif

      endif

      return
      end


