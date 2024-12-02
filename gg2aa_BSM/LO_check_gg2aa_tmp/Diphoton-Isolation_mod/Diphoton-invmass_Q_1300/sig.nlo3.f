c. 
c. 2 -> 3 nlo cross section
c.
      subroutine sig_nlo3(f1,f2,p1,p2,p3,p4,p5,tot)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      dimension f1(-6:6),f2(-6:6)
      dimension xl(15)
      parameter (cf=4d0/3d0)

      common/machine/mid
      common/slice/deltas,deltac
      common/isub/io,is
      common/amf/am3,am4,am5

      s12=2d0*dot(p1,p2)

      t13=-2d0*dot(p1,p3)+am3**2.0d0
      t14=-2d0*dot(p1,p4)+am4**2.0d0
      t15=-2d0*dot(p1,p5)+am5**2.0d0

      t23=-2d0*dot(p2,p3)+am3**2.0d0
      t24=-2d0*dot(p2,p4)+am4**2.0d0
      t25=-2d0*dot(p2,p5)+am5**2.0d0

      s34=2d0*dot(p3,p4)+am3**2.0d0+am4**2.0d0
      s35=2d0*dot(p3,p5)+am3**2.0d0+am5**2.0d0
      s45=2d0*dot(p4,p5)+am4**2.0d0+am5**2.0d0

      rs12=dsqrt(s12)

      soft=0.5d0*deltas*rs12
      coll=deltac*s12

      i15=0
      i25=0

      is5=0

c     collinear
      if(dabs(t15).le.coll) i15=1
      if(dabs(t25).le.coll) i25=1

c     soft
      e5=0.5d0*(s12-s34)/rs12
      if(e5.le.soft) is5=1

      tot=0d0

      call setlum(f1,f2,xl)

      if (io.eq.1) then                 ! SM

       
        if (is.eq.1) then               ! qqb NLO3

c           write(*,*)'Input parameters',io,is,tot
           itest=i15+i25+is5
           if(itest.eq.0)then         
           call SMQQB(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMQQB)
           qqb_sm=SSMQQB
           tot=tot+xl(1)*qqb_sm
           endif

        elseif (is.eq.2) then           ! qg NLO3

           itest=i25
           if(itest.eq.0)then
           call SMQG(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMQG)
           qg_sm=SSMQG
           tot=tot+xl(7)*qg_sm
           endif

           itest=i15
           if(itest.eq.0)then
           call SMGQ(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMGQ)
           gq_sm=SSMGQ
           tot=tot+xl(8)*gq_sm
           endif

        endif


      elseif (io.eq.2) then             ! BSM direct

        if (is.eq.1) then               ! qqb
        itest=i15+i25+is5

        if(itest.eq.0)then         
        call GRQQB(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SGRQQB)
        qqb_bsm=SGRQQB
        tot=tot+xl(3)*qqb_bsm
        endif

        elseif (is.eq.2) then             ! qg
        itest=i15+i25

        if(itest.eq.0)then
        call GRQG(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SGRQG)
        qg_bsm=SGRQG
        tot=tot+xl(9)*qg_bsm
        endif

        itest=i15+i25
        if(itest.eq.0)then
        call GRGQ(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SGRGQ)
        gq_bsm=SGRGQ
        tot=tot+xl(10)*gq_bsm
        endif

        elseif (is.eq.3) then           ! gg
        itest=i15+i25+is5

        if(itest.eq.0)then         
        call GRGG(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SGRGG)
        gg_bsm=SGRGG
        tot=tot+xl(4)*gg_bsm
        endif

        endif

      elseif (io.eq.3) then             ! SM*BSM interference

        if (is.eq.1) then               ! qqb
        itest=i15+i25+is5

        if(itest.eq.0)then         
        call SMGRQQB(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMGRQQB)
        qqb_intf_bsm=SSMGRQQB
        tot=tot+xl(5)*qqb_intf_bsm
        endif

        elseif (is.eq.2) then           ! qg
        itest=i25+i15

        if(itest.eq.0)then
        call SMGRQG(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMGRQG)
        qg_intf_bsm=SSMGRQG
        tot=tot+xl(11)*qg_intf_bsm
        endif

        itest=i15+i25
        if(itest.eq.0)then
        call SMGRGQ(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,SSMGRGQ)
        gq_intf_bsm=SSMGRGQ
        tot=tot+xl(12)*gq_intf_bsm
        endif

        endif

      endif

      return
      end


