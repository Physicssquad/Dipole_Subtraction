      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),zz(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2),all_PS(0:3,1:4,0:2)
      dimension SumP(1:2),SumK(1:2)
      dimension pa1(0:3),pa2(0:3),pa3(0:3),pa4(0:3)
      dimension pb1(0:3),pb2(0:3),pb3(0:3),pb4(0:3)
      dimension xinvmass(0:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/leg_choice/leg
      common/usedalpha/AL,ge
      common/distribution/xq
      external PqqP,Pqqreg,Born_uU2eE
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(c)
      xmin   = 0.0d0
      xmax   = 1.0d0 - 1d-3
      xjac   = (xmax-xmin)
      x      = xmin+ xjac*yy(4)
      yy(4)  = x
      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac = xjac*2.0d0
      xnorm=hbarc2

      eps = 1.0d0
      xlow = xq - eps
      xhigh = xq + eps


c     all_PS -> (4momenta,particleID,kinematic type)
c        call kinvar2_PK(yy,xinvmass,all_PS)

        
        do k = 1,2

        if (k .eq. 1) then
        call kinvar2_PK(x*xa,xb,xc,xxinvmass,p1,p2,p3,p4)
        elseif (k .eq. 1) then
        call kinvar2_PK(xa,x*xb,xc,xxinvmass,p1,p2,p3,p4)
        endif

        call kinvar2_PK(xa,xb,xc,xxinvmass,p1,p2,p3,p4)


c        print*,"4momenta,particleID,kinematic type"
c        do k=0,2
c          do j=1,4
c            do i=0,3
c              print*,i,j,k,all_PS(i,j,k)
c            enddo
c          enddo
c          print*," "
c        enddo
c        print*,"invmass:",xinvmass
c

c        momenta_giver(kin_type,all_PS,p1,p2,p3,p4)
c        call momenta_giver(0,all_PS,p1,p2,p3,p4)
c        call momenta_giver(1,all_PS,pa1,pa2,pa3,pa4)
c        call momenta_giver(2,all_PS,pb1,pb2,pb3,pb4)

c        call kinvar2(yy,xinv,p1,p2,p3,p4)
c            do i=0,3
c             xp1(i) = x*p1(i)
c             xp2(i) = x*p2(i)
c            enddo
c            print*,"xp1:",xp1
c            print*,"xp2:",xp2
c        call momenta_giver(1,all_PS,p1,p2,p3,p4)
c        print*,"p1:",p1
c        call momenta_giver(2,all_PS,p1,p2,p3,p4)
c        print*,"p1:",p2
c        stop


        call p1d_to_p2d_4(p1,p2,p3,p4,p)

c        scalex2 = 2.0d0*x*dot(p1,p2)
c        scalex = dsqrt(scalex2)
        scalex = xinvmass(2) 

c        scale2 = 2.0d0*dot(p1,p2)
c        scale = dsqrt(scale2)
        scale = xinvmass(0) 

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp
        flux_x = 4.0d0*x*pin*rsp
        coef = Born_uU2eE(0,p1,p2,p3,p4)


        flo2_PK  = 0.0d0
        PKplus_x = 0.0d0
        PKplus_1 = 0.0d0
        PKRegDel = 0.0d0
        PKplus   = 0.0d0

c         goto 101
c         Block for + distribution x = x
        if (scalex .ge. xlow .and. scalex .le. xhigh) then   

             xmuf = scalex
            xmuf2 = xmuf*xmuf 
             xmur = scalex
               AL = alphasPDF(xmur)
              ALP = AL/2d0/Pi
              s12 = 2d0*dot(p1,p2)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)


            call getPK(1,x,xmuf,p,xp1,xp2,SumP,SumK)
                
c           coefPterm   = (-1.0d0)*dlog(xmuf2/s12/x)
c           SumP(1) = coefPterm*PqqP(x) 

c           SumK(1) = AKbarP_qq(x)+AKtilP_qq(x)

c            sig1 = xl(1)*( SumP(1) + SumK(1) ) 
            sig1 = xl(1)*( SumP(1) + SumK(1) ) 
  
            sig = Alp*sig1*coef

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt_x = sig/flux_x*ps2*xjac*vwgt
            PKplus_x = xnorm*wgt_x/vwgt/2d0/eps

         endif

c 101    continue

c        goto 102

        if (scale .ge. xlow .and. scale .le. xhigh)  then   

             xmuf = scale
            xmuf2 = xmuf*xmuf 
             xmur = scale
               AL = alphasPDF(xmur)
              ALP = AL/2d0/Pi
              s12 = 2d0*dot(p1,p2)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)



            call getPK(1,x,xmuf,p,xp1,xp2,SumP,SumK)

c           coefPterm   = (-1.0d0)*dlog(xmuf2/s12/1d0)
c           SumP(1) = coefPterm*PqqP(x) 
c
c           SumK(1) = AKbarP_qq(x)+AKtilP_qq(x)

c            sig1 = xl(1)*( SumP(1) + SumK(1) ) 
            sig1 = xl(1)* SumP(1) 
  
            sig = Alp*sig1*coef

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt = sig/flux_1*ps2*xjac*vwgt
            PKplus_1 = xnorm*wgt/vwgt/2d0/eps

         endif
         PKplus = ( PKplus_x - PKplus_1)

c  102     continue
c         goto 192

c       This is for regular and delta terms

        if (scalex .ge. xlow .and. scalex .le. xhigh) then

             xmuf = scalex
            xmuf2 = xmuf*xmuf 
             xmur = scalex
               AL = alphasPDF(xmur)
              ALP = AL/2d0/Pi
              s12 = 2d0*dot(p1,p2)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)


            call getPK(0,x,xmuf,p,xp1,xp2,SumP,SumK)

c           Areg = (AKbarreg_qq(x)+AKtilreg_qq(x))
c           ADel = (AKbarD_qq(x)+AKtilD_qq(x))

c           sig = Alp*(Areg + Adel)

            sig1 = xl(1)* SumP(1) 
            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt_x = sig/flux_x*ps2*xjac*vwgt
            PKReg = xnorm*wgt_x/vwgt/2d0/eps

         endif
c         
c        if (scale .ge. xlow .and. scale .le. xhigh)  then   
c
c           ADel = (AKbarD_qq(x)+AKtilD_qq(x))
c
c           sig = Alp*ADel
c
c           azmth = 2.0d0*pi
c            pf = 0.5d0*rsp
c            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
c
c            wgt = sig/flux_1*ps2*xjac*vwgt
c            PKDel = xnorm*wgt/vwgt/2d0/eps
c
c         endif
192      continue
         flo2_PK = PKplus + PKReg !+ PKDel

      return
      end

c---------------------------------------------------------------------
      subroutine momenta_giver(kin_type,all_PS,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
c     all_PS(4momenta,particleID,kinematic type)
       dimension all_PS(0:3,1:4,0:2)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
         do j=0,3
           p1(j) = all_PS(j,1,kin_type) 
           p2(j) = all_PS(j,2,kin_type) 
           p3(j) = all_PS(j,3,kin_type) 
           p4(j) = all_PS(j,4,kin_type) 
         enddo
      end
c---------------------------------------------------------------------
 

c---------------------------------------------------------------------
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       common/usedalpha/AL,ge
c       ge=1d0/128d0
       e= DSQRT(ge*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p1,p4) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      qu2 = 1d0!4d0/9d0

      Born_uU2eE= CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c---------------------------------------------------------------------
