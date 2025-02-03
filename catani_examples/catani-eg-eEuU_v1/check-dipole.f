C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%  This is main driver routine Edit accordingly  %%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------------
       subroutine printoutput(i,a,b)
       implicit double precision (a-h,o-z)
       if(i .eq. 1) then
       write(*,*)"            helas                     ours       
     .      (helas-our)/helas*1000        ratio"
       endif
       write(*,"(3e27.15,1x,3e27.15,i3)")a,b,dabs(a-b)/a*1000,a/b
c       write(*,*)a,b,b/a
       end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%     ALL DIPOLES ARE COLLECTED HERE   %%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------------
C                 %%%   Dipole Initial State    %%%
c                  !! g-g splitting with g spec !!
c                 "k" : arg "1" leg1 & arg "2" leg2
c                        initial state dipoles
c---------------------------------------------------------------------
      function dipole_gg_g(k,p)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
      Pi=3.1412d0
      AL=0.118d0
      call p2dtop1d_5(p,p1,p2,p3,p4,p5)

      s12=2.d0*dot(p1,p2) 
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s15=2.d0*dot(p1,p5)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s25=2.d0*dot(p2,p5)
      s34=2.d0*dot(p3,p4)
      s35=2.d0*dot(p3,p5)
      s45=2.d0*dot(p4,p5)

      if (k .eq. 1 ) then ! g 2 gg at leg1
       call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)
       Born = born_gg2uux(1,p6,p7,p8,p9)
       anum2 = (s12**2 -s12*(s15+s25)+(s15+s25)**2)**2
       aden2 = s12*s15*(s15+s25)*(s15+s25-s12)**2
       dipole_gg_g= -16d0*Pi*AL*anum2/aden2*Born

      else if (k .eq. 2) then      ! g 2 gg at leg2
       call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
       Born = born_gg2uux(2,p6,p7,p8,p9)

       anum2 = (s12**2 -s12*(s15+s25)+(s15+s25)**2)**2
       aden2 = s12*s25*(s15+s25)*(s15+s25-s12)**2
       dipole_gg_g= -16d0*Pi*AL*anum2/aden2*Born
      endif

      return
      end
c---------------------------------------------------------------------
c           FinalStateDipole(e E > u U g)
c---------------------------------------------------------------------
      function dipole_uu_g(k,p1,p2,p3,p4,p5)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
     .         ,p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:4)
      parameter(Pi=3.141592653589793238D0)
      AL=0.118d0
      call reducemomenta1(k,p1,p2,p3,p4,p5,p6,p7,p8,p9)
c      call p1dtop2d_4(p1,p2,p3,p4,p)
c      call smatrix(p,msq)
c      print*,"BOOOOOOOORN:   ",msq
      Born= Born_eE2uU_m(k,p6,p7,p8,p9)
c      Born= Born_eE2uU(k,p6,p7,p8,p9)
      s12=2.d0*dot(p1,p2)
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s15=2.d0*dot(p1,p5)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s25=2.d0*dot(p2,p5)
      s34=2.d0*dot(p3,p4)
      s35=2.d0*dot(p3,p5)
      s45=2.d0*dot(p4,p5)
      IF( k .eq. 1d0) then 
       anum1 = (2.0d0*s34**2 +2.0d0*s34*s45 +s45**2 +s45*s35)
       aden1 = (s34+s45)*s35*(s35+s45)
       dipole_uu_g = -8.0d0*Pi*AL*anum1/aden1*Born
      ELSE IF( k .eq. 2d0) THEN
       anum1 = (2.0d0*s34**2 +2.0d0*s34*s35 +s35**2 +s35*s45)
       aden1 = (s35+s45)*s45*(s34+s35)
       dipole_uu_g = -8.0d0*Pi*AL*anum1/aden1*Born
      ENDIF
      return
      end
c---------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%L=
C%%%%%%%%%%%%%%%      THIS IS THE END OF DIPOLES      %%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%     ALL BORN EXPRESSIONS ARE  HERE   %%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------------
       function Born_gg2uux(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       ge=0.007547169811320755d0
       e=DSQRT(ge*4.d0*PI)
       AS=0.118d0
       gs=DSQRT(AS*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0               !Leading Order K=0 
      IF(k .eq. 1)  CF =  1d0               !leg 1
      IF(k .eq. 2)  CF =  1d0               !Leg 2    
       s13 = - 2.0d0*dot(p1,p3) ! t
       s23 = - 2.0d0*dot(p2,p3) ! u
       s12 =   2.0d0*dot(p3,p4) ! s
       qu2 = 4.0d0/9.0d0
       s=s12
       t=s13
       u=s23
c       qd2 = 1.0d0/9.0d0
       qd2 = 0d0
       CF  = 1d0
       XNC = 1d0
       xnorm =1d0 

c      BornN = -gs**4*(t**2+u**2)*(s**2-9d0*t**2-9d0*u**2)
c      BornD = 48d0*s**2*t*u
c      Born_gg2uux  = BornN/BornD*CF*XNORM

      Born = 1d0*(s13**2+s23**2)/(s13*s23*6d0) -
     .        3d0*(s13**2+s23**2)/(8d0*s12**2)
      Born_gg2uux = Born*CF*XNC*gs**4

       return
       end
c--------------------------------------------------------------------o
c     [e E -> u U]  Born with mass-less electron pair
c--------------------------------------------------------------------o
       function Born_eE2uU(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       ge=0.007547169811320755d0
       AS=0.118d0
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(AS*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0 
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2    
      s13 = - 2.0d0*dot(p1,p3) ! t
      s23 = - 2.0d0*dot(p2,p3) ! u
      s12 =   2.0d0*dot(p3,p4) ! s
      qu2 =   4.0d0/9.0d0
      s=s12
      t=s13
      u=s23
c       qd2 = 1.0d0/9.0d0
      qd2 = 0.0d0
      CF  = CF * 3d0
      XNC = 1/4d0
      xnorm =1d0 
      Born = (u**2+t**2)/s**2
      Born_eE2uux = 8d0*e**4*Born*CF*XNC*qu2
       return
       end
c---------------------------------------------------------------------
c     [e E -> u U]  born with massive electron pair
c--------------------------------------------------------------------o
       function Born_eE2uU_m(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       double precision me,am1,am2,am3,am4,am5
      common/amass/am1,am2,am3,am4,am5
       me = am1
       ge=0.007547169811320755d0
       AS=0.118d0
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(AS*4.d0*PI)
c       me = 0.51099895000d-3
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0 
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2    
       s13 = - 2.0d0*dot(p1,p3) ! t
       s23 = - 2.0d0*dot(p2,p3) ! u
       s12 =   2.0d0*dot(p3,p4) ! s
       qu2 =   4.0d0/9.0d0
       s=s12
       t=s13
       u=s23
c        qd2 = 1.0d0/9.0d0
       qd2 = 0.0d0
       CF  = CF * 3d0
       XNC = 1/4d0

       Born=2d0*(-8d0*me**4+4d0*me**2*s+t**2+u**2)/6/s**2
       Born_eE2uU_m = e**4*Born*CF*XNC
      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
       function Born_uU2eE(p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       ge=0.007547169811320755d0
       e=DSQRT(ge*4.d0*PI)

       s12 = 2.0d0*dot(p1,p2)
       qu2 = 4.0d0/9.0d0
c       qd2 = 1.0d0/9.0d0
       qd2 = 0.0d0
       CF = 4.0d0/3.0d0
       XNC = 3.0d0

        Born=-8.0d0*e**4*(dot(p1,p3)*dot(p2,p4)
     &            + dot(p1,p4)*dot(p2,p3))/s12**2
        Born= Born*(1.0d0*qu2 + 3.0d0*qd2)
        Born = Born*CF*XNC

       return
       end
c---------------------------------------------------------------------

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%   ALL BORN EXPRESSIONS ENDS HERE   %%%%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c---------------------------------------------------------------------

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%% ALL real emission |M^2| expression %%%%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c---------------------------------------------------------------------
       function eE2uU_r(p1,p2,p3,p4,p5)
       implicit double precision(a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       double precision me2
       double precision me,am1,am2,am3,am4,am5
       common/amass/am1,am2,am3,am4,am5
       parameter(Pi=3.141592653589793238D0)
       me = am1
       s12=2.d0*dot(p1,p2)
       s13=2.d0*dot(p1,p3)
       s14=2.d0*dot(p1,p4)
       s15=2.d0*dot(p1,p5)
       s23=2.d0*dot(p2,p3)
       s24=2.d0*dot(p2,p4)
       s25=2.d0*dot(p2,p5)
       s34=2.d0*dot(p3,p4)
       s35=2.d0*dot(p3,p5)
       s45=2.d0*dot(p4,p5)

c       me2=(0.51099895000d-3)**2
       me2=me**2
       AL=0.118d0
       Prop1=1d0/(s12 + 2d0*me2)
       Prop2=1d0/(2d0*me2 + s12 - s13 - s23)
       Prop3=1d0/(s12 - s14 - s24 + 2d0*me2)
       ge=0.007547169811320755d0
       qe=DSQRT(ge*4.d0*PI)
       gs=DSQRT(AL*4.d0*PI)
       ge4gs2=qe**4*gs**2
       qu2 = 4d0/9d0
c       ge4gs2=qe**4**gs**2*qu2
c       write(*,*)'me^2 =',me2
c       write(*,*)'ge4gs2 =',ge4gs2
       CL= 4d0
       SP=1d0/4d0
           sigma=  Prop1**2*Prop2**2*(64*ge4gs2*me2*s13*s14 + 
     -     16*ge4gs2*s12*s13*s14 + 32*ge4gs2*me2*s14*s23 - 
     -     16*ge4gs2*s13*s14*s23 + 16*ge4gs2*s14*s23**2 + 
     -     32*ge4gs2*me2*s13*s24 + 16*ge4gs2*s13**2*s24 + 
     -     64*ge4gs2*me2*s23*s24 + 16*ge4gs2*s12*s23*s24 - 
     -     16*ge4gs2*s13*s23*s24 - 64*ge4gs2*me2**2*s34 - 
     -     32*ge4gs2*me2*s12*s34 - 32*ge4gs2*me2*s13*s34 - 
     -     16*ge4gs2*s12*s13*s34 - 32*ge4gs2*me2*s23*s34 - 
     -     16*ge4gs2*s12*s23*s34 + 32*ge4gs2*s13*s23*s34) + 
     -       Prop1**2*Prop3**2*(64*ge4gs2*me2*s13*s14 + 
     -     16*ge4gs2*s12*s13*s14 + 32*ge4gs2*me2*s14*s23 + 
     -     16*ge4gs2*s14**2*s23 + 32*ge4gs2*me2*s13*s24 - 
     -     16*ge4gs2*s13*s14*s24 + 64*ge4gs2*me2*s23*s24 + 
     -     16*ge4gs2*s12*s23*s24 - 16*ge4gs2*s14*s23*s24 + 
     -     16*ge4gs2*s13*s24**2 - 64*ge4gs2*me2**2*s34 - 
     -     32*ge4gs2*me2*s12*s34 - 32*ge4gs2*me2*s14*s34 - 
     -     16*ge4gs2*s12*s14*s34 - 32*ge4gs2*me2*s24*s34 - 
     -     16*ge4gs2*s12*s24*s34 + 32*ge4gs2*s14*s24*s34) + 
     -      Prop1**2*Prop2*Prop3*
     -     (-16*ge4gs2*s13*s14*s23 + 16*ge4gs2*s14**2*s23 + 
     -     16*ge4gs2*s14*s23**2 + 16*ge4gs2*s13**2*s24 - 
     -     16*ge4gs2*s13*s14*s24 - 16*ge4gs2*s13*s23*s24 - 
     -     16*ge4gs2*s14*s23*s24 + 16*ge4gs2*s13*s24**2 + 
     -     256*ge4gs2*me2**2*s34 + 192*ge4gs2*me2*s12*s34 + 
     -     32*ge4gs2*s12**2*s34 - 32*ge4gs2*me2*s13*s34 - 
     -     16*ge4gs2*s12*s13*s34 - 32*ge4gs2*me2*s14*s34 - 
     -     16*ge4gs2*s12*s14*s34 - 32*ge4gs2*me2*s23*s34 - 
     -     16*ge4gs2*s12*s23*s34 - 32*ge4gs2*s14*s23*s34 - 
     -     32*ge4gs2*me2*s24*s34 - 16*ge4gs2*s12*s24*s34 - 
     -     32*ge4gs2*s13*s24*s34 + 64*ge4gs2*me2*s34**2 + 
     -     32*ge4gs2*s12*s34**2)
           eE2uU_r=sigma*CL*SP*qu2
c           write(*,*)'Matrix elements =', eE2uU_r
      return
      end
c---------------------------------------------------------------------

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%   ALL real Emission |m^2|  ends    %%%%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%    All momenta modification are here  %%%%%%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c---------------------------------------------------------------------
c            Final state gluon emission both legs 
c---------------------------------------------------------------------
       subroutine reducemomenta1(k,p1,p2,p3,p4,p5,p6,p7,p8,p9)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
     .           ,p8(0:3),p9(0:3),p35til(0:3),p4til(0:3),p45til(0:3)
     .           ,p3til(0:3)

      s12=2.d0*dot(p1,p2) 
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s15=2.d0*dot(p1,p5)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s25=2.d0*dot(p2,p5)
      s34=2.d0*dot(p3,p4)
      s35=2.d0*dot(p3,p5)
      s45=2.d0*dot(p4,p5)
      IF (k .eq. 1) THEN
        Do i=0,3
         p35til(i) = p3(i) + p5(i) -s35*p4(i)/(s45+s34)
         p4til(i)  = (s35 + s45 + s34)*p4(i)/(s45+s34)

         p6(i) = p1(i)
         p7(i) = p2(i)
         p8(i) = p35til(i)
         p9(i) = p4til(i)
        End Do
      ELSE IF( k .eq. 2 ) Then
        Do i=0,3
         p3til(i)  = (s45 + s35 + s34)*p3(i)/(s35+s34)
         p45til(i) = p4(i) + p5(i) -s45*p3(i)/(s35+s34)

         p6(i) = p1(i)
         p7(i) = p2(i)
         p8(i) = p3til(i)
         p9(i) = p45til(i)
        End Do
      End If
      End
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c       Initial State splitting with initial spectator
c---------------------------------------------------------------------
       subroutine reducemomenta2(k,p1,p2,p3,p4,p5,p15,p2til,p3til,p4til)
       implicit double precision (a-h,o-z)

       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .   p15(0:3),p2til(0:3),p3til(0:3),p4til(0:3),ak(0:3),aktil(0:3)
     .   ,akadd(0:3),p25(0:3),p1til(0:3)
      s12=2.d0*dot(p1,p2) 
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s15=2.d0*dot(p1,p5)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s25=2.d0*dot(p2,p5)
      s34=2.d0*dot(p3,p4)
      s35=2.d0*dot(p3,p5)
      s45=2.d0*dot(p4,p5)

        x512=(s12-s15-s25)/s12
        if (k .eq .1) then   ! this choice is for leg1
        do i=0,3

        p15(i)=x512*p1(i)
        p2til(i)=p2(i)        
        ak(i)=p1(i)+p2(i)-p5(i)
        aktil(i)= p15(i)+p2(i) 
        akadd(i)=ak(i)+aktil(i)
        enddo
        do i=0,3
        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
        enddo
        else if (k .eq. 2) then  ! this choice is for leg2

        do i=0,3

        p25(i)=x512*p2(i)
        p1til(i)=p1(i)        
        ak(i)=p1(i)+p2(i)-p5(i)
        aktil(i)= p25(i)+p1(i) 
        akadd(i)=ak(i)+aktil(i)
        p15(i) = p1til(i) ! in output p6,p7,p8,p9 numbered accordingly
        p2til(i)= p25(i)
        enddo
        do i=0,3
        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
        enddo
        endif

        end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine p1dtop2d_4(p1,p2,p3,p4,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
        enddo
         end
c---------------------------------------------------------------------

c---------------------------------------------------------------------
         subroutine p2dtop1d_5(p,p1,p2,p3,p4,p5)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
        do i=0,3
         p1(i)=p(i,1)
         p2(i)=p(i,2)
         p3(i)=p(i,3)
         p4(i)=p(i,4)
         p5(i)=p(i,5)
        enddo
         end 
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine p1dtop2d_5(p1,p2,p3,p4,p5,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
         p(i,5)=p5(i)
        enddo
         end 
c---------------------------------------------------------------------
         subroutine p2dtop1d_4(p,p1,p2,p3,p4)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p1(i)=p(i,1)
         p2(i)=p(i,2)
         p3(i)=p(i,3)
         p4(i)=p(i,4)
        enddo
         end 
c---------------------------------------------------------------------

c---------------------------------------------------------------------
         subroutine momentareader(ptil)
                 implicit double precision (a-h,o-z)
                 dimension ptil(0:3,1:4)
                 write(*,*)"momenta tilde= ",ptil
         end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine printmomenta(p1,p2,p3,p4,p5)
                 implicit double precision (a-h,o-z)
                 dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
                 write(*,*)"p1= ",p1
                 write(*,*)"p2= ",p2
                 write(*,*)"p3= ",p3
                 write(*,*)"p4= ",p4
                 write(*,*)"p5= ",p5
         end
c---------------------------------------------------------------------

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%    All momenta modification END here  %%%%%%%%%%%%%%%%%%%%     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
