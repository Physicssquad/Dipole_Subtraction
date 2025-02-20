c---------------------------------------------------------------------
c     Initial state dipole for the case of Drell- Yan.
c---------------------------------------------------------------------

       function dipole_uU_g(k,p)
      implicit double precision (a-h,o-z)
      parameter(PI=3.141592653589793238D0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
      common/set/set1
      common/usedalpha/AL,ge
      call p2dtop1d_5(p,p1,p2,p3,p4,p5)
c      AL=0.118d0
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

      if ( k .eq. 1 ) then ! Dipole leg 1 
        call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)
        Born= born_uu2ee(1,p6,p7,p8,p9) 
       dipole_uU_g=(8*AL*Pi*(s15**2 - 2*s15*s12 + 2*s12**2 + 2*s15*s25 -
     .      2*s12*s25 + s25**2)*Born)/
     .  (s15*(s15 + s25)*(s15 - s12 + s25))

      else if ( k .eq. 2 ) then ! Diople leg 2
        call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
        Born= born_uu2ee(2,p6,p7,p8,p9) 

        dipole_uU_g=(8*AL*Pi*(s15**2 -2*s15*s12 + 2*s12**2 + 2*s15*s25 -
     .      2*s12*s25 + s25**2)*Born)/
     .  (s25*(s15 + s25)*(s15 - s12 + s25))
      endif  
      return
      end 
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     Initial state dipole for the case of Drell- Yan gq channel

       function dipole_gq_q(k,p)
      implicit double precision (a-h,o-z)
      parameter(PI=3.141592653589793238D0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
      common/usedalpha/AL,ge
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

      Tr = 0.5d0

      if ( k .eq. 1 ) then ! Dipole leg 1 
        call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)
        Born= born_uu2ee(1,p6,p7,p8,p9) 

        dipole_gq_q=
     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*s15*(s12 - s15 - s25))

      else if ( k .eq. 2 ) then ! Diople leg 2
        call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
        Born= born_uu2ee(2,p6,p7,p8,p9) 

        dipole_gq_q=
     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*(s12 - s15 - s25)*s25)

      endif  
      return
      end 
c---------------------------------------------------------------------


c---------------------------------------------------------------------c
c       Initial State splitting with initial spectator                c
c            REDUCED MOMENTA                                          c
c---------------------------------------------------------------------c
       subroutine reducemomenta2(k,p1,p2,p3,p4,p5,p15,p2til,p3til,p4til)
       implicit double precision (a-h,o-z)

       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .   p15(0:3),p2til(0:3),p3til(0:3),p4til(0:3),ak(0:3),aktil(0:3)
     .   ,akadd(0:3),p25(0:3),p1til(0:3),diff(0:3)
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
