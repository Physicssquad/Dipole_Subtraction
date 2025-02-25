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
cc---------------------------------------------------------------------
cc---------------------------------------------------------------------
cc     Initial state dipole for the case of Drell- Yan gq channel
c
c       function dipole_gq_q(k,p)
c      implicit double precision (a-h,o-z)
c      parameter(PI=3.141592653589793238D0)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
c     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
c      common/usedalpha/AL,ge
c      call p2dtop1d_5(p,p1,p2,p3,p4,p5)
c
c      s12=2.d0*dot(p1,p2) 
c      s13=2.d0*dot(p1,p3)
c      s14=2.d0*dot(p1,p4)
c      s15=2.d0*dot(p1,p5)
c      s23=2.d0*dot(p2,p3)
c      s24=2.d0*dot(p2,p4)
c      s25=2.d0*dot(p2,p5)
c      s34=2.d0*dot(p3,p4)
c      s35=2.d0*dot(p3,p5)
c      s45=2.d0*dot(p4,p5)
c
c      Tr = 0.5d0
c
c      if ( k .eq. 1 ) then ! Dipole leg 1 
c        call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)
c        Born= born_uu2ee(1,p6,p7,p8,p9) 
c
c        dipole_gq_q=
c     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
c     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*s15*(s12 - s15 - s25))
c
c      else if ( k .eq. 2 ) then ! Diople leg 2
c        call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
c        Born= born_uu2ee(2,p6,p7,p8,p9) 
c
c        dipole_gq_q=
c     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
c     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*(s12 - s15 - s25)*s25)
c
c      endif  
c      return
c      end 
cc---------------------------------------------------------------------

c     Initial state dipole for the case of gg --> Higgs [2 ~~> 1+1jet]

       function dipole_type_1_gg_g(k,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      parameter(PI=3.141592653589793238D0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
     .         ,p5(0:3),p6(0:3),p7(0:3)
      common/usedalpha/AL,ge

      s12=2.d0*dot(p1,p2) 
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s34=2.d0*dot(p3,p4)

      Tr = 0.5d0
      CA = 3.0d0

      if ( k .eq. 1 ) then ! Dipole leg 1 
        call reducemomenta_type_1(1,p1,p2,p3,p4,p5,p6,p7) 
        Born= Born_gg2h(1,p5,p6,p7) 

        dipole_type_1_gg_g=
     -          (-16*AL*Born*CA*Pi*(s12**4 - 2*s12**3*(s14 + s24) + 
     -      s12**2*(s14 + s24)**2 - 2*s12*(s14 + s24)**3 + 
     -      (s14 + s24)**4))/(s12*s14*(s12 - s14 - s24)**2*(s14 + s24))

      else if ( k .eq. 2 ) then ! Diople leg 2
        call reducemomenta_type_1(2,p1,p2,p3,p4,p5,p6,p7)
        Born= born_gg2h(2,p5,p6,p7) 

        dipole_type_1_gg_g=
     -          (-16*AL*Born*CA*Pi*(s12**4 - 2*s12**3*(s14 + s24) + 
     -      s12**2*(s14 + s24)**2 - 2*s12*(s14 + s24)**3 + 
     -      (s14 + s24)**4))/(s12*(s12 - s14 - s24)**2*s24*(s14 + s24))

      endif  
      return
      end 
c---------------------------------------------------------------------

c---------------------------------------------------------------------c
c       Initial State splitting with initial spectator                c
c            REDUCED MOMENTA                                          c
c---------------------------------------------------------------------c
       subroutine reducemomenta_type_1(k,p1,p2,p3,p4,p1til,p2til,p3til)
       implicit double precision (a-h,o-z)

       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
     .   p1til(0:3),p2til(0:3),p3til(0:3),ak(0:3),aktil(0:3)
     .   ,akadd(0:3),diff(0:3)

      s12=2.d0*dot(p1,p2)
      s13=2.d0*dot(p1,p3)
      s14=2.d0*dot(p1,p4)
      s23=2.d0*dot(p2,p3)
      s24=2.d0*dot(p2,p4)
      s34=2.d0*dot(p3,p4)
 
        x412=(s12-s14-s24)/s12
        if (k .eq .1) then   ! this choice is for leg1

          do i=0,3

           p1til(i) = x412*p1(i)
           p2til(i) = p2(i)

           ak(i)    = p1(i)+p2(i)-p4(i)
           aktil(i) = p1til(i)+p2(i)
           akadd(i) = ak(i)+aktil(i)
          enddo

        do i=0,3
        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
        enddo


        else if (k .eq. 2) then  ! this choice is for leg2

        do i=0,3

        p1til(i) = p1(i)
        p2til(i) = x412*p2(i)
           ak(i) = p1(i)+p2(i)-p4(i)
        aktil(i) = p2til(i)+p1(i)
        akadd(i) = ak(i)+aktil(i)

        enddo
        do i=0,3
        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
        enddo

        endif

        end
c---------------------------------------------------------------------

cc---------------------------------------------------------------------
cc     Initial state dipole for the case of Drell- Yan gg type of channel [2>2+1jet]
cc
c       function dipole_gg_g(k,p)
c      implicit double precision (a-h,o-z)
c      parameter(PI=3.141592653589793238D0)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
c     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
c      common/usedalpha/AL,ge
c
c      call p2dtop1d_5(p,p1,p2,p3,p4,p5)
c
c      s12=2.d0*dot(p1,p2) 
c      s13=2.d0*dot(p1,p3)
c      s14=2.d0*dot(p1,p4)
c      s15=2.d0*dot(p1,p5)
c      s23=2.d0*dot(p2,p3)
c      s24=2.d0*dot(p2,p4)
c      s25=2.d0*dot(p2,p5)
c      s34=2.d0*dot(p3,p4)
c      s35=2.d0*dot(p3,p5)
c      s45=2.d0*dot(p4,p5)
c
c      Tr = 0.5d0
c      CA = 3.0d0
c
c      if ( k .eq. 1 ) then ! Dipole leg 1 
c        call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)  ! reduce momenta2 is for initial splitting 
c        Born= Born_gg2h(1,p6,p7,p8,p9) 
c
c        dipole_gg_g=
c     -          (-16*AL*Born*CA*Pi*(s12**4 - 2*s12**3*(s15 + s25) + 
c     -      s12**2*(s15 + s25)**2 - 2*s12*(s15 + s25)**3 + 
c     -      (s15 + s25)**4))/(s12*s15*(s12 - s15 - s25)**2*(s15 + s25))
c
c
c      else if ( k .eq. 2 ) then ! Diople leg 2
c        call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
c        Born= born_gg2h(2,p6,p7,p8,p9) 
c
c        dipole_gg_g=
c     -          (-16*AL*Born*CA*Pi*(s12**4 - 2*s12**3*(s15 + s25) + 
c     -      s12**2*(s15 + s25)**2 - 2*s12*(s15 + s25)**3 + 
c     -      (s15 + s25)**4))/(s12*(s12 - s15 - s25)**2*s25*(s15 + s25))
c
c      endif  
cc      print*,"Born",Born
c      return
c      end 
cc---------------------------------------------------------------------
c
cc---------------------------------------------------------------------c
cc       Initial State splitting with initial spectator                c
cc            REDUCED MOMENTA                                          c
cc---------------------------------------------------------------------c
c       subroutine reducemomenta2(k,p1,p2,p3,p4,p5,p15,p2til,p3til,p4til)
c       implicit double precision (a-h,o-z)
c
c       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
c     .   p15(0:3),p2til(0:3),p3til(0:3),p4til(0:3),ak(0:3),aktil(0:3)
c     .   ,akadd(0:3),p25(0:3),p1til(0:3),diff(0:3)
c      s12=2.d0*dot(p1,p2)
c      s13=2.d0*dot(p1,p3)
c      s14=2.d0*dot(p1,p4)
c      s15=2.d0*dot(p1,p5)
c      s23=2.d0*dot(p2,p3)
c      s24=2.d0*dot(p2,p4)
c      s25=2.d0*dot(p2,p5)
c      s34=2.d0*dot(p3,p4)
c      s35=2.d0*dot(p3,p5)
c      s45=2.d0*dot(p4,p5)
c
c        x512=(s12-s15-s25)/s12
c        if (k .eq .1) then   ! this choice is for leg1
c
c        do i=0,3
c
c        p15(i)=x512*p1(i)
c        p2til(i)=p2(i)
c        ak(i)=p1(i)+p2(i)-p5(i)
c        aktil(i)= p15(i)+p2(i)
c        akadd(i)=ak(i)+aktil(i)
c        enddo
c        do i=0,3
c        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
c     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
c        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
c     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
c        enddo
c        else if (k .eq. 2) then  ! this choice is for leg2
c
c        do i=0,3
c
c        p25(i)=x512*p2(i)
c        p1til(i)=p1(i)
c        ak(i)=p1(i)+p2(i)-p5(i)
c        aktil(i)= p25(i)+p1(i)
c        akadd(i)=ak(i)+aktil(i)
c
c        p15(i) = p1til(i) ! in output p6,p7,p8,p9 numbered accordingly
c        p2til(i)= p25(i)
c
c        enddo
c        do i=0,3
c        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
c     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
c        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
c     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
c        enddo
c        endif
c
c        end
cc---------------------------------------------------------------------
