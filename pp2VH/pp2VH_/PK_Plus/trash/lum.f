      subroutine setlum(fa,fb,f)
      implicit double precision (a-h,o-z)
      dimension fa(-6:6),fb(-6:6)
      dimension f(15)
      common/machine/mid
       mid = 1
C----------------------------------------
C    SM             qqb=1     gg=2      C
C    BSM            qqb=3     gg=4      C
C    intf: SM-BSM   qqb=5     gg=6      C

C    SM             qg=7      gq=8      C
C    BSM            qg=9      gq=10     C
C    intf: SM-BSM   qg=11     gq=12     C
C----------------------------------------

C     DIGAMMA
C qqb -sm
C---------

      If(mid.eq.1)Then          ! LHC
         f(1)= !(fa(-1)*fb(+1)+fa(+1)*fb(-1))/81d0
     .        + (fa(-2)*fb(+2)+fa(+2)*fb(-2))*4d0/9d0
c     .        + (fa(-3)*fb(+3)+fa(+3)*fb(-3))/81d0
c     .        + (fa(-4)*fb(+4)+fa(+4)*fb(-4))*16d0/81d0
c     .        + (fa(-5)*fb(+5)+fa(+5)*fb(-5))/81d0
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(1)= (fa(-1)*fb(-1)+fa(+1)*fb(1))/81d0
     .        + (fa(-2)*fb(-2)+fa(+2)*fb(2))*16d0/81d0
     .        + (fa(-3)*fb(-3)+fa(+3)*fb(3))/81d0
     .        + (fa(-4)*fb(-4)+fa(+4)*fb(4))*16d0/81d0
     .        + (fa(-5)*fb(-5)+fa(+5)*fb(5))/81d0
      ENDIF

C gg -sm
C---------

      If(mid.eq.1)Then          ! LHC
         f(2) =fa(0)*fb(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(2) =fa(0)*fb(0)
      ENDIF


C qqb- bsm
C---------
      IF(mid.eq.1)THEN          ! LHC         
         f(3)=  (fa(-1)*fb(+1)+fa(+1)*fb(-1))
     .        +  (fa(-2)*fb(+2)+fa(+2)*fb(-2))
     .        +  (fa(-3)*fb(+3)+fa(+3)*fb(-3))
     .        +  (fa(-4)*fb(+4)+fa(+4)*fb(-4))
     .        +  (fa(-5)*fb(+5)+fa(+5)*fb(-5))
      ELSEIF(mid.eq.0)then      ! Tevatron
         f(3)=  (fa(-1)*fb(-1)+fa(+1)*fb(+1))
     .        +  (fa(-2)*fb(-2)+fa(+2)*fb(+2))
     .        +  (fa(-3)*fb(-3)+fa(+3)*fb(+3))
     .        +  (fa(-4)*fb(-4)+fa(+4)*fb(+4))         
     .        +  (fa(-5)*fb(-5)+fa(+5)*fb(+5))
      endif

C gg - bsm
C---------

         f(4) =fa(0)*fb(0)

C qqb -int
C---------

      If(mid.eq.1)Then          ! LHC
         f(5)= (fa(-1)*fb(+1)+fa(+1)*fb(-1))/9d0
     .        + (fa(-2)*fb(+2)+fa(+2)*fb(-2))*4d0/9d0
     .        + (fa(-3)*fb(+3)+fa(+3)*fb(-3))/9d0
     .        + (fa(-4)*fb(+4)+fa(+4)*fb(-4))*4.d0/9d0
     .        + (fa(-5)*fb(+5)+fa(+5)*fb(-5))/9d0
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(5)= (fa(-1)*fb(-1)+fa(+1)*fb(1))/9d0
     .        + (fa(-2)*fb(-2)+fa(+2)*fb(2))*4d0/9d0
     .        + (fa(-3)*fb(-3)+fa(+3)*fb(3))/9d0
     .        + (fa(-4)*fb(-4)+fa(+4)*fb(4))*4.d0/9d0
     .        + (fa(-5)*fb(-5)+fa(+5)*fb(5))/9d0
      ENDIF

C gg-int
C-------
      If(mid.eq.1)Then          ! LHC
         f(6) =fa(0)*fb(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(6) =fa(0)*fb(0)
      ENDIF

C qg sm
C--------

      If(mid.eq.1)Then          ! LHC
      f(7)=(fa(1)+16.0d0*fa(2)+fa(3)+16.0d0*fa(4)+fa(5))/81.0d0*fb(0)
     .+(fa(-1)+16.0d0*fa(-2)+fa(-3)+16.0d0*fa(-4)+fa(-5))/81.0d0*fb(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
      f(7)=(fa(1)+16.0d0*fa(2)+fa(3)+16.0d0*fa(4)+fa(5))/81.0d0*fb(0)
     .+(fa(-1)+16.0d0*fa(-2)+fa(-3)+16.0d0*fa(-4)+fa(-5))/81.0d0*fb(0)
      ENDIF


C gq sm
C--------

      If(mid.eq.1)Then          ! LHC
      f(8)=(fb(1)+16.0d0*fb(2)+fb(3)+16.0d0*fb(4)+fb(5))/81.0d0*fa(0)
     .+(fb(-1)+16.0d0*fb(-2)+fb(-3)+16.0d0*fb(-4)+fb(-5))/81.0d0*fa(0)
      ELSEIF(mid.eq.0)Then      !TEVATRON
      f(8)=(fb(1)+16.0d0*fb(2)+fb(3)+16.0d0*fb(4)+fb(5))/81.0d0*fa(0)
     .+(fb(-1)+16.0d0*fb(-2)+fb(-3)+16.0d0*fb(-4)+fb(-5))/81.0d0*fa(0)
      ENDIF


C qg bsm
C-------
 
      If(mid.eq.1)Then          ! LHC
         f(9)= ( fa(1)+fa(2)+fa(3)+fa(4)+fa(5) )*fb(0)
     .        +( fa(-1)+fa(-2)+fa(-3)+fa(-4)+fa(-5) )*fb(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(9)= ( fa(1)+fa(2)+fa(3)+fa(4)+fa(5) )*fb(0)
     .        +( fa(-1)+fa(-2)+fa(-3)+fa(-4)+fa(-5) )*fb(0)
      ENDIF
 
C gq bsm
C-------

      If(mid.eq.1)Then          ! LHC
         f(10)=  ( fb(1)+fb(2)+fb(3)+fb(4)+fb(5) )*fa(0)
     .        +( fb(-1)+fb(-2)+fb(-3)+fb(-4)+fb(-5) )*fa(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
         f(10)=  ( fb(1)+fb(2)+fb(3)+fb(4)+fb(5) )*fa(0)
     .        +( fb(-1)+fb(-2)+fb(-3)+fb(-4)+fb(-5) )*fa(0)
      ENDIF
      
C qg intf
C--------

      If(mid.eq.1)Then          ! LHC
      f(11)=(fa(1)+4.0d0*fa(2)+fa(3)+4.0d0*fa(4)+fa(5))/9.0d0*fb(0)
     .+( fa(-1)+4.0d0*fa(-2)+fa(-3)+4.0d0*fa(-4)+fa(-5))/9.0d0*fb(0)
      ELSEIF(mid.eq.0)Then      ! TEVATRON
      f(11)=(fa(1)+4.0d0*fa(2)+fa(3)+4.0d0*fa(4)+fa(5))/9.0d0*fb(0)
     .+( fa(-1)+4.0d0*fa(-2)+fa(-3)+4.0d0*fa(-4)+fa(-5))/9.0d0*fb(0)
      ENDIF
      
      
C gq intf
C--------
      
      If(mid.eq.1)Then          ! LHC
      f(12)=(fb(1)+4.0d0*fb(2)+fb(3)+4.0d0*fb(4)+fb(5))/9.0d0*fa(0)
     .+(fb(-1)+4.0d0*fb(-2)+fb(-3)+4.0d0*fb(-4)+fb(-5))/9.0d0*fa(0)
      ELSEIF(mid.eq.0)Then      !TEVATRON
      f(12)=(fb(1)+4.0d0*fb(2)+fb(3)+4.0d0*fb(4)+fb(5))/9.0d0*fa(0)
     .+(fb(-1)+4.0d0*fb(-2)+fb(-3)+4.0d0*fb(-4)+fb(-5))/9.0d0*fa(0)
      ENDIF

      return
      end
      

