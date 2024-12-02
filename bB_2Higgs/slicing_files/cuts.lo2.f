      subroutine cuts0(p1,p2,p3,p4,ipass)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)

c      common/machine/mid

      mid=1
      ipass=0
      
c     yy1: rapidity of photons 3
      yrpda  = (p3(0)+p3(3))/(p3(0)-p3(3))
      y1    = 0.5*dlog(yrpda)

c     yy2: rapiditiy of photon 4
      yrpdb  = (p4(0)+p4(3))/(p4(0)-p4(3))
      y2    = 0.5*dlog(yrpdb)

c     ppt1 ppt2
      pt1=dsqrt(p3(1)**2+p3(2)**2)
      pt2=dsqrt(p4(1)**2+p4(2)**2)

      IF(mid.eq.0)THEN          ! Tevatron
        if( (dabs(y1).ge.0.9d0 .or. dabs(y2).ge.0.9d0) ) return
        if( (pt1.lt.14.d0 .or. pt2.lt.14.0d0) )return
        
        ipass =1
        return          
      ELSEIF(mid.eq.1)THEN      !LHC     
         if( (dabs(y1).ge.2.50d0 .or. dabs(y2).ge.2.50d0) ) return
         if( (pt1.lt.40.d0 .or. pt2.lt.40.d0) ) return
         
         ipass =1
         return
      ELSE
         write(*,*)'UNKNOWN MACHINE'
         stop
      ENDIF
      
      end
