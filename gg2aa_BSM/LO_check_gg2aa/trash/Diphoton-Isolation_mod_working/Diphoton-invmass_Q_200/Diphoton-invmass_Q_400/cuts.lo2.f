
      subroutine cuts0(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,ipass)
      implicit double precision (a-h,o-z)

      common/machine/mid

      ipass=0
      
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
      
      

