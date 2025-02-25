
      subroutine cuts3(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
     &     r34,r35,r45,Et5,QT,ipass)      
      implicit double precision (a-h,o-z)

      common/machine/mid
      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      
      ipass=0

      IF(mid.eq.0)THEN         
c         r0=0.4d0
c         rgg=0.3d0
c         ET_iso=1.0d0
         
         if( (abs(y1).ge.0.9d0 .or. abs(y2).ge.0.9d0) ) return          
         if( (pt1.lt.14d0 .or. pt2.lt.13d0) 
     &        .and. (pt1.lt.13d0 .or. pt2.lt.14d0) )return
         
         if(r34.lt.rgg)return         
         if(r35.gt.r0 .and. r45.gt.r0)then
            ipass=1
            return
         endif                  

         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
         IF(r45.le.r0 .and. Et5.gt.ET_iso) return                 
         

        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return
               
         
         ipass=1               
         return

      ELSEIF(mid.eq.1)THEN
c         r0=0.4d0
c         rgg=0.4d0
c         ET_iso=15d0

         if( (abs(y1).ge.2.50d0 .or. abs(y2).ge.2.50d0) ) return
         if( (pt1.lt.40d0 .or. pt2.lt.25d0) 
     &        .and. (pt1.lt.25d0 .or. pt2.lt.40d0) )return
         
         if(r34.lt.rgg)return         
         if(r35.gt.r0 .and. r45.gt.r0)then
            ipass=1
            return
         endif            

         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
         IF(r45.le.r0 .and. Et5.gt.ET_iso) return                 
         

        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return
        
        ipass=1               
         return
      ENDIF
      end
      
