      subroutine cuts3(p1,p2,p3,p4,p5,rsp,ipass,inf_PS)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)

      common/machine/mid
      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      common/amf/am3,am4,am5
      common/t_cuts/e_cut,t_cut


c-----------------------------------------------------------------------c
c       Technical cut off to avoid extreme singular phase space points
c-----------------------------------------------------------------------c

         e1 = p1(0)
         e2 = p2(0)
         e5 = p5(0)

         e5f = e5/rsp

         p1m = p1(1)**2.0d0 + p1(2)**2.0d0 + p1(3)**2.0d0
         p2m = p2(1)**2.0d0 + p2(2)**2.0d0 + p2(3)**2.0d0
         p5m = p5(1)**2.0d0 + p5(2)**2.0d0 + p5(3)**2.0d0

         ctheta15 = (e1*e5 - dot(p1,p5))/p1m/p5m
         ctheta25 = (e2*e5 - dot(p2,p5))/p1m/p5m

         cut1 = e_cut
         cut2 = 1.0d0 - t_cut

          inf_PS = 0
          if (e5f .le. cut1)       inf_PS =  1 
          if (ctheta15 .ge. cut2)  inf_PS =  1 
          if (ctheta25 .ge. cut2)  inf_PS =  1 
c-----------------------------------------------------------------------c



c-----------------------------------------------------------------------c
      ipass=0

c     yy1: rapidity of photons 3 
      y1 = rapid(p3)

c     yy2: rapiditiy of photon 4
      y2 = rapid(p4)

c     ppt1 ppt2
      pt1 = eperp(p3)
      pt2 = eperp(p4)

c     cone
      r34 = rjet(p3,p4)
      r35 = rjet(p3,p5)
      r45 = rjet(p4,p5)

c     transverse energy
      Et5 = eperp(p5)  ! for a massless particle Et=pt

c-----------------------------------------------------------------------c
      IF(mid.eq.0)THEN

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

c-----------------------------------------------------------------------c
      ELSEIF(mid.eq.1)THEN
c         r0=0.4d0
c         rgg=0.4d0
c         ET_iso=15d0

c         if( (abs(y1).ge.2.50d0 .or. abs(y2).ge.2.50d0) ) return
c         if( (pt1.lt.40d0 .or. pt2.lt.25d0)
c     &        .and. (pt1.lt.25d0 .or. pt2.lt.40d0) )return

         pt_hard = max(pt1,pt2)
         pt_soft = min(pt1,pt2)


         if (dabs(y1) .ge. 2.50d0) return
         if (dabs(y2) .ge. 2.50d0) return

         if (pt_hard .le. 40.0d0) return
         if (pt_soft .le. 25.0d0) return

c         if(r34.lt.rgg)return
c         if(r35.gt.r0 .and. r45.gt.r0)then
c            ipass=1
c            return
c         endif
c
c         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
c         IF(r45.le.r0 .and. Et5.gt.ET_iso) return
c
c
c        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
c        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return
c-----------------------------------------------------------------------c

        ipass=1
         return
      ENDIF
      end
