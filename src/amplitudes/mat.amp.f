c---------------------------------------------------------------------
c       This file contains      
c     1.  Born / reduced_Born gg2h in EFT(mass of top infinite approxomation) 
c     2.  1- real matrix amplitude. 
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     [g g -> h]  Born Effective field theory
c--------------------------------------------------------------------o
       function Born_gg2h(k,p1,p2,p3)
       use openloops
         implicit double precision (a-h,o-z)
         dimension p1(0:3),p2(0:3),p3(0:3),p_ex(0:3,3)
         parameter(PI=3.141592653589793238D0)
         common/usedalpha/AL,ge
         common/param2/xmur
         common/amass/am1,am2,amh,am4,am5
         common/energy/s
      common/caller/icall,id
      common/prc_id/id_LO,id_NLO_1r

         s12 = 2d0*dot(p1,p2)
         rp34  = dsqrt(s12)
         
         IF(k .eq. 0)  CF =  1d0               !Leading Order k=0
         IF(k .eq. 1)  CF = - 1d0               !leg 1 reduced born k=1 
         IF(k .eq. 2)  CF = - 1d0               !Leg 2 reduced born k=2
c CF is the colour factor in comes up in reduced born.
           NA = 8
            v = 246d0
           AS = AL/4d0/PI
           ch = -4d0*AS/3d0/v 
          ch2 = ch * ch
           avg_pol = 4d0

          Born_gg2H =  CF*s12**2*ch2/64d0 ! this factor 1.002008 is a missing constant factor.
      ! Call openloops born using these routines.        
c        call p1d_to_p2d_3(p1,p2,p3,p_ex) 
c        call evaluate_tree(id_LO, p_ex,answer)
c        Born_gg2H  = CF*answer

       return
       end
c---------------------------------------------------------------------

       subroutine amp_mat_r(p1,p2,p3,p4,sig)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
     . ,a(1:100)
       parameter(Pi=3.141592653589793238D0)
       double precision msq(5),msq1,msq2,lambda
       common/usedalpha/AL,ge
       common/scales/xinvmass
       common/amass/am1,am2,amh,am4,am5


               v = 246d0
               AS = AL/4d0/PI
              ch = -4d0*AS/3d0/v
              ch2 = ch*ch

              gs = DSQRT(AL*4.d0*PI)
          gs2CH2 = gs**2*ch2

           s12 = 2d0*dot(p1,p2)
           t11 = 2d0*dot(p1,p4)
           t21 = 2d0*dot(p2,p4)


         sig = (48d0*gs2CH2*(s12**4 - 2d0*s12**3*(t11 + t21) + 
     -        3d0*s12**2*(t11 + t21)**2 - 2d0*s12*(t11 + t21)**3 + 
     -        (t11**2 + t11*t21 + t21**2)**2))/(s12*t11*t21)/256d0

      return
      end
c---------------------------------------------------------------------
