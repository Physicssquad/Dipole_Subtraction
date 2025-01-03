c---------------------------------------------------------------------
c       This file contains      
c     1.  Born / reduced_Born gg2h in EFT(mass of top infinite approxomation) 
c     2.  1- real matrix amplitude. 
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     [g g -> h]  Born Effective field theory
c--------------------------------------------------------------------o
       function Born_gg2h(k,p1,p2,p3)
         implicit double precision (a-h,o-z)
         dimension p1(0:3),p2(0:3),p3(0:3)
         parameter(PI=3.141592653589793238D0)
         common/usedalpha/AL,ge
         common/param2/xmur
         common/amass/am1,am2,amh,am4,am5
         common/energy/s

         s12 = 2d0*dot(p1,p2)
         rp34  = dsqrt(s12)
         
c         IF(k .eq. 1)  CF = -4d0/3d0               !leg 1 reduced born k=1 
c         IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2 reduced born k=2

         IF(k .eq. 0)  CF =  1d0               !Leading Order k=0
         IF(k .eq. 1)  CF = -1d0               !leg 1 reduced born k=1 
         IF(k .eq. 2)  CF = -1d0               !Leg 2 reduced born k=2

           NA = 8
c          AS = alphasPDF(amh)
c          AS = AL/4d0/PI

            v = 246d0
c           AL = DSQRT(AL*4d0*PI)
           AS = AL/4d0/PI
           ch = -4d0*AS/3d0/v 
          ch2 = ch * ch
           avg_pol = 4d0
c          Born_gg2h= 16d0*AL**2*amh**4d0/(3d0*PI*v)**2/NA
c           Born_gg2h= AS**2/72d0/PI/v**2/8d0

c 3manually.frm :amp=  1/36*PI^-4*v^-2*s12^2*AL^2    .or.  amp= 1d0/2d0*s12**2*ch2
c          Born_gg2H =  CF*0.5d0*s12**2*ch2/16d0
c          Born_gg2H =  CF*0.5d0*s12**2*ch2/avg_pol/8d0
          Born_gg2H =  CF*s12**2*ch2/64d0
! For gluon average of spin polarization will be 1/(d+epsilon)**2 
c           Born_gg2H =  CF*0.5d0*s12**2*ch2/16d0
c          Born_gg2H = CF * 0.5d0*amh**4*ch2/16d0

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


c       qe=DSQRT(ge*4.d0*PI)
cc      gs=DSQRT(AL*4.d0*PI)
c       qu =1d0! 2d0/3d0
c
c       aem = ge
c       model = 1
c       lambda = 0.226d0 
c       rp34  = dsqrt(P34)
c
c       N=3
c       Cf=(N*N-1.0D0)/2.D0/N
c       Tf=1.0D0/2.0D0
c       N2m1=N*N-1.0D0
c       e=DSQRT(4.0D0*PI*aem)
c
c       nf = 5
c       xnf=nf
c
c       as = AL/4.0D0/PI
           
c
c          print*,"amH:",dsqrt(dot(p3,p3))
c           AS = alphasPDF(amh)
c           AS = AS/4d0/PI

               v = 246d0
               AS = AL/4d0/PI
              ch = -4d0*AS/3d0/v
c              ch = -4d0*AS/3d0/v - 44d0*AS**2/3d0/v 
c              ch =  - 44d0*AS**2/3d0/v 
              ch2 = ch*ch

c             ch2 = ch * ch + (352d0*as**3)/9d0/v/v
c             ch2 = 16d0*AS**2/9d0 + (352d0*AS**3)/9d0
c             ch2 =  (352d0*AS**3)/9d0
              gs = DSQRT(AL*4.d0*PI)
c              CH = -4d0*AS/3d0*(1d0 + AS*(11d0 + 1d0/6d0))
c             ch2 = ch**2/v/v
c              gs = AL 
          gs2CH2 = gs**2*ch2

           s12 = 2d0*dot(p1,p2)
           t11 = 2d0*dot(p1,p4)
           t21 = 2d0*dot(p2,p4)

c          sig = (48d0*gs2CH2*(s12**4 - t11**4 + t11**2*t21**2 - t21**4- 
c     -      s12**3*(t11 + t21) + 
c     -      2d0*s12**2*(t11**2 + 3d0*t11*t21 + t21**2) + 
c     -      s12*(t11**3 + 4d0*t11**2*t21 + 4d0*t11*t21**2 + t21**3)))/
c     -      (s12*t11*t21)
c	print*,"sig1:",sig

         sig = (48d0*gs2CH2*(s12**4 - 2d0*s12**3*(t11 + t21) + 
     -        3d0*s12**2*(t11 + t21)**2 - 2d0*s12*(t11 + t21)**3 + 
     -        (t11**2 + t11*t21 + t21**2)**2))/(s12*t11*t21)/256d0
c	print*,"sig2:",sig

c
c         sig = (24d0*gs2CH2*(2d0*s12**4*t11*t21-
c     -          3d0*t11**3*t21+2d0*t21**2+
c     -      t11**2*(2d0-4d0*t21**2)+t11*(4d0*t21 - 3d0*t21**3d0) - 
c     -    s12**3*(-t11+t11**3 - t21 + 2d0*t11**2*t21 + 2d0*t11*t21**2 + 
c     -         t21**3) + s12**2*
c     -       (2d0 + 5d0*t11**3*t21 + 10d0*t11**2*t21**2 + 
c     -         t11*t21*(-2d0 + 5d0*t21**2)) + 
c     -      s12*(t11 + t21 + 10d0*t11*t21**2 + 3d0*t21**3 + 
c     -       t11**3*(3d0-4d0*t21**2)+t11**2*(10d0*t21-4d0*t21**3))))/s12
c
      return
      end
c---------------------------------------------------------------------
