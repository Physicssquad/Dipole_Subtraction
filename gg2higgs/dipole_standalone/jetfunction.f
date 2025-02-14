c_______________________________________________________c
c... jet function for m+1 (real radiation) 
      SUBROUTINE Fjm_plus_one(p1,p2,p3,p4,jet_cut)
      implicit double precision(a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      common/distribution/rapidity,p_T
      eps = 50d0
      obsvl = p_T
      xlow = obsvl - eps
      xhigh= obsvl + eps
      amH = 125d0 !... in Gev

      s12 = 2d0*dot(p1,p2) 
      s13 = 2d0*dot(p1,p3) 
      s23 = 2d0*dot(p2,p3) 
      s14 = 2d0*dot(p1,p4) 
      s24 = 2d0*dot(p2,p4) 

        jet_cut = 0

      xjet_pt = dmax1(higgs_pt, gluon_pt)

c      if (xjet_pt .ge. xlow .and. xjet_pt .le. xhigh ) then
c        jet_cut = 1
c      else
c      endif

      return
      end

c_______________________________________________________c
c... jet function for dipoles
      SUBROUTINE Fjm(ptilde,jet_cut)
      implicit double precision(a-h,o-z)
      dimension ptilde(0:3,3,2)
      common/distribution/rapidity,p_T
      eps = 50d0
      obsvl = p_T
c      xlow = obsvl - eps
      xlow = 0d0
      xhigh= obsvl + eps

      s12 = 2d0*dot2(ptilde,1,2,1) 
      s13 = 2d0*dot2(ptilde,1,3,1) 
      s23 = 2d0*dot2(ptilde,2,3,1) 
c      print*,"m-body",s13/2d0 + s23/2d0
 

      s12 = 2d0*dot2(ptilde,1,2,1) 
      s14 = 2d0*dot2(ptilde,1,4,1) 
      s24 = 2d0*dot2(ptilde,2,4,1) 
c      print*,"m-body",s14 + s24

      xjet_pt = dmax1(higgs_pt1, higgs_pt2)

      if (xjet_pt .ge. xlow .and. xjet_pt .le. xhigh ) then
        jet_cut = 1
      else
        jet_cut = 0
      endif

      return
      end
c_______________________________________________________c



































































































cc_______________________________________________________c
cc... jet function for m+1 (real radiation) 
c      SUBROUTINE Fjm_plus_one(p1,p2,p3,p4,jet_cut)
c      implicit double precision(a-h,o-z)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c      common/distribution/rapidity,p_T
c      eps = 0.5d0
c      obsvl = p_T
c      xlow = obsvl - eps
c      xhigh= obsvl + eps
c
c      higgs_pt = dsqrt(p3(1)**2 + p3(2)**2)
c      if (higgs_pt .ge. xlow .and. higgs_pt .le. xhigh ) then
c        jet_cut = 1
c      else
c        jet_cut = 0
c      endif
c
c      return
c      end
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
cc... jet function for dipoles
c      SUBROUTINE Fjm(ptilde,jet_cut)
c      implicit double precision(a-h,o-z)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptilde(0:3,3,2)
c      common/distribution/rapidity,p_T
c      eps = 0.5d0
c      obsvl = p_T
c      xlow = obsvl - eps
c      xhigh= obsvl + eps
c
c      higgs_pt1 = dsqrt(ptilde(1,3,1)**2 + ptilde(2,3,1)**2)
c      higgs_pt2 = dsqrt(ptilde(1,3,2)**2 + ptilde(2,3,2)**2)
c
c      if (higgs_pt .ge. xlow .and. higgs_pt .le. xhigh ) then
c        jet_cut = 1
c      else
c        jet_cut = 0
c      endif
c
c      return
c      end
c
