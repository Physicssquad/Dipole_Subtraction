c_______________________________________________________c
c... jet function for m+1 (real radiation) 
      subroutine Fjm_plus_one(p1,p2,p3,p4,p5,jet_cut)
      implicit double precision(a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      common/distribution/xq
      eps = 0.5d0 
      xlow = xq - eps
      xhigh= xq + eps
      
      xinvmass = dsqrt(2d0*dot(p3,p4))
      if (xinvmass .ge. xlow .and. xinvmass .le. xhigh ) then
        jet_cut = 1
      else 
        jet_cut = 0
      endif

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c... jet function for dipoles
      subroutine Fjm(p1,p2,p3,p4,jet_cut)
      implicit double precision(a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common/distribution/xq
      eps = 0.5d0 
      xlow = xq - eps
      xhigh= xq + eps
      
      xinvmass = dsqrt(2d0*dot(p3,p4))
      if (xinvmass .ge. xlow .and. xinvmass .le. xhigh ) then
        jet_cut = 1
      else 
        jet_cut = 0
      endif

      return
      end
c_______________________________________________________c


