c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      subroutine getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
      common /usedalpha/ AL,ge 
      external Born_gg2aa
      external PggP,Pggreg
      external AKbarP_gg,AKbarreg_gg,AKbarD_gg
      external AKtilP_gg,AKtilreg_gg,AKtilD_gg

       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
        xmuf2 = xmuf*xmuf

        Born = Born_gg2aa(0,p1,p2,p3,p4)
        coef = Born

        P_Reg = PggReg(x)*(-1.0d0)*dlog(xmuf2/s12)
       Areg = (P_Reg + AKbarreg_gg(x) + AKtilreg_gg(x) )*coef

       SumReg = Areg
      return
      end
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
c      subroutine getPKDel(x,xmuf,p,xp1,xp2,SumDel)
c      implicit double precision (a-h,o-z)
c      parameter (pi=3.14159265358979d0)
c      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
c      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
c      dimension  xp1(0:3),xp2(0:3)
cc      common /usedalpha/ AL,ge 
c      external Born_uU2eE
c      external PggP,Pggreg
c      external AKbarP_gg,AKbarreg_gg,AKbarD_gg
c      external AKtilP_gg,AKtilreg_gg,AKtilD_gg
c
c       SumP(1) = 0d0
c       SumK(1) = 0d0
c       SumP(2) = 0d0
c       SumK(2) = 0d0
c
c      call p2d_to_p1d_4(p,p1,p2,p3,p4)
c        s12 = 2d0*dot(p1,p2)
c        Cf = 4d0/3d0
c        Tr = 0.5d0
cc        Alp = Al/2.0d0/pi
c        Alp = 1.0d0
c        xmuf2 = xmuf*xmuf
c
c      do k = 1,2
c
c        Born = Born_uU2eE(0,p1,p2,p3,p4)
c        coef = Born
c
c
c         Pdel = PggD(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
c      AllP(k) = (Pdel + AKbarD_gg(x) + AKtilD_gg(x))*coef
c      enddo
c
c       SumDel = AllP(1)+AllP(2)
c
c      return
c      end
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
