c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
      subroutine getPKPlus(iplus,x,xmuf,p1,p2,p3,SumPlus)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common /usedalpha/ AL,ge 
c~~~~[ These functions are taken from misc.f ]~~~~~~~~c
      external PggP,PggReg,PggDel
      external AKbarP_gg,AKbarReg_gg,AKbarD_gg
      external AKtilP_gg,AKtilreg_gg,AKtilD_gg
      external aKbar_gq,aKtil_gq,Pgq_reg

            s12 = 2d0*dot(p1,p2)
             Cf = 4d0/3d0                      
             Tr = 0.5d0

          xmuf2 = xmuf*xmuf

          Pplus = 0.0d0
       SumPlus  = 0.0d0

c... s12 = x * [2*p1 \times p2]

          Pplus = PggP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! kinematics depends on the PS generation 
        SumPlus = Pplus + AKbarP_gg(x) + AKtilP_gg(x)
c	print*,"Plus:",Pplus,AKbarP_gg(x),AKtilP_gg(x),x
c	SumPlus = PqqP(x)
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      subroutine getPKReg(x,xmuf,p1,p2,p3,SumReg)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(0:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
      common /usedalpha/ AL,ge 

      external PggPlus,PggReg,PggDel
      external AKbarP_gg,AKbarReg_gg,AKbarD_gg
      external AKtilP_gg,AKtilreg_gg,AKtilD_gg
      external aKbar_gq,aKtil_gq,Pgq_reg



      xmuf2 = xmuf*xmuf
      s12   = 2d0*dot(p1,p2)

       Preg = PggReg(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
       Areg =  PReg + AKbarReg_gg(x) + AKtilReg_gg(x)      

       SumReg = Areg
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      subroutine getPKDel(x,xmuf,p1,p2,p3,SumDel)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
c      common/scales/xmuf,xmur
      external Born_gg2h

      external PggPlus,PggReg,PggDel
      external AKbarP_gg,AKbarReg_gg,AKbarD_gg
      external AKtilP_gg,AKtilreg_gg,AKtilD_gg
      external aKbar_gq,aKtil_gq,Pgq_reg

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
        xmuf2 = xmuf*xmuf
         xmur = xmuf
        AL = alphasPDF(xmur)

      do k = 1,2
        Born = Born_gg2h_(0,AL,p1,p2,p3)
        coef = Born
       xmuf2 = xmuf*xmuf

        Pdel = PggD(x)*(-1d0)*dlog(xmuf2/s12)    ! here x=1
      AllP(k)= ( Pdel+AKbarD_gg(x)+AKtilD_gg(x) )*coef
      enddo

       SumDel = AllP(1)+AllP(2)
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[** END **]
