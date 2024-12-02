c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
      subroutine getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common /usedalpha/ AL,ge 

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

          Pplus = PggP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! kinematics depends on the PS generation 
        SumPlus = Pplus + AKbarP_gg(x) + AKtilP_gg(x)
c	SumPlus = PqqP(x)

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      subroutine getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
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
      subroutine getPKDel(x,xmuf,p,xp1,xp2,SumDel)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
c      common /usedalpha/ AL,ge 
      external Born_gg2aa

      external PggPlus,PggReg,PggDel
      external AKbarP_gg,AKbarReg_gg,AKbarD_gg
      external AKtilP_gg,AKtilreg_gg,AKtilD_gg
      external aKbar_gq,aKtil_gq,Pgq_reg



       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
        xmuf2 = xmuf*xmuf
        

      do k = 1,2

        Born = Born_gg2aa(0,p1,p2,p3,p4)
        coef = Born


       xmuf2 = xmuf*xmuf
       s12   = 2d0*dot(p1,p2)

        Pdel = PggD(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1

      AllP(k)= (Pdel + AKbarD_gg(x) + AKtilD_gg(x))*coef
      enddo

       SumDel = AllP(1)+AllP(2)


      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[** END **]
