c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
      subroutine getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common /usedalpha/ AL,ge 
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg 

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0                      
        Tr = 0.5d0
        Alp = Al/2.0d0/pi
        Alp = 1.0d0
        xmuf2 = xmuf*xmuf

          Pplus = 0.0d0
       SumPlus  = 0.0d0
         SumReg = 0.0d0
         SumDel = 0.0d0

        if ( iplus .eq. 1 ) then
          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)
c          write(*,*)'iplus, s12 =', iplus, s12
        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)

        elseif( iplus .eq. 0 ) then
c          write(*,*)'iplus, s12 =', iplus, s12
          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)
        endif

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
      external Born_uU2eE
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
        xmuf2 = xmuf*xmuf

        Born = Born_uU2eE(0,p1,p2,p3,p4)
        coef = Born

       Areg = (AKbarreg_qq(x)+AKtilreg_qq(x))*coef

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
      external Born_uU2eE
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
c        Alp = Al/2.0d0/pi
        Alp = 1.0d0
        xmuf2 = xmuf*xmuf

      do k = 1,2

        Born = Born_uU2eE(0,p1,p2,p3,p4)
        coef = Born

       AllP(k)= (AKbarD_qq(x)+AKtilD_qq(x))*coef
      enddo

       SumDel = AllP(1)+AllP(2)

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
