c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
      function dot(p,q)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end

      function dot2(p,i,j,l)
      implicit double precision (a-h,o-z)
      dimension p(0:3,3,2)
      dot2= p(0,i,l)*p(0,j,l)-p(1,i,l)*p(1,j,l)
     .     -p(2,i,l)*p(2,j,l)-p(3,i,l)*p(3,j,l)
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Functions used in PK_Plus PK_Delta and PK_Regular]

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ For qq initiated process with qg splliting ]
      double precision function PqqP(x)
      implicit double precision (a-h,o-z)
      data PI/3.141592653589793238462643D0/

      CF = 4.0d0/3.0d0
c      PqqP = CF*(1.0d0+x*x)/(1.0d0-x)

      XM=1.0D0-X
      DLOGXM= DLOG(XM)
      XD1=DLOGXM/XM
      CQQB1PLUS=16.0D0*XD1
      PqqP=CF*CQQB1PLUS
      return
      end
c________________________________________________________________________c 
      double precision function Pqqreg(x)
      implicit double precision (a-h,o-z)
      CF = 4.0d0/3.0d0
      pqqreg = -CF*(1.0d0+x)
      return
      end
c________________________________________________________________________c 

      double precision function AKbarP_qq(x)
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarP_qq = CF*2.0d0/(1.0d0-x)*dlgx
      return
      end
c________________________________________________________________________c 

      double precision function AKbarreg_qq(x)
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarreg_qq = -(1.0d0+x)*dlgx + (1-x)
      AKbarreg_qq = CF*AKbarreg_qq
      return
      end
c________________________________________________________________________c 

      double precision function AKbarD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      CF = 4.0d0/3.0d0
      AKbarD_qq = -CF*(5.0d0-pi*pi)
      return
      end
c________________________________________________________________________c 

      double precision function AKtilP_qq(x)
      implicit double precision (a-h,o-z)

      CF = 4.0d0/3.0d0
      dlgx = dlog(1.0d0-x)
      AKtilP_qq = CF*2.0d0/(1.0d0-x)*dlgx
      return
      end
c________________________________________________________________________c 


      double precision function AKtilreg_qq(x)
      implicit double precision (a-h,o-z)
      external Pqqreg

      dlgx = dlog(1.0d0-x)
      AKtilreg_qq = Pqqreg(x)*dlgx
      return
      end
c________________________________________________________________________c 

      double precision function AKtilD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      CF = 4.0d0/3.0d0

      AKtilD_qq = -CF*pi*pi/3.0d0
      return
      end
c________________________________________________________________________c 


      double precision function Pgq_reg(x) !same as Pgqb(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      Tr = 0.5d0
      Pgq_reg = Tr*(x**2 + (1-x)**2)

      return
      end
c________________________________________________________________________c 

      double precision function aKbar_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      Tr = 0.5d0
      dlgx = dlog((1d0-x)/x)
      aKbar_gq = Pgq_reg(x)*dlgx +Tr*2d0*x*(1d0-x)
      return
      end
c________________________________________________________________________c 


      double precision function aKtil_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      dlgx = dlog(1d0 - x)
      aKtil_gq = Pgq_reg(x)*dlgx

      return
      end
c________________________________________________________________________c 

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ For gg initiated process with gg splliting ]

      double precision function PggP(x)                            !      [---> P_gg_Plus]
      implicit double precision (a-h,o-z)
      data PI/3.141592653589793238462643D0/

      CA = 3d0
      one_minus_x = (1.0d0-x)
      PggP = 2d0*CA*1d0/one_minus_x

      return
      end
c________________________________________________________________________c 

      double precision function PggReg(x)                    !            [--->P_gg_Regular terms]
      implicit double precision (a-h,o-z)

      CF = 4.0d0/3.0d0
      CA = 3d0

      Pggreg = (1d0-x)/x -1d0 +x*(1d0-x) 
      Pggreg = 2d0*CA*Pggreg   

      return
      end
c________________________________________________________________________c 

      double precision function PggD(x)                    !              [--->P_gg_delta terms]
      implicit double precision (a-h,o-z)

      CA = 3d0
      TR = 0.5d0
      NF = 5d0   !...[?] number of quark flavours?

      PggD = 11d0/6d0*CA - 2d0/3d0*NF*TR 
      return
      end
c________________________________________________________________________c 
c________________________________________________________________________c 

      double precision function AKbarP_gg(x)            !                  [----> K_Bar_Plus terms  for gg channel
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CA = 3d0
      AKbarP_gg = 2d0*CA*1d0/(1d0-x)*dlgx

      return
      end
c________________________________________________________________________c 

      double precision function AKbarReg_gg(x)            !               [----> K_Bar_Regular terms  for gg channel
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CA = 3d0
      AKbarReg_gg = 2d0*CA*((1d0-x)/x-1d0+x*(1d0-x))*dlgx
      return
      end
c________________________________________________________________________c 

      double precision function AKbarD_gg(x)            !                 [----> K_Bar_Delta terms  for gg channel
      implicit double precision (a-h,o-z)

      CA = 3d0
      TR = 0.5d0
      NF = 5d0  !...[?] is this no of quark flavours?
c...for checking I am taking 5 flavour as in the matrix.[11 feb]

      AKbarD_gg = (50d0/9d0-Pi**2)*CA - 16d0/9d0*TR*NF

      return
      end
c________________________________________________________________________c 
c________________________________________________________________________c 

      double precision function AKtilP_gg(x)
      implicit double precision (a-h,o-z)

      CA = 3.0d0
      dlgx = dlog(1.0d0-x)
      AKtilP_gg = CA*2.0d0/(1.0d0-x)*dlgx
      return
      end
c________________________________________________________________________c 


      double precision function AKtilreg_gg(x)                  !         [----> K_Tilde_Regular terms  for gg channel
      implicit double precision (a-h,o-z)

      CA = 3.0d0

      dlgx = dlog(1.0d0-x)
      AKtilreg_gg = 2d0*CA*((1d0-x)/x - 1d0 + x*(1d0-x))
      AKtilreg_gg = AKtilreg_gg * dlgx 

      return
      end
c________________________________________________________________________c 

      double precision function AKtilD_gg(x)                   !         [----> K_Tilde_Delta terms  for gg channel
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      CA = 3.0d0

      AKtilD_gg = -CA*pi*pi/3.0d0

      return
      end
c________________________________________________________________________c 


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[End of PK functions]
      SUBROUTINE MATDIP_CHECK(RESULT, SIG, SUMD, MODE)
      DOUBLE PRECISION SIG, SUMD
      CHARACTER*(*) MODE
      INTEGER RESULT

      ! Initialize RESULT to 0 by default
      RESULT = 0

      ! Check for 'flip' mode
      IF (MODE .EQ. 'flip') THEN
          IF ((SIG .GT. 0D0 .AND. SUMD .LT. 0D0) .OR.
     &        (SIG .LT. 0D0 .AND. SUMD .GT. 0D0)) THEN
              RESULT = 1
          ENDIF

      ! Check for 'noflip' mode
      ELSE IF (MODE .EQ. 'noflip') THEN
          IF ((SIG .GT. 0D0 .AND. SUMD .GT. 0D0) .OR.
     &        (SIG .LT. 0D0 .AND. SUMD .LT. 0D0)) THEN
              RESULT = 1
          ENDIF
      ENDIF

      RETURN
      END
      
