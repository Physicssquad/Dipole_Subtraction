c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
      function dot(p,q)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Functions used in PK_Plus PK_Delta and PK_Regular]
      double precision function PqqP(x)
      implicit double precision (a-h,o-z)
      data PI/3.141592653589793238462643D0/

      CF = 4.0d0/3.0d0
      PqqP = CF*(1.0d0+x*x)/(1.0d0-x)

c      XM=1.0D0-X
c      DLOGXM= DLOG(XM)
c      XD1=DLOGXM/XM
c      CQQB1PLUS=16.0D0*XD1
c      PqqP=CF*CQQB1PLUS
      return
      end

      double precision function Pqqreg(x)
      implicit double precision (a-h,o-z)
      CF = 4.0d0/3.0d0
      pqqreg = -CF*(1.0d0+x)
      return
      end


      double precision function AKbarP_qq(x)
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarP_qq = CF*2.0d0/(1.0d0-x)*dlgx
      return
      end


c      double precision function AKbarP_qq(x)
c      implicit double precision (a-h,o-z)
c      CF = 4.0d0/3.0d0
c
cc      dlgx = dlog((1.0d0-x)/x)
c      dlgx1 = dlog((1.0d0-x))
c      dlgx2 = dlog(x)
c      AKbarP_qq =  2d0*CF*( dlgx1/(1d0-x)-dlgx2/(1-x) )
cc      AKbarP_qq = CF*2.0d0/(1.0d0-x)*(dlgx1 - dlgx2) 
cc      AKbarP_qq = CF*2.0d0/(1.0d0-x)*dlgx 
c      return
c      end

      double precision function AKbarreg_qq(x)
      implicit double precision (a-h,o-z)

      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarreg_qq = -(1.0d0+x)*dlgx + (1-x)
      AKbarreg_qq = Cf *AKbarreg_qq 
      return
      end

      double precision function AKbarD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      CF = 4.0d0/3.0d0
      AKbarD_qq = -CF*(5.0d0-pi*pi)
      return
      end

      double precision function AKtilP_qq(x)
      implicit double precision (a-h,o-z)

      CF = 4.0d0/3.0d0
      dlgx = dlog(1.0d0-x)
      AKtilP_qq = CF*2.0d0/(1.0d0-x)*dlgx
      return
      end


      double precision function AKtilreg_qq(x)
      implicit double precision (a-h,o-z)
      external Pqqreg

      dlgx = dlog(1.0d0-x)
      AKtilreg_qq = Pqqreg(x)*dlgx
      return
      end

      double precision function AKtilD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      CF = 4.0d0/3.0d0

      AKtilD_qq = -CF*pi*pi/3.0d0
      return
      end


      double precision function Pgq_reg(x) !same as Pgqb(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      Tr = 0.5d0
      Pgq_reg = Tr*(x**2 + (1-x)**2)

      return
      end

      double precision function aKbar_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      Tr = 0.5d0
      dlgx = dlog((1d0-x)/x)
      aKbar_gq = Pgq_reg(x)*dlgx +Tr*2d0*x*(1d0-x)
      return
      end


      double precision function aKtil_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      dlgx = dlog(1d0 - x)
      aKtil_gq = Pgq_reg(x)*dlgx

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[End if PK functions]
      
