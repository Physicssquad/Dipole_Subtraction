CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO RETURN STRONG COUPLING ALPHA_S=G^2/4/PI         C
C  ---->  WITH INTERPOLATION ACROSS THRESHOLDS  <----          C
C  INPUTS:                                                     C
C         IALPHAS                                              C
C                 = 0: CONSTANT (0.1)                          C
C                 = 1: ONE LOOP                                C
C                 = 2: TWO LOOP                                C
C         MU2                                                  C
C                 = RENORMALIZATION SCALE SQUARED              C
C         LQCD5                                                C
C                 = QCD SCALE FOR NF=5                         C
C  OUTPUT:                                                     C
C         DOUBLE PRECISION ALPHAS                              C
C         NF = EFFECTIVE NUMBER OF FLAVORS                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION ALPHAS(IALPHAS,XMU2,XLQCD5,NEFF)
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z )
      ALPHAS = 0.1D0
      NEFF=5
      IF( IALPHAS.EQ.0 )RETURN
      ALPHAS = ALFAS5(XMU2,XLQCD5,IALPHAS,NEFF)
      RETURN
      END
C.------------------------------------------------------------.C
C. FUNCTION TO RETURN ALPHA_S GIVEN LAMBDA_QCD^5 AND SCALE Q2 .C
C. AND ILOOP=1 (ONE LOOP) OR ILOOP=2 (TWO LOOP)               .C
C. SEE G. ALTARELLI ET AL. NUCL. PHYS. B308 (1988) 724.       .C
C. EQNS (7)->(13)                                             .C
C.------------------------------------------------------------.C
      FUNCTION ALFAS5(Q2,LQCD5,ILOOP,NEFF)
      IMPLICIT DOUBLE PRECISION(A-Z)
      INTEGER ILOOP,NEFF
      A  = 1.0D0
      MB = 4.75D0
      MC = 1.5D0
      AMB=A*MB
      AMC=A*MC
      MU=DSQRT(Q2)
      LQCD52=LQCD5*LQCD5
      IF      (MU.GT.AMB) THEN
         ALFAS5 = ASLOOP(Q2,LQCD52,5,ILOOP)
C. MUST RETURN NF=5 FOR RENORMALIZATION GROUP TO WORK
         NEFF=5
      ELSE IF (MU.GT.AMC) THEN
         AMB2=AMB*AMB
         ALFAS5 = 1D0 / ( 1D0 / ASLOOP(Q2  ,LQCD52,4,ILOOP) + 
     &                    1D0 / ASLOOP(AMB2,LQCD52,5,ILOOP) -
     &                    1D0 / ASLOOP(AMB2,LQCD52,4,ILOOP) )
C. MUST RETURN NF=4 FOR RENORMALIZATION GROUP TO WORK
         NEFF=4
      ELSE
         AMB2=AMB*AMB
         AMC2=AMC*AMC
         ALFAS5 = 1D0 / ( 1D0 / ASLOOP(Q2  ,LQCD52,3,ILOOP) + 
     &                    1D0 / ASLOOP(AMC2,LQCD52,4,ILOOP) +
     &                    1D0 / ASLOOP(AMB2,LQCD52,5,ILOOP) -
     &                    1D0 / ASLOOP(AMB2,LQCD52,4,ILOOP) - 
     &                    1D0 / ASLOOP(AMC2,LQCD52,3,ILOOP) )
C. MUST RETURN NF=3 FOR RENORMALIZATION GROUP TO WORK
         NEFF=3
      ENDIF
      RETURN
      END
C.---------------------------------------------------------------.C
C. FUNCTION TO RETURN ALPHAS GIVEN MU2,LQCD2,NF AND NORDER       .C
C.---------------------------------------------------------------.C
      FUNCTION ASLOOP(XMU2,XLQCD2,NF,NORDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.14159265358979D0)
      IF    ( NORDER.EQ.1 ) THEN
         B1=(33D0-2D0*NF)/(12D0*PI)
         T=XMU2/XLQCD2
         ASLOOP=1D0/(B1*DLOG(T))
      ELSEIF( NORDER.EQ.2 ) THEN
         B1=(33D0-2D0*NF)/(12D0*PI)
         B2=(153D0-19D0*NF)/(2D0*PI*(33D0-2D0*NF))
         T=XMU2/XLQCD2
         F1=B1*DLOG(T)
         ASLOOP=(1D0-B2*DLOG(DLOG(T))/F1)/F1
      ELSE
         WRITE(*,*) 'ERROR IN ASLOOP: NORDER = ',NORDER
         STOP
      ENDIF
      RETURN
      END
