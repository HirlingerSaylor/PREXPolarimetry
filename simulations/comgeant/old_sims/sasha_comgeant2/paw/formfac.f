      REAL FUNCTION FORMFAC(Q2)
*
*  Formfactor of a nucleus or nucleon
*
      IMPLICIT NONE
      REAL Q2,t
C
      VECTOR PKIN_PSI(10)
      REAL am,qr,rnucl,ge
C
      am=PKIN_PSI(2)
      IF(am.GT.1.5) THEN
         rnucl=1.07*am**0.3333/0.197 ! raduius in GeV**-1, (1.07fm*a**3)
         qr=rnucl*SQRT(ABS(Q2))
         ge=3.*1./qr**3*(SIN(qr)-qr*COS(qr))
C         ge=ge/(1./(1.+Q2/0.71)**2)
         FORMFAC=ge**2
      ELSE
         ge=1./(1.+Q2/0.71)**2
         FORMFAC=ge**2
      ENDIF
      END
