      REAL FUNCTION T_MOLL(X)
C
C---      Various Moller calculations
C ===     V_MOLL(1)=0. : X-theta CM, output - COS(Theta Lab) 
C ===     V_MOLL(2)=E beam GeV
C
      IMPLICIT NONE
      REAL X
      VECTOR V_MOLL(10)
C
      REAL am,e,gam,bet,thetcm,thetlab,ecm,ctl,stl,ctcm,stcm,ep,pp
C
      am=0.511
      e=V_MOLL(2)*1000.
      ecm=SQRT(2.*am**2+2.*am*e)
      gam=(e+am)/ecm
      bet=SQRT(1.-1./gam**2)
      ep=ecm/2.
      pp=SQRT(ep**2-am**2)
C
      IF(ABS(V_MOLL(1)-0.).LT.0.1) THEN
         thetcm=X*3.14159/180.
         ctcm=COS(thetcm)
         stcm=SIN(thetcm)
         ctl=pp*stcm/gam/(pp*ctcm+bet*ep)
         T_MOLL=ctl
      ENDIF
C
      END
