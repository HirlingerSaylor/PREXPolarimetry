      REAL FUNCTION SOLEN_P(XX)
*
* ===   Plots the B solenoid field along X at Z defined in the code
*
      IMPLICIT NONE
      REAL XX
      REAL SOLEN_V
      COMMON/CSOLEN/ PARSOL(12,3),BSOL(3)
      REAL           PARSOL,BSOL
C
      REAL yy,zz,b
C
C     -----------------------------------------------------------------
C
      yy=0.
      zz=0.
C
      b=SOLEN_V(XX,yy,zz)
      SOLEN_P=BSOL(3)
C
      RETURN
      END
C
      INCLUDE 'solen_v.f'

