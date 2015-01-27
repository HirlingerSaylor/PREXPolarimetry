      REAL FUNCTION SOLEN_V(XX,YY,ZZ)
*
* ===   Solenoid field
* ===   Transfer the vector parameters to a COMMON and call solen_f.f
*
      IMPLICIT NONE
      REAL XX,YY,ZZ
      REAL SOLEN_F
C
      VECTOR SOL_PAR(12,3)
      VECTOR SOL_B(3)
C
      COMMON/CSOLEN/ PARSOL(12,3),BSOL(3)
      REAL           PARSOL,BSOL
C
      INTEGER i,j
C
C
C     -----------------------------------------------------------------
C
      DO i=1,3
         DO j=1,12
            PARSOL(j,i)=SOL_PAR(j,i)
         ENDDO
      ENDDO
C
      SOLEN_V=SOLEN_F(XX,YY,ZZ)
C
      DO i=1,3
         SOL_B(i)=BSOL(i)
      ENDDO
C
      RETURN
      END
C
      INCLUDE 'solen_f.f'
