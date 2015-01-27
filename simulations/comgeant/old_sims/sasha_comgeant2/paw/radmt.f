      REAL FUNCTION RADMT(THET)
      IMPLICIT NONE
C
      REAL THET,ELOSS
      REAL RADMOTT
C
C===        Returns the Rad. Mott cross-section in barn/ster 
C
      VECTOR VRADMOT(4)
C
      COMMON/CRADMOT/ EBEAMRM,ATARGRM,ZTARGRM
      REAL            EBEAMRM    !   - beam energy
     +               ,ATARGRM    !   - target A
     +               ,ZTARGRM    !   - target Z
C
      EBEAMRM=VRADMOT(1)
      ATARGRM=VRADMOT(2)
      ZTARGRM=VRADMOT(3)
      ELOSS=VRADMOT(4)
C
      RADMT=RADMOTT(THET,ELOSS)
C
      END
C
      INCLUDE 'radmott.f'
