      REAL FUNCTION ENCAL(IDET,ICHA)
C
C===    Returns the full energy of a calorimeter channel
C===                 IDET - detector number
C===                 ICHA < 1 - all the channels
C===                      >=1 - the channel number
C
      IMPLICIT NONE
      INTEGER  IDET,ICHA
      INCLUDE ?
C
      INTEGER ihi,k,idi,idg1,idg2,iddd
C
C     -----------------------------------------------------------------
C
      ENCAL=0.
C
C
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
C
      IF(npdig.GT.0) THEN
C
         DO idi=1,ndig
            k=ip1dig(idi)
            iddd=IBITS(k,0,16)
            IF(IDET.EQ.0.OR.iddd.EQ.IDET) THEN
               k=ip2dig(idi)
               idg1=IBITS(k,0,16)
               idg2=IBITS(k,16,16)
               IF(ICHA.EQ.0.OR.idg1.EQ.ICHA) THEN
                  ENCAL=ENCAL+REAL(idg2)
               ENDIF
            ENDIF
         ENDDO
C
      ENDIF
C
C
 999  CONTINUE
      END



