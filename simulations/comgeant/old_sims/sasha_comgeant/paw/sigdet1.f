      REAL FUNCTION SIGDET1(IDET,ICH,ITR)
C
C===    Returns the PHIT in IDET,ICH
C===            for the track ITR>0 or full if ITR=0 
C
      IMPLICIT NONE
      INTEGER  IDET,ICH,ITR
      INCLUDE ?
C
      INTEGER i,k,idi,i1,i2,ihi,iddd,itrh,iorh,idg1,idg2
C
C     -----------------------------------------------------------------
C
      SIGDET1=0.
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
C
      IF(npdig.GT.0) THEN
         DO idi=1,ndig
            k=ip1dig(idi)
            iddd=IBITS(k,0,16)
            IF(iddd.EQ.IDET) THEN
               k=ip2dig(idi)
               idg1=IBITS(k,0,16)
               idg2=IBITS(k,16,16)
               IF(idg1.EQ.ICH) THEN
                  i1=1
                  IF(idi.GT.1) THEN
                     k=ip1dig(idi-1)
                     i1=IBITS(k,16,16)
                  ENDIF
                  k=ip1dig(idi)
                  i2=IBITS(k,16,16)
                  DO ihi=i1,i2-1
                     k=ip1hit(ihi)
                     itrh=IBITS(k,0,16)
                     iorh=IBITS(k,16,16)
C
                     IF((itrh.EQ.ITR.AND.ITR.GT.0).OR.ITR.EQ.0) THEN
C                        WRITE(6,*) 
C     +                        IDNEVT,nfind,idi,ihi,itrh,iorh,i1,i2
                        SIGDET1=SIGDET1+phit(ihi)
                     ENDIF
C
                  END DO
               ENDIF
            ENDIF
         END DO
      ENDIF
C
 999  RETURN
C
      END



