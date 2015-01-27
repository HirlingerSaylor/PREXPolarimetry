      REAL FUNCTION IHITR(IDET,ICH,ITR)
C
C===    Returns the hit number of the track ITR (the original)
C===            in the counter IDET, channel ICH
C===            =0 - no hit
C
      IMPLICIT NONE
      INTEGER  IDET,ICH,ITR
      INCLUDE ?
C
      INTEGER i,k,idi,i1,i2,ihi,iddd,itrh,iorh,nfind,idg1,idg2
C
C     -----------------------------------------------------------------
C
      IHITR=0
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      nfind=0
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
                     IF(itrh.EQ.ITR.AND.iorh.EQ.0) THEN
C                        WRITE(6,*) idi,ihi,itrh,iorh,i1,i2
                        nfind=nfind+1
                        IHITR=ihi
                     ENDIF
C
                  END DO
               ENDIF
            ENDIF
         END DO
      ENDIF
C
      IF(nfind.GT.1) WRITE(6,*) ' == Too many hits found=',nfind
C
 999  RETURN
C
      END



