      REAL FUNCTION HITRAXY(IXY,IDET,ICH,ITR)
C
C===    Returns the X/Y-coordinate MRS/DRS of the hit originated from
C===            IXY=1,2  - X,Y DRS
C                  -1,-2 - X,Y MRS
C                  -3    - Z MRS
C                  -4    - P mom
C                  -5    - Slope X DRS
C                  -6    - Slope Y DRS
C===            the track ITR (the original)
C===            in the counter IDET, channel ICH
C===            ICH=0 -  1-st hit
C
      IMPLICIT NONE
      INTEGER  IXY,IDET,ICH,ITR
      INCLUDE ?
C
      INTEGER i,k,idi,i1,i2,ihi,iddd,itrh,iorh,nfind,idg1,idg2,j
C
C     -----------------------------------------------------------------
C
      HITRAXY=-100.
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      IF(IABS(IXY).LT.1.OR.IABS(IXY).GT.6) GO TO 999
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
C                     write(6,*) idi,iddd,idg1
               IF(idg1.EQ.ICH.OR.ICH.EQ.0) THEN
                  i1=1
                  IF(idi.GT.1) THEN
                     k=ip1dig(idi-1)
                     i1=IBITS(k,16,16)
                  ENDIF
                  k=ip1dig(idi)
                  i2=IBITS(k,16,16)
                  DO j=i1,i2-1
                     ihi=jpdig(j)
                     k=ip1hit(ihi)
                     itrh=IBITS(k,0,16)
                     iorh=IBITS(k,16,16)
C
                     IF(itrh.EQ.ITR.AND.iorh.EQ.0) THEN
                        nfind=nfind+1
C                        WRITE(6,*) 
C     +                        IDNEVT,nfind,idi,ihi,itrh,iorh,i1,i2
                        IF(IXY.GT.0) THEN
                           HITRAXY=hitd(IXY,ihi)
                        ELSE 
                           IF(IABS(IXY).LT.3) THEN
                              HITRAXY=hit(-IXY,ihi)
                           ELSE IF(IABS(IXY).EQ.3) THEN
                              HITRAXY=zhit(ihi)
                           ELSE IF(IABS(IXY).EQ.4) THEN
                              HITRAXY=phit(ihi)
                           ELSE IF(IABS(IXY).EQ.5) THEN
                              HITRAXY=hitsld(1,ihi)
                           ELSE IF(IABS(IXY).EQ.6) THEN
                              HITRAXY=hitsld(2,ihi)
                           ENDIF
                        ENDIF
                     ENDIF
C
                  END DO
               ENDIF
            ENDIF
         END DO
      ENDIF
C
      IF(nfind.GT.1) WRITE(6,*) ' == ev=',IDNEVT
     +              ,' Too many hits found=',nfind
     +              ,' IXY,IDET,ICH,ITR=',IXY,IDET,ICH,ITR
C
 999  RETURN
C
      END


