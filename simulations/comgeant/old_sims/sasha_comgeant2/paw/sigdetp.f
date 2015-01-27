      REAL FUNCTION SIGDETP(IDET,ICH,ITR,IPAR,ID)
C
C===    Returns the PHIT in IDET,ICH
C===            for the track ITR>0 or full if ITR=0 (only the original track)
C===            Fills the histogram ID with all tracks originated from ITRA of the type IPAR:
C===            only the original track if IPAR=0 
C===            or for particles IPAR>0 (ex: IPAR=1/2/3 - gamme/electrons/positrons)
C===
C
      IMPLICIT NONE
      INTEGER  IDET,ICH,ITR,IPAR,ID
      INCLUDE ?
C
      LOGICAL  HEXIST
      EXTERNAL HEXIST
C
      INTEGER i,k,idi,i1,i2,ihi,iddd,itrh,iorh,idg1,idg2,ifirst
      REAL wg
      DATA ifirst/1/
C
C     -----------------------------------------------------------------
C
      SIGDETP=0.
C
      IF(ifirst.EQ.1) THEN
         IF(ID.NE.0.AND.HEXIST(id)) THEN
            CALL HRESET(ID,'    ')
         ELSE
            WRITE(6,*) ' *** ERROR: no histogram ID=',ID
         ENDIF
      ENDIF
      ifirst=0
C
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
C
      wg=1.
      IF(WGCROS.GT.0.) wg=WGCROS
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
                     iorh=IBITS(k,16,16)  ! =0 - the original particle, =1 - gamma etc...
C
                     IF((itrh.EQ.ITR.AND.ITR.GT.0).OR.ITR.EQ.0) THEN
C                        WRITE(6,*) 
C     +                        IDNEVT,nfind,idi,ihi,itrh,iorh,i1,i2
                        SIGDETP=SIGDETP+phit(ihi)
C
                        IF(iorh.EQ.0.OR.iorh.EQ.IPAR) THEN
                           CALL HF1(ID,phit(ihi),wg)
                        ENDIF
C
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



