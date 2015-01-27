      REAL FUNCTION HITPAR(IDET,IPAR,ID)
C
C===    Fills the histogram ID with the energies of the particles type IPAR
C===              hitting the detector IDET
C
      IMPLICIT NONE
      INTEGER  IDET,IPAR,ID
      INCLUDE ?
C
      INTEGER ihi,nfind,k,itrh,iorh,iddh,itmh
C
C     -----------------------------------------------------------------
C
      HITPAR=0.
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      nfind=0
C
      IF(npdig.GT.0) THEN
         DO ihi=1,nhit
C
            k=ip1hit(ihi)
            itrh=IBITS(k,0,16)
            iorh=IBITS(k,16,16)
            k=ip2hit(ihi)
            iddh=IBITS(k,0,16)
            itmh=IBITS(k,16,16)-32768
C
            IF(itra(itrh).EQ.IPAR.AND.iddh.EQ.IDET) THEN
               nfind=nfind+1
               HITPAR=HITPAR+1.
               CALL HF1(ID,phit(ihi),1.)
C               WRITE(6,*) HITPAR,phit(ihi),itrh,itra(itrh)
            ENDIF
         ENDDO
      ENDIF
C
 999  CONTINUE
      END
