      REAL FUNCTION ASYMLG(RESOL)
C
C===    Returns asymmetry : (L-R)/(L+R) 
C===    RESOL - energy resolution
C      
      IMPLICIT NONE
      REAL  RESOL
      INCLUDE ?
C
      INTEGER idetlg,i1,i2,iddd,idg1,idg2,i,k
      REAL aa,bb,sl,sr
C
C     -----------------------------------------------------------------
C
      ASYMLG=0.
      idetlg=2
      sl=0.
      sr=0.
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
C
      DO i=1,ndig
C
C--           Get the digitisations
C	
         i1=0
         i2=0
         IF(npdig.GT.0) THEN
            i1=1
            IF(i.GT.1) THEN
               k=ip1dig(i-1)
               i1=IBITS(k,16,16)
            ENDIF
            k=ip1dig(i)
            i2=IBITS(k,16,16)
         ENDIF
         k=ip1dig(i)
         iddd=IBITS(k,0,16)
         k=ip2dig(i)
         idg1=IBITS(k,0,16)
         idg2=IBITS(k,16,16)
         IF(iddd.EQ.idetlg) THEN
            IF(idg1.GE.1.AND.idg1.LE.4) THEN 
               sl=sl+REAL(idg2)
            ELSE
               sr=sr+REAL(idg2)
            ENDIF
         ENDIF
C     
      END DO
      CALL RANNOR(aa,bb)
      sl=sl+aa*RESOL*SQRT(1000.*sl)
      sr=sr+bb*RESOL*SQRT(1000.*sr)
      ASYMLG=(sl-sr)/(sl+sr)
C
 999  RETURN
C
      END





