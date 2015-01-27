      REAL FUNCTION SIGLG(IDET,N1,N2)
C
C===    Returns a sum of signals in detector IDET , channels N1-N2 
C
      IMPLICIT NONE
      INTEGER  IDET,N1,N2
      INCLUDE ?
C
      INTEGER i1,i2,iddd,idg1,idg2,i,k
      REAL aa,bb
C
C     -----------------------------------------------------------------
C
      SIGLG=0.
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
         IF(iddd.EQ.IDET.AND.idg1.GE.N1.AND.idg1.LE.N2) 
     +          SIGLG=SIGLG+REAL(idg2)
C     
      END DO
C      CALL RANNOR(aa,bb)
C      SIGLG=SIGLG*0.5
C      IF(SIGLG=SIGLG*(1.+aa*0.05)
C
 999  RETURN
C
      END

