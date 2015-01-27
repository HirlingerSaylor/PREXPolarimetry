      REAL FUNCTION DAPNUM(ICOUNT,ICH)
C
C===    Returns the number of digitizings in a counter
C===    ICOUNT - the counter (=1 - apperture)
C===    ICH - the number of the channel (0-all, or 1 or 2)
C===          If ICH>0 the return is the sum of the track numbers which 
C===                   contributed to the digitizing
C
C

      IMPLICIT NONE
      INTEGER  ICOUNT,ICH
      INCLUDE ?
C
      INTEGER idetlg,i1,i2,iddd,idg1,idg2,i,k
C
C     -----------------------------------------------------------------
C
      DAPNUM=0.
      idetlg=ICOUNT
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
            IF(ICH.EQ.0.OR.ICH.EQ.idg1) THEN
               DAPNUM=DAPNUM+1.
            ENDIF
         ENDIF
C     
      END DO

C
 999  RETURN
C
      END



