      REAL FUNCTION APNUM(IARM,IFL,ITHRESH)
C
C===    Returns the number of digitizings in the apperture counters
C            
C     IARM=0 - in both arms
C          1 - in the left arm
C          2 - in the right arm
C     IFL=0  - returns the number of hits (amplit>ITHRESH)
C        =1  - returns the sum of amplitudes
C
C
      IMPLICIT NONE
      INTEGER  IARM,IFL,ITHRESH
      INCLUDE ?
C
      INTEGER idetlg,i1,i2,iddd,idg1,idg2,i,k
C
C     -----------------------------------------------------------------
C
      APNUM=0.
      idetlg=1
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
            IF(IARM.EQ.0.OR.IARM.EQ.idg1) THEN
               IF(idg2.GE.ITHRESH) THEN
                  IF(IFL.NE.0) THEN
                     APNUM=APNUM+REAL(idg2)
                  ELSE
                     IF(idg2.GE.ITHRESH) THEN
                        APNUM=APNUM+1.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C     
      END DO

C
 999  RETURN
C
      END

