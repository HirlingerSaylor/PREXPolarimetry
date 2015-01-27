      REAL FUNCTION RAN_TST(T)
C     
      CALL RANNOR(a,b)
      WRITE(6,*) a,b
      RAN_TST=a
      END
