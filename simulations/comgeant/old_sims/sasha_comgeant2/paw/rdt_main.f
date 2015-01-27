      PROGRAM RDT_MAIN
C
C---    Test the OMSOLMAP
C
      INTEGER mag
      REAL BZR(2)
C
      mag=2
 10   WRITE(6,*) 'Enter Z,R'
C
      READ(5,*) z,r
      CALL OMSOLMAP(mag,z,r,bzr)
      WRITE(6,*) ' bzr ',bzr
      GO TO 10
      END
C
      INCLUDE 'omsolmap.f'
