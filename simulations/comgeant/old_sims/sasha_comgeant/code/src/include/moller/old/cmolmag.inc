* $Header:$
* $Log:$
*
#ifndef CERNLIB_OMCUNIT_INC
* #define CERNLIB_OMCUNIT_INC
C
C---           Moller magnets
C
      INTEGER    MXMOLMAG,MXMOLMW
      PARAMETER (MXMOLMAG=9,MXMOLMW=10)
      COMMON/CMOLMAG/ NMOLMAG,MOLMAGP(MXMOLMW,MXMOLMAG)
C
      INTEGER         NMOLMAG,MOLMAGP
      REAL                    PMOLMAGP(MXMOLMW,MXMOLMAG)
      EQUIVALENCE (MOLMAGP(1,1),PMOLMAGP(1,1))
C
#if defined CERNLIB_COMMENT
C
C===   NMOLMAG - number of defined elements, molm=1,NMOLMAG
C=F=   MOLMAGP(1,molm) - type =0 - not defined, =1 -quad, =2 - dipole, =3 - solenoid
C     PMOLMAGP(2,molm) - X - size (Q - diameter) of the field   
C             (3,molm) - Y - size (Q - diameter)   
C             (4,molm) - Z - size 
C            .......     not used at the moment
#endif
C
#endif
