* $Header:$
* $Log:$
*
#ifndef CERNLIB_OMCUNIT_INC
* #define CERNLIB_OMCUNIT_INC
C
C---           Degtyarenko's NTUPLE
C
      INTEGER    MXDBUF
      PARAMETER	(MXDBUF=10)
      COMMON/CDEGTNT/  LUNDNT,IDDNT,NVARDNT
     +                ,BUFDNT(MXDBUF)
C
      INTEGER          LUNDNT,IDDNT,NVARDNT
      REAL             BUFDNT
C
      COMMON/CDEGTNTC/ CNAMDNT(MXDBUF)
      CHARACTER        CNAMDNT*8
C
C
#if defined CERNLIB_COMMENT
C
C=F=   LUNDNT          = LUN for Degtyarenko's input NTUPLE, =0 - do not open/read the NTUPLE
C
C      IDDNT           - ID (=20)
C      NVARDNT         - dimension
C      BUFDNT(1)       - GEANT particle type +100 or 200
C            (2)       - kinetic energy
C            (3)       - theta angle
C            (4)       - phi angle
C            (5-7)     - the X,Y,Z coordinates of the origin
#endif
C
#endif






