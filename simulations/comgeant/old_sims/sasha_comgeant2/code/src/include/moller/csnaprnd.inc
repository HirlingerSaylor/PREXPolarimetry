* $Header:$
* $Log:$
*
#ifndef CERNLIB_CSNAPRND_INC
* #define CERNLIB_CSNAPRND_INC
C
C---       Snapshot of all random number generators' seeds, stored at 
C                             the beginning of each event 
C=== Important!! :          ranma2.inc should be included before this include
C
      COMMON/CSNAPRND/ IEVSEED,ISEEDEV(103,MXRMSEQ)
      INTEGER          IEVSEED,ISEEDEV
C
#if defined CERNLIB_COMMENT
C
C===   IEVSEED   - event number
C=F=   ISEEDEV   - values
C
#endif
C
#endif
