* $Header:$
* $Log:$
*
      SUBROUTINE OUNTRES(NTFLA)
C
C     ******************************************************************
C     *   Used for e- scattering NTUPLE:                               *
C     *                                                                *
C     *   Users routine                                                *
C     *   Fill the user's part of the NTUPLE                           *
C     *                                                                *
C     *   INPUT :  NTFLA - flag (optional)                             *
C     *                                                                *
C     *   Usage:                                                       *
C     *   1) fill the appropriate COMMON (see ountini.F)               *
C     *                                                                *
C     *    ==>Called by : OMNTUPEV                                     *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER  NTFLA
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
#include "moller/cntescat.inc"
C
C     ------------------------------------------------------------------
C
      NEWOR=0
C
 999  RETURN
      END

