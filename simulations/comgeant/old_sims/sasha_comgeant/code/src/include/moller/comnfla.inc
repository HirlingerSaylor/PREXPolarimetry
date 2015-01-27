* $Header:$
* $Log:$
*
#ifndef CERNLIB_OMCUNIT_INC
* #define CERNLIB_OMCUNIT_INC
C
C---           Various flags
C
      COMMON/COMNFLA/ CHERFAC,ISAVRND
C
      REAL            CHERFAC
      INTEGER         ISAVRND
C
#if defined CERNLIB_COMMENT
C
C===   CHERFAC   - a factor for the number of photons per 1cm, (used in GNCKOV) Default=369.81E9
C=F=   ISAVRND   - >0  - write down a direct access file with the random number generator (RMMAQ) 
C                        seeds for each event, written to NTUPLEs 
#endif
C
#endif
