* $Header:$
* $Log:$
*
      SUBROUTINE OMMTOD(XM,XD,IFLAG,LEV)
C.
C.    ******************************************************************
C     *                                                                *
C.    *  Copy of GMTOD (from GEANT) with an extra parameter LEV - level*
C.    *       Computes coordinates XD (in DRS)                         *
C.    *       from known coordinates XM in MRS                         *
C.    *       The local reference system can be initialized by         *
C.    *         - the tracking routines and GMTOD used in GUSTEP       *
C.    *         - a call to GMEDIA(XM,NUMED)                           *
C.    *         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER)            *
C.    *             (inverse routine is GDTOM)                         *
C.    *                                                                *
C.    *        If IFLAG=1  convert coordinates                         *
C.    *           IFLAG=2  convert direction cosinus                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun ,M Hansroul  *********                  *
C.    *                                                                *
C.    ******************************************************************
C
      IMPLICIT NONE
      REAL XD(3),XM(3)
      INTEGER IFLAG,LEV
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
#include "geant321/gcvolu.inc"
C
      INTEGER n
      REAL t1,t2,t3
C.
C.    ------------------------------------------------------------------
C.
      IF(LEV.LT.1.OR.LEV.GT.NLEVEL) THEN
         LEV=0
         GO TO 999
      ENDIF
C
      n=LEV
      IF(IFLAG.EQ.1)THEN
         IF(GRMAT(10,N).NE.0.)THEN
            T1 = XM(1) - GTRAN(1,N)
            T2 = XM(2) - GTRAN(2,N)
            T3 = XM(3) - GTRAN(3,N)
            XD(1) = GRMAT(1,N)*T1+GRMAT(2,N)*T2+GRMAT(3,N)*T3
            XD(2) = GRMAT(4,N)*T1+GRMAT(5,N)*T2+GRMAT(6,N)*T3
            XD(3) = GRMAT(7,N)*T1+GRMAT(8,N)*T2+GRMAT(9,N)*T3
*
         ELSE
            XD(1) = XM(1) - GTRAN(1,N)
            XD(2) = XM(2) - GTRAN(2,N)
            XD(3) = XM(3) - GTRAN(3,N)
*
         ENDIF
*
      ELSE
         IF(GRMAT(10,N).NE.0.)THEN
            XD(1) = GRMAT(1,N)*XM(1)+GRMAT(2,N)*XM(2)+GRMAT(3,N)*XM(3)
            XD(2) = GRMAT(4,N)*XM(1)+GRMAT(5,N)*XM(2)+GRMAT(6,N)*XM(3)
            XD(3) = GRMAT(7,N)*XM(1)+GRMAT(8,N)*XM(2)+GRMAT(9,N)*XM(3)
*
         ELSE
            XD(1) = XM(1)
            XD(2) = XM(2)
            XD(3) = XM(3)
*
         ENDIF
*
      ENDIF
 999  RETURN
      END
