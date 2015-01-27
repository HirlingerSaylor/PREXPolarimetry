* $Header:$
* $Log:$
*
      SUBROUTINE ODEGNTIN
C
C     ******************************************************************
C     *                                                                *
C     *   Initializes Degtyarenko's RW NTUPLE:                         *
C     *                                                                *
C     *    ==>Called by : OUINIT                                       *
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
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "moller/cdegtnt.inc"
C
      INTEGER lrec,ierr,i,nentr
      CHARACTER ctit*128
      REAL      rmin(MXDBUF),rmax(MXDBUF)
C
C     ------------------------------------------------------------------
C
      lrec=0
      ierr=0
      IF(LUNDNT.GT.30.AND.LUNDNT.LT.99) THEN
         CALL HROPEN(LUNDNT,'DEGNT','degt.nt',' ',lrec,ierr)
         IF(lrec.EQ.0.OR.ierr.NE.0) THEN
            WRITE(LOUT,*) ' *** ODEGNTIN error: failed to open NTUPLE'
     +               ,' degt.nt , lrec,irec=',lrec,ierr
            ierr=2
         ENDIF
      ELSE
         WRITE(LOUT,*) ' *** ODEGNTIN error: LUNTDEGT=',LUNDNT
     +               ,' is out of range of 31-99'
     +               ,' job terminated'
         ierr=1
      ENDIF
C
      IF(ierr.NE.0) THEN
         IEORUN=1
         GO TO 999
      ENDIF
C
      IDDNT=20
      WRITE(LOUT,1000)
 1000 FORMAT(/10X,'Input NTUPLE of Pavel Degtyarenko format is opened')
      CALL HLDIR(' ','A')
C
      IERR=0
      CALL HRIN(IDDNT,ierr,0)
C
      NVARDNT=MXDBUF
      CALL HGIVEN(IDDNT,ctit,NVARDNT,CNAMDNT,rmin,rmax)
      DO i=1,NVARDNT
         WRITE(LOUT,2000) i,CNAMDNT(i),rmin(i),rmax(i)
 2000    FORMAT(I3,2X,A8,2X,2E15.3)
      ENDDO
C
      IF(NVARDNT.LT.7) THEN
         WRITE(LOUT,*) ' *** ODEGNTIN error: ID=',IDDNT
     +               ,' has too few variables =',NVARDNT
         IEORUN=1
         GO TO 999
      ELSE IF(NVARDNT.GT.7) THEN
         WRITE(LOUT,*) ' === ODEGNTIN message: ID=',IDDNT
     +               ,' has more than 7 variables =',NVARDNT
      ENDIF
C
      CALL HNOENT(IDDNT,nentr)
      WRITE(LOUT,*) '  ID=',IDDNT,'  has ',nentr,' entries'
C
 999  RETURN
      END
