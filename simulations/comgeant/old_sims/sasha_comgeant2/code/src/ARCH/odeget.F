* $Header:$
* $Log:$
*
      SUBROUTINE ODEGET
C
C     ******************************************************************
C     *                                                                *
C     *       Read P.Degtyarenko's NTUPLE                              *
C     *                                                                *
C     *    ==>Called by : GUKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "omgeant/omcfla.inc"
#include "moller/cdegtnt.inc"
C
      INTEGER   i,nbuf,itpar,nevnt,ierr,itypp,nvtx,ntrk
      CHARACTER cnampar*20,cntdir*5
      REAL      ampar,chpar,tlpar,ubuf(20),pp(3),ptot
C
C     ------------------------------------------------------------------
C
C      WRITE(LOUT,*) 'Get event ...'
      cntdir='DEGNT'
      CALL HCDIR('//'//cntdir,' ')
C
      CALL HNOENT(IDDNT,nevnt)
      IF(nevnt.LT.IEVENT) THEN
         WRITE(LOUT,1000) nevnt,IEVENT
 1000    FORMAT(/2X,'=== ODEGET: no events left in the input NTUPLE '
     +         ,' NEV=',I8,' at the GEANT event=',I8)
         IEORUN=1
         GO TO 999
      ENDIF
C
      CALL HGNPAR(IDDNT,'ODEGET')
      CALL HGNF(IDDNT,IEVENT,BUFDNT,ierr)
      IF(IDEBUG.NE.0) THEN
         WRITE(LOUT,FMT='(A12,I8,7E12.4)') 'Event :',IEVENT
     +                  ,(BUFDNT(i),i=1,7)
      ENDIF
      IF(ierr.NE.0) THEN
         WRITE(LOUT,1100) ierr,IEVENT
 1100    FORMAT(/2X,'*** ODEGET: error in the input NTUPLE '
     +         ,' ierr=',I8,' at the GEANT event=',I8)
         IEORUN=1
         GO TO 999
      ENDIF
C
      itypp=MOD(INT(BUFDNT(1)),100)
      ampar=-100.
      CALL GFPART(itypp,cnampar,itpar,ampar,chpar,tlpar
     +              ,ubuf,nbuf) 
C
      IF(ampar.LT.0..OR.ampar.GT.100.) THEN
         WRITE(LOUT,1200) itypp,ampar,IEVENT
 1200    FORMAT(/2X,'*** ODEGET: input NTUPLE '
     +      ,'unknown particle type ',I8,' at the GEANT event=',I8)
         IEOTRI=1
         GO TO 999
      ENDIF
C
      ptot=SQRT((ampar+BUFDNT(2))**2-ampar**2)
      pp(1)=ptot*SIN(BUFDNT(3))*COS(BUFDNT(4))
      pp(2)=ptot*SIN(BUFDNT(3))*SIN(BUFDNT(4))
      pp(3)=ptot*COS(BUFDNT(3))
C
C---     Store the particle
C
      CALL GSVERT(BUFDNT(5),0,0,0,0,nvtx)
      CALL GSKINE(pp(1),itypp,nvtx,0,0,ntrk)
C
 999  RETURN
      END
