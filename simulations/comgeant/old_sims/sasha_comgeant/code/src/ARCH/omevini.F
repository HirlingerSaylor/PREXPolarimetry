* $Header: /afs/cern.ch/exp/compass/src/cvs/comgeant/code/src/omgbatch/ompro/omevini.F,v 1.1.1.1 1997/06/02 17:39:53 fdr Exp $
* $Log: omevini.F,v $
* Revision 1.1.1.1  1997/06/02 17:39:53  fdr
* Comgeant Monte-Carlo
*
* Revision 3.4  1997/05/05 17:11:04  gen
*  Chudakov: define the include name is removed from all includes
*
* Revision 3.3  1997/05/02 15:20:41  gen
* Chudakov: pipeup, halo, NTUPLE output, some improvements
*
* Revision 3.2.0.1  1996/11/07 19:23:26  las
* First CVS version.
*
*CMZ :          05/10/96  14.24.29  by  E.Chudakov
*-- Author :    E.Chudakov   26/02/93
C
      SUBROUTINE OMEVINI
C
C     ******************************************************************
C     *                                                                *
C     *       INITIALIZE the event                                     *
C     *                                                                *
C     *    ==>Called by : GUKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
C
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omceven.inc"
#include "omgeant/omctrig.inc"
#include "omgeant/omclg.inc"
#include "omgeant/ranma2.inc"
#include "moller/csnaprnd.inc"
#include "moller/comnfla.inc"
#if defined CERNLIB_UNIX 
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id: omevini.F,v 1.1.1.1 1997/06/02 17:39:53 fdr Exp $'/
#endif
      INTEGER IABORT,IERCOD
      COMMON/PABORT/ IABORT,IERCOD
#endif
      INTEGER i,iseq
C
C     ------------------------------------------------------------------
C
      NTRAJEC=0
      KTMASK=0
      KBMASK=0
      EVWEIT=1.
      IBEAMTRA=0
      PBEAM(1)=0.
      PBEAM(2)=0.
      PBEAM(3)=0.
      IVRTPRIM=0
      ITARGET=0
      NTRPOL=0
      NTRUSF=0
      NTRRES=0
      NPAINCL=0
      NENTXYZ=0
C
C---  Trigger part
C
      DO i=1,NTRIGDET
        NTRIGHIT(i)=0
        ETALLHIT(i)=0
      END DO
C
      IHITLG=0
C
C---  NTUPLE data reset
C
      CALL OMNTUPRS
C
#if defined CERNLIB_UNIX
C
C---    Check the abort flag
C
      IF(IABORT.GE.2) THEN
        WRITE(LOUT,FMT='(//'' *** ABORT is called , IABORT='',I2)')
     +       IABORT
        IEORUN=1
        GO TO 999
      ENDIF
C
C---   Print the current status
C
      IF(IABORT.GE.1) THEN
        WRITE(LOUT,FMT='(//'' *** Number of simulated events '',I9)')
     +       IEVENT
#if ! defined CERNLIB_HPUX
        CALL FLUSH(LOUT)
#endif
        IABORT=0
      ENDIF
#endif
C
C---    Prepare a snapshot of random number seeds at the beginning of this event
C
      IF(ISAVRND.NE.0) THEN
         IEVSEED=IEVENT
         DO iseq=1,MXRMSEQ
            DO i=1,103
               ISEEDEV(i,iseq)=0
            END DO
            CALL RMMAQ(ISEEDEV(1,iseq),iseq,'RV')
         END DO
      ENDIF
C
  999 CONTINUE
      RETURN
      END
