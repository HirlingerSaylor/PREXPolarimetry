* $Header: /afs/cern.ch/exp/compass/src/cvs/comgeant/code/src/omgbatch/ompro/gukine.F,v 1.1.1.1 1997/06/02 17:39:52 fdr Exp $
* $Log: gukine.F,v $
* Revision 1.1.1.1  1997/06/02 17:39:52  fdr
* Comgeant Monte-Carlo
*
* Revision 3.5  1997/05/11 19:45:04  gen
* *** empty log message ***
*
* Revision 3.4  1997/05/02 15:20:29  gen
* Chudakov: pipeup, halo, NTUPLE output, some improvements
*
* Revision 3.3  1997/03/16 22:11:29  gen
* E.Chudakov: modifications #1
*
* Revision 3.2.0.1  1996/11/07 19:23:24  las
* First CVS version.
*
*CMZ :          23/01/96  01.11.18  by  E.Chudakov
*-- Author :    E.Chudakov   01/02/93
C
      SUBROUTINE GUKINE
C
C     ******************************************************************
C     *                                                                *
C     *       Read or Generates Kinematics for primary tracks          *
C     *                                                                *
C     *    ==>Called by : GTRIG                                        *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
#include "geant321/gclist.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omckin.inc"
#include "omgeant/omcgen.inc"
#include "omgeant/omcbeam.inc"
#include "omgeant/omchalo.inc"
#include "omgeant/omceven.inc"
#include "omgeant/omctmp.inc"
#include "omgeant/omcran.inc"
#include "omgeant/omcrun.inc"
#include "moller/cdegtnt.inc"
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id: gukine.F,v 1.1.1.1 1997/06/02 17:39:52 fdr Exp $'/
#endif
C
      REAL     RNDM
      EXTERNAL RNDM
C
      INTEGER ityp,j,it,nver,nusr,ifla(1000),npile,ierr
      REAL    pilem,tpile
      EQUIVALENCE (ifla(1),TMP(201))
C
C     ------------------------------------------------------------------
C
C---    Event initialization
C
      NCOUNTR(1)=NCOUNTR(1)+1
C
      CALL OMEVINI
      IF(IEORUN.NE.0) GO TO 999
C
C---    Random sequence for the beam/inclusive processes
C
      ISEQRAN=1
C
      IF(IDEBUG.NE.0) THEN
        WRITE(LOUT,2000) IEVENT
      ENDIF
C
C---     Read the kinematic data?
C
      IF(KINGET.GT.0) THEN
        CALL OMKINGET
      ENDIF
C
C---     Read the Degtyarenko's NTUPLE
C
      IF(LUNDNT.NE.0) THEN
        CALL ODEGET
      ENDIF
      IF(IEOTRI.NE.0) GO TO 999
C
C---    Normal kinematics
C
      IF(KINTST(1).EQ.0) THEN
C
C---    Simulate the beam and the interaction point
C
        CALL OMBEAM(0.)
C
        IF(IEOTRI.NE.0) GO TO 999
        NCOUNTR(2)=NCOUNTR(2)+1
C
C---      Inclusive production
C
        IF(NINCL.GT.0) THEN
          IF(KINCHARM.NE.0.AND.NINCL.EQ.2) THEN
            CALL OKINCHAR
          ELSE
            CALL OMINCL
          ENDIF
        ENDIF
C
        IF(IEOTRI.NE.0) GO TO 999
        NCOUNTR(3)=NCOUNTR(3)+1
C
C---       Multiplicity generator (FRITIOF...)
C
        ISEQRAN=2
        CALL OMKINE
C
      ELSE
C
C---      Special kinematics: test of material using one track
C
        CALL GSVERT(PKINTST(4),0,0,0,0,IVRTPRIM)
        CALL GSKINE(PKINTST(7),KINTST(3),IVRTPRIM,0,0,it)
        IBEAMTRA=0
C
      ENDIF
C
C---    Restore the random sequence for GEANT
C
      ISEQRAN=3
C
      IF(IEOTRI.NE.0) GO TO 999
      NCOUNTR(4)=NCOUNTR(4)+1
C
C---       Check the particle types
C
      IF(NTYPAR.GT.0) THEN
C
        DO it=1,NTRACK
          ifla(it)=0
        END DO
C
        DO j=1,NTYPAR
          IF(KTYPAR(j).NE.0) THEN
            DO it=1,NTRACK
              IF(ifla(it).EQ.0) THEN
                CALL GFKINE(it,TMP(1),TMP(4),ityp,nver,TMP(10),nusr)
                IF(ityp.EQ.KTYPAR(j)) THEN
                  ifla(it)=1
                  NTRUSF=NTRUSF+1
                  ITRUSF(NTRUSF)=it
                  GO TO 40
                ENDIF
              ENDIF
            END DO
C
            IEOTRI=1
            GO TO 50
C
          ENDIF
C
   40     CONTINUE
C
        END DO
C
   50   CONTINUE
C
      ENDIF
C
C---   Add the beam pileup?
C
      pilem=PILEUP(1)*2.*PILEUP(2)*1.E-9
      IF(pilem.GT.0.001) THEN
        CALL POISSN(pilem,npile,ierr)
        IF(ierr.EQ.0) THEN
          DO j=1,npile
            tpile=PILEUP(2)*(-1.+RNDM(j)*2.)
            CALL OMBEAM(tpile)
          END DO
        ELSE
          WRITE(LOUT,2100) pilem,ierr
        ENDIF
      ENDIF
C
C---   Add the halo pileup?
C
      IF(IHALOPAR.NE.0) CALL OMHALO
      tpile=HALOORIG(1,1)
C
      IF(IDEBUG.NE.0.AND.ISWIT(2).NE.0) THEN
        CALL GPVERT(0)
        CALL GPKINE(0)
      ENDIF
C
  999 RETURN
 2000 FORMAT(//1X,20('='),5X,' Start the simulation of the event ',I6
     +      ,5X,20('='))
 2100 FORMAT(' *** GUKINE error: PILEUP POISSN mean=',E12.3,' err=',I5)
      END

