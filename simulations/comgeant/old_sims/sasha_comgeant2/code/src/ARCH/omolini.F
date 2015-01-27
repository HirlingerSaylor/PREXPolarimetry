* $Header:$
* $Log:$
*
      SUBROUTINE OMOLINI
C
C     ******************************************************************
C     *                                                                *
C     *       Moller/Bhabha/Annih initialization                       *
C     *                                                                *
C     *    ==>Called by : OMGINIT                                      *
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
C #include "geant321/gcflag.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omcbeam.inc"
C #include "moller/cmolmag.inc"
#include "omgeant/cescatt.inc"
C #include "moller/cdegtnt.inc"
C
C
      COMMON/CKINTMP/ KEMPROC,EEMPROC,THCMEMPROC
      INTEGER KEMPROC
      REAL EEMPROC,THCMEMPROC
C
      INTEGER  IDHNEW
      REAL     HSUM
      REAL     SIGEMPR
      EXTERNAL HSUM
      EXTERNAL SIGEMPR
      LOGICAL  HEXIST
C
      INTEGER nbin
      CHARACTER cproc*16
      REAL    r0,csfac,ebeam,crosec
     +       ,ctstep,phmx,x1,x2,ratem,rateh,ratefac
C
      DATA r0/2.818/  ! electron radius in fm
C
C     ------------------------------------------------------------------
C
      ratefac=6.02*1.E-1*1.E-6/1.602E-19 ! Avog, 1uA
C
      IF(IBEAMPAR.EQ.3) THEN
         cproc='        Moeller '
         KEMPROC=1
      ELSE IF(IBEAMPAR.EQ.2) THEN
         IF(IOMEXKIN.EQ.1) THEN
            cproc='        Bhabha '
            KEMPROC=2
         ELSE IF(IOMEXKIN.EQ.6) THEN
            cproc='  Annihilation '
            KEMPROC=3
         ENDIF
      ELSE
         WRITE(LOUT,1020) IBEAMPAR
         GO TO 999
      ENDIF
C
      IF(IDMOLTET.NE.0) THEN
         IF(HEXIST(IDMOLTET)) CALL HDELET(IDMOLTET)
      ENDIF
C
      IDMOLTET=IDHNEW(-1)
      nbin=1200
      x1=COS(THETMOL(2)*PI/180.)
      x2=COS(THETMOL(1)*PI/180.)
      ebeam=(BEAMOLIM(1)+BEAMOLIM(2))/2.
      EEMPROC=ebeam
      CALL HBFUN1(IDMOLTET,'Moller Theta$',nbin,x1,x2,SIGEMPR)
      csfac=1.
      ctstep=(x2-x1)/nbin
      phmx=(PHIMOL(2)-PHIMOL(1))*PI/180.
      crosec=csfac*HSUM(IDMOLTET)*ctstep*phmx
      ratem=crosec*0.003*7.87/55.85*26.*ratefac
      rateh=crosec*100.0*0.0708/1.*1.*ratefac
      WRITE(LOUT,1000) cproc,THETMOL(1),THETMOL(2)
     +                ,PHIMOL(1),PHIMOL(2),crosec
      WRITE(LOUT,1200) cproc,ratem,cproc,rateh
C
C
 999  RETURN
 1020 FORMAT(/' *** OMOLINI Error: Moller/Bhabha not defined for beam'
     +     ,'  IBEAMPAR=',I3,'  Exit!')
 1000 FORMAT(/' === ',A16,' cross-section '
     +     ,F7.2,'<theta CM<',F7.2,2X,F7.2,'<phi CM<',F7.2,' deg'
     +     ,' is ',E12.4,' barn per one target electron')
 1200 FORMAT(' === ',A16,' event rate/sec =',E12.4
     +     ,' for 0.003cm Fe target, 1.microA'
     +     /,' === ',A16,' event rate/sec =',E12.4
     +     ,' for 100.0cm H2 target, 1.microA')
      END
C
      REAL FUNCTION SIGEMPR(COSTHET)
C
C===      Moller/Bhabha/annih.. COS(Theta) distrubution
C
      IMPLICIT NONE
      REAL COSTHET,ct,thet
      REAL     OSIG_EM
      EXTERNAL OSIG_EM
C
      ct=COSTHET
      IF(ct.LT.-1.) ct=-1.
      IF(ct.GT.1.)  ct=1.
      thet=ACOS(ct)*180./3.14159
      SIGEMPR=OSIG_EM(thet)
C
      RETURN
      END
