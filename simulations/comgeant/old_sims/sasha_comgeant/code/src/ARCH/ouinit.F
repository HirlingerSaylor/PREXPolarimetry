* $Header:$
* $Log:$
*
      SUBROUTINE OUINIT
C
C     ******************************************************************
C     *                                                                *
C     *       USER's routine called at the end of initialization       *
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
#include "geant321/gcflag.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omcbeam.inc"
#include "omgeant/omcunit.inc"
#include "omgeant/ranma2.inc"
#include "omgeant/cescatt.inc"
#include "moller/cdegtnt.inc"
#include "moller/comnfla.inc"
C
      COMMON/CCSMOTT/ CSDMOTT(4),ALMOTT
      REAL            CSDMOTT,ALMOTT
C---                  CSMOTT(1:3) = E0 electron (LAB), Z target, A target.
C---                  ALMOTT      - exponential slope of the variable
C
      COMMON/CPHOTPR/ EBEAMPH   ! Bremstr. photoproduction: beam energy
     +               ,IPARTPH   ! Produced particle 
     +               ,TARGLPH   ! H2 target length (cm) 
     +               ,BFLUXPH   ! e- beam flux (particles per time_interval/pulse...
      REAL            EBEAMPH,TARGLPH,BFLUXPH
      INTEGER         IPARTPH
C
C
      COMMON/CKINTMP/ KEMPROC,EEMPROC,THCMEMPROC
      INTEGER KEMPROC
      REAL EEMPROC,THCMEMPROC
C
      INTEGER  IDHNEW
      REAL     HSUM,CSMOTT,PHOTSPEC,RADMOTT
      EXTERNAL HSUM,CSMOTT,PHOTSPEC,RADMOTT
C,ANPWEMPR
C      REAL HRNDM1,RADMOT1
C      EXTERNAL    RADMOT1
C
      INTEGER i,itart,nbuf,itrttg,itrtut,nbin,nbinl,lrec,mag
      CHARACTER namtarg*20,cproc*16,cn*4
      REAL    r0,atarg,chtarg,tltarg,ubuf(20),csfac,gam,ebeam,crosec
     +       ,ctstep,phmx,x1,x2,ratem,rateh,t1,t2,y1,y2,crosec1
     +       ,xl1,xl2,aa,bb,ratefac
C
      DATA r0/2.818/  ! electron radius in fm
C
C     ------------------------------------------------------------------
C
      ratefac=6.02*1.E-1*1.E-6/1.602E-19 ! Avog, 1uA
C
C      WRITE(6,*) ' IESCATT=',IESCATT
C
      IF(IOMEXKIN.EQ.1.OR.IOMEXKIN.EQ.6) THEN
C
C---    Moller/Bhabha/Annih scattering
C
         CALL OMOLINI
C
      ELSE IF(IOMEXKIN.EQ.2) THEN
C
C---    Mott     
C
         CALL OMOTTINI
C
      ELSE IF(IOMEXKIN.EQ.3) THEN
C
         cproc='Bremst. photopr.'
         IDPHOTPR=IDHNEW(-1)
         EBEAMPH=(BEAMOLIM(1)+BEAMOLIM(2))/2.
         IPARTPH=IPHOTPR
         BFLUXPH=1.E11
         TARGLPH=150.
         CALL HBFUN2(IDPHOTPR,'Photo-prod theta/pmom'
     +       ,100,THETPHPR(1),THETPHPR(2)
     +       ,100,PMOMPHPR(1),PMOMPHPR(2),PHOTSPEC)
C         
         rateh=HSUM(IDPHOTPR)*(THETPHPR(2)-THETPHPR(1))/100
     +                       *(PMOMPHPR(2)-PMOMPHPR(1))/100
         WRITE(LOUT,1300) cproc,IPARTPH,rateh
     +                         ,THETPHPR(1),THETPHPR(2)
     +                         ,PMOMPHPR(1),PMOMPHPR(2)
     +                         ,TARGLPH,BFLUXPH
C
      ELSE IF(IOMEXKIN.EQ.4) THEN
C
         cproc='VMeson photopr.'
C
         IF(IPHOTVMP.LE.0.OR.IPHOTVMP.GT.99999.OR.
     +      IPHOTVMT.LE.0.OR.IPHOTVMT.GT.99999.OR.
     +      PHOTVMTT.LE.0..OR.PHOTVMTT.GT.9999.) THEN
            WRITE(LOUT,1500) cproc,IPHOTVMP,IPHOTVMT,PHOTVMTT
            IEORUN=1
         ENDIF
C
      ELSE IF(IOMEXKIN.EQ.5) THEN
C
         cproc='Mott on lead'
C
C!         CALL NSETUP
C
         WRITE(LOUT,1100) cproc,THETMOT(1),THETMOT(2)
     +                   ,PHIMOT(1),PHIMOT(2),0.,0.
C
C
      ENDIF
C
      IF(LUNDNT.NE.0) CALL ODEGNTIN
C
C---      File with event RNDM seeds
C
      IF(ISAVRND.NE.0) THEN
         lrec=(MXRMSEQ*103+1)*4
         OPEN(LUNOMG(5),FILE='cntevseeds.dat',ACCESS='DIRECT'
     +       ,FORM='UNFORMATTED',RECL=lrec)
         WRITE(LOUT,1700) LUNOMG(5),lrec
      ENDIF
C
 999  RETURN
 1050 FORMAT( ' === ',A16,' cross-section '
     +     ,' x=exp(-',E12.4,'*theta), x1,x2=',4E12.4)
 1100 FORMAT(/' === ',A16,' cross-section '
     +     ,F11.6,'<theta LAB<',F11.6,2X,F7.2,'<phi LAB<',F7.2,' deg'
     +     ,' are ',E12.4,2X,E12.4,' barn per one target nucleus')
 1200 FORMAT(' === ',A16,' event rate/sec =',E12.4
     +     ,' for 0.003cm Fe target, 1.microA'
     +     /,' === ',A16,' event rate/sec =',E12.4
     +     ,' for 100.0cm H2 target, 1.microA')
 1300 FORMAT(/' === ',A16,' particle',I2,' rate=',E12.3,'  for '
     +     ,F11.6,'<theta LAB<',F11.6,2X,F7.3,'<P/Pbeam LAB<',F7.3
     +     ,' Liquid H2 target ',F5.1,'cm',3X,E10.3,' beam particles')
 1500 FORMAT(/' *** OUINIT ',A16,' Error: error in IPHOTVMP'
     +     ,',IPHOTVMT,PHOTVMTT',3X,3I11)
 1600 FORMAT(/' *** OUINIT Error: Annihilation: beam type .NE. 2'
     +     ,' IBEAMPAR=',I3,'  Exit!')
 1700 FORMAT(//' === FILE OPENED cntevseeds.dat on lun=',I3
     +  ,'  for random generator seeds for events written to NTUPLE'
     +  ,' lrec=',I6)
 1800 FORMAT(/4X,'magnet type mother     Xcen      Ycen      Zcen     '
     +      ,'rot  '
     +      ,' 1/2 size: X       Y       Z   field/scale     flags ')
 1820 FORMAT(6X,I2,3X,I3,3X,A4,2X,3F10.2,3X,I3,8X,3F8.2,2X,F10.4
     +      ,2X,6F7.1)
      END
