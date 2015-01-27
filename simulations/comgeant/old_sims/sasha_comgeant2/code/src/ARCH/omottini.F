* $Header:$
* $Log:$
*
      SUBROUTINE OMOTTINI
C
C     ******************************************************************
C     *                                                                *
C     *       Mott initialization                                      *
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
#include "geant321/gcflag.inc"
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
C #include "omgeant/omcfla.inc"
#include "omgeant/omcbeam.inc"
C #include "omgeant/omcunit.inc"
#include "omgeant/cescatt.inc"
C
      COMMON/CCSMOTT/ CSDMOTT(4),ALMOTT
      REAL            CSDMOTT,ALMOTT
C---                  CSMOTT(1:3) = E0 electron (LAB), Z target, A target.
C---                  ALMOTT      - exponential slope of the variable
C
C
      INTEGER  IDHNEW
      REAL     HSUM,CSMOTT,RADMOTT
      EXTERNAL IDHNEW,HSUM,CSMOTT,RADMOTT
C
      INTEGER nbuf,itrttg,itrtut,nbin,nbinl
      CHARACTER namtarg*20,cproc*16
      REAL    atarg,chtarg,tltarg,ubuf(20),ebeam,crosec
     +       ,phmx,x1,x2,ratem,rateh,t1,t2,y1,y2,crosec1
     +       ,xl1,xl2,aa,bb,ratefac
C
C     ------------------------------------------------------------------
C
      ratefac=6.02*1.E-1*1.E-6/1.602E-19 ! Avog, 1uA
C
      cproc='         Mott   '
C
C---       Find the exp slope 
C

      ebeam=(BEAMOLIM(1)+BEAMOLIM(2))/2.
      CSDMOTT(1)=ebeam
      CALL GFPART(MOTTARG,namtarg,itrttg,atarg,chtarg,tltarg
     +              ,ubuf,nbuf) 
      IF(itrttg.EQ.0) THEN
         WRITE(LOUT,*) '*** OUINIT error: MOTTARG is not defined'
     +                    ,MOTTARG
         atarg=1.
         chtarg=1.
         IEORUN=1
      ENDIF
      CSDMOTT(2)=chtarg
      CSDMOTT(3)=atarg
      ALMOTT=0.
C
      phmx=(PHIMOT(2)-PHIMOT(1))*PI/180.
      crosec1=0.
C
      t1=THETMOT(1)*PI/180.
      t2=THETMOT(2)*PI/180.
      IDMOTT=IDHNEW(-1)
      IF(IRADMOT.EQ.0) THEN
         nbin=1000
         CALL HBFUN1(IDMOTT,'Mott vs theta',nbin,t1,t2,CSMOTT)
         crosec=HSUM(IDMOTT)*(t2-t1)/nbin*phmx
      ELSE
         cproc=' Radiative Mott '
C          nbin=200
         nbin=20
         nbinl=10
         xl1=RADMOTL(1)
         xl2=RADMOTL(2)
         CALL HBFUN2(IDMOTT,'Rad. Mott vs theta vs Eloss'
     +                 ,nbin,t1,t2,nbinl,xl1,xl2,RADMOTT)
         crosec=HSUM(IDMOTT)*(t2-t1)/nbin*phmx*(xl2-xl1)/nbinl*ebeam
C            CALL HBFUN1(102,'Mott '
C     +                 ,nbin,t1,t2,RADMOT1)
      ENDIF

      IF(IMOTTWG.NE.2) THEN
         IF(IRADMOT.EQ.0) THEN
            y1=CSMOTT(t1)
            y2=CSMOTT(t2)
         ELSE
            y1=RADMOTT(t1,0.5)
            y2=RADMOTT(t2,0.5)
         ENDIF
         ALMOTT=-(ALOG(y2)-ALOG(y1))/(t2-t1)
         ALPHMOTT=ALMOTT
         x1=EXP(-ALMOTT*t2)
         x2=EXP(-ALMOTT*t1)
C     
         IF(IRADMOT.EQ.0) THEN
C             CALL HCOPY(IDMOTT,102,'    ')
            CALL HDELET(IDMOTT)
            CALL HBFUN1(IDMOTT,'Mott vs exp(-aph*thet)'
     +           ,nbin,x1,x2,CSMOTT)
            crosec1=HSUM(IDMOTT)*(x2-x1)/nbin*phmx
         ELSE
            IDRADMOT=IDHNEW(-1)
            CALL HBFUN2(IDRADMOT
     +           ,'Rad. Mott vs exp(-alp*thet)theta vs Eloss'
     +           ,nbin,x1,x2,nbinl,xl1,xl2,RADMOTT)
         ENDIF
         WRITE(LOUT,1050) cproc,ALMOTT,x1,x2,y1,y2
      ENDIF          
C     
      ratem=crosec*0.003*7.87/55.85*ratefac
      IF(CSDMOTT(2).LT.1.1) ratem=ratem*26.
      rateh=crosec*100.0*0.07080/1.*ratefac
      WRITE(LOUT,1100) cproc,THETMOT(1),THETMOT(2)
     +                   ,PHIMOT(1),PHIMOT(2),crosec,crosec1
      WRITE(LOUT,1200) cproc,ratem,cproc,rateh
C
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
      END
