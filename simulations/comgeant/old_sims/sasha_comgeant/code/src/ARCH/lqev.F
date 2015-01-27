
C **********************************************************************

      SUBROUTINE LQEV

C...Generate an ordinary 2-jet event, q-event.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
#include "jetset74/lujets.inc"                                          
CCC      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
#include "jetset74/ludat1.inc"
CCC      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
#include "jetset74/ludat2.inc"
CCC      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C
      INTEGER itry

      itry=0

      LST(24)=1
      W=SQRT(W2)

C...Choose flavour of scattered quark and target remnant.
  200 CALL LFLAV(IFL,IFLR)
      itry=itry+1
      IF(itry.GT.1000) THEN
         WRITE(6,2000) itry,W,AMIFL+AMIFLR+0.9+PARJ(32)
 2000    FORMAT(' *** LQEV (LEPTO): infinite loop broken at ',I5
     +   ,' try. May be the W**2 limit is too low. '
     +   ,' Current W=',F8.3,' , required at least ',F8.3)
         GO TO 998
      ENDIF
C
C      WRITE(6,*) 'LQEV 1',IFL,IFLR,LST(21)
      IF(LST(21).NE.0) GOTO 200

      GOTO 210
C...Entry used for Ariadne
      ENTRY LQEVAR(IFLAR,IFLRAR)
      IFL=IFLAR
      IFLR=IFLRAR
      LST(24)=1
      W=SQRT(W2)

  210 CONTINUE
      MSTJ(93)=1
C      WRITE(6,*) 'LQEV 2',IFL,IFLR
      AMIFL=ULMASS(IFL)
      MSTJ(93)=1
      AMIFLR=ULMASS(IFLR)
      IF(LST(14).EQ.0.OR.IFLR.GT.10
     &.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
C...Check if energy in jet system is enough for fragmentation.
C      WRITE(6,*) 'LQEV 3',W,AMIFL+AMIFLR+PARJ(32),AMIFL,AMIFLR,PARJ(32)
        IF(W.LT.AMIFL+AMIFLR+PARJ(32)) GOTO 200
        CALL LU2ENT(MSTU(1),IFL,IFLR,W)
        K(MSTU(1)+1,3)=2
      ELSE
C...Target remnant is not a simple diquark, special treatment needed.
C      WRITE(6,*) 'LQEV 4',W,AMIFL+AMIFLR+0.9+PARJ(32)
        IF(W.LT.AMIFL+AMIFLR+0.9+PARJ(32)) GOTO 200
        IFLRO=IFLR
        NREMH=0
  300   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 999
C...Give balancing pt to IFLQ and IFLQQ.
        CALL LPRIKT(PARL(14),PT,PHI)
        CALL LREMH(IFLRO,PT,IFLR,K2,XT)
        MSTJ(93)=1
        AMIFLR=ULMASS(IFLR)
        MSTJ(93)=1
        AMK2=ULMASS(K2)
        PT2=PT**2
        TM2K2=AMK2**2+PT2
        EK2=.5*(XT*W+TM2K2/XT/W)
        PZK2=-.5*(XT*W-TM2K2/XT/W)
        EPZ=W-TM2K2/XT/W
        WT=(1.-XT)*W*EPZ-PT2
C...Check if energy in jet system is enough for fragmentation.
        IF(WT.LT.(AMIFL+AMIFLR+PARJ(32))**2) GOTO 300
        WT=SQRT(WT+PT2)
        TMIFLR=AMIFLR**2+PT2
        EIFL=.5*(WT+(AMIFL**2-TMIFLR)/WT)
        EIFLR=.5*(WT+(-AMIFL**2+TMIFLR)/WT)
        THER=ULANGL(-SQRT(EIFLR**2-TMIFLR),PT)
C...Form jet system.
        CALL LU1ENT(-MSTU(1),IFL,EIFL,0.,0.)
        CALL LU1ENT(MSTU(1)+1,IFLR,EIFLR,THER,PHI)
        CALL LUDBRB(MSTU(1),0,0.,0.,0.D0,0.D0,
     &  (DBLE(EPZ)-(1.D0-DBLE(XT))*DBLE(W))/
     &  (DBLE(EPZ)+(1.D0-DBLE(XT))*DBLE(W)))
        THEK2=ULANGL(PZK2,PT)
C...Add formed "target" particle.
        MSTU(10)=1
        P(MSTU(1)+2,5)=AMK2
        CALL LU1ENT(MSTU(1)+2,K2,EK2,THEK2,PHI+3.1415927)
        MSTU(10)=2
        K(MSTU(1)+1,3)=2
        K(MSTU(1)+2,3)=2
CIC...Target remnants required to go backwards in hadronic cms
        IF(P(MSTU(1)+1,3).GT.0..OR.P(MSTU(1)+2,3).GT.0.) GOTO 300
      ENDIF
C
 998  CONTINUE
C
      LST(21)=0
      RETURN

  999 LST(21)=3
      RETURN
      END
