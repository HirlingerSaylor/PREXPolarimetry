* $Header:$
* $Log:$
*
      SUBROUTINE OATOMIC(SHELL,E00,CTCM,STCM,CPP,SPP,PLAB)
C
C========  Scattering on a bound electron
C========  Author: L. Levchuk
C
      IMPLICIT REAL*8(A-H,O-Z)

#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      PARAMETER (PI=3.14159D0)
      CHARACTER*1 SHELL
      REAL*8 ELFE(7)/524.0D 00,62.3D 00,53.4D 00,52.4D 00,7.1D 00,
     *4.2D 00,0.D 00/
      INTEGER NELFE(7)/2,2,2,4,2,6,8/
      REAL*8 ELCO(7)/568.2D 00,68.6D 00,59.0D 00,57.8D 00,7.7D 00,
     *4.7D 00,0.D 00/
      INTEGER NELCO(7)/2,2,2,4,2,6,9/
      DIMENSION E(7),PT(3)
      REAL*4 E00,CTCM,STCM,CPP,SPP,PLAB(4,2),RNDM
C
      IF(SHELL.NE.'K'.AND.SHELL.NE.'L') RETURN 
C
      W=0.511D-3
      E0=E00
      P0MOD=DSQRT(E0**2-W**2)
      SP=SPP
      CP=CPP
      CTC=CTCM
      STC=STCM
      ACO=58.D0
      ZCO=27.D0
      ELPOLC=1.720D0
      AFE=56.D0
      ZFE=26.D0
      ELPOLF=2.216D0
      A=0.5D0*(AFE+ACO)
      Z=0.5D0*(ZFE+ZCO)
      ELPOL=(ELPOLF+ELPOLC)/(ZFE+ZCO)
      AT=.9385D0*A
      AT=AT-.007D0*A
      R=13.6D-9
      DO 1 J=1,4
      E(J)=R*0.5D0*(ELFE(J)+ELCO(J))
    1 CONTINUE
      IF(SHELL.EQ.'K') EN=E(1)
      IF(SHELL.EQ.'L') EN=.25D0*E(2)+.25D0*E(3)+.50D0*E(4)
      E1 = .5D0*(E0+W) + .5D0*(E0-W)*CTC
***************************************************************
C     WHAT FOLLOWS IS A SIMPLIFIED BUT REASONABLE CONSIDERATION
***************************************************************
      ANGMOL=DSQRT(2.D0*W/(E0+W))*STC/(1.D0+CTC)
      CTATOM = 1.D0 - 2.D0*RNDM(E00)
      PTMOD=DSQRT(2.D0*W*EN)
      ARG = 1.D0 - EN/W + (PTMOD/W)*CTATOM
***      PRINT 1000,ARG
*** 1000 FORMAT(10X,'ARG=',E12.4)
      TANANG = DSQRT(ARG)*ANGMOL
      CT=1.D0/(1.D0+TANANG**2)
      ST=TANANG/(1.D0+TANANG**2)
      PHIT=2.*PI*RNDM(SPP)
      STATOM=SQRT(1.D0 - CTATOM**2)
      PT(1)=PTMOD*STATOM*COS(PHIT)
      PT(2)=PTMOD*STATOM*SIN(PHIT)
      PT(3)=PTMOD*CTATOM
***************************************************************
      P1MOD=SQRT(E1**2-W**2)
      PLAB(1,1)=P1MOD*ST*CPP
      PLAB(2,1)=P1MOD*ST*SPP
      PLAB(3,1)=P1MOD*CT
      PLAB(4,1)=E1
      ER = E0 + W - EN - E1 - .5D0*PTMOD**2/AT
      P2MOD=SQRT(ER**2-W**2)
      PLAB(1,2) = -PLAB(1,1) -PT(1)
      PLAB(2,2) = -PLAB(2,1) -PT(2)
      PLAB(3,2) =  P0MOD -PLAB(3,1) -PT(3)
      PLAB(4,2) = ER
      END
