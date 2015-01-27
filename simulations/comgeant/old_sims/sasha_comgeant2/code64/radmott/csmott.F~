* $Header:$
* $Log:$
*
      REAL FUNCTION CSMOTT(THETX)
C
C---     Mott cross-section d(sigma)/d(phi)/d(theta) in LAB in barn for:
C---         THETX                 - scattering angle in Lab
C                                 
C
      IMPLICIT NONE
      REAL THETX
C
      COMMON/CCSMOTT/ CSDMOTT(4),ALMOTT
      REAL            CSDMOTT,ALMOTT
C---                  CSMOTT(1:3) = E0 electron (LAB), Z target, A target.
C---                  ALMOTT      - exponential slope of the variable
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
C      VECTOR VMOTT(3)
C
      DOUBLE PRECISION costh,sinth,cos2,sin2,q2,e0,p0,e,p,ep,pp
     +    ,cross,ztar,atar,ecme,pcm,s,sqs,gam,bet,costhcm,cosp,sinp
     +    ,thetp,ge,gm,crmott,t,rnucl,qr,cang,dthetx
      REAL ame,alpha,amp,pmu
      REAL gevfm               !  1/1GeV in fm 
      DATA ame/0.000511/
      DATA amp/0.938/
      DATA gevfm/0.1973/
      DATA pmu/2.79/
C
C     ------------------------------------------------------------------
C
      CSMOTT=-1.
C
C
C      dthetx=1./THETX**2
      IF(ABS(ALMOTT).LT.1.E-6) THEN
         dthetx=THETX
      ELSE
         dthetx=-ALOG(THETX)/ALMOTT
      ENDIF
      costh=DCOS(dthetx)
      alpha=1./137.
C
      e0=CSDMOTT(1)
      ztar=CSDMOTT(2)
      atar=CSDMOTT(3)
      p0=SQRT(e0**2-ame**2)
C
      s=ame**2+amp**2+2.*amp*e0
      sqs=SQRT(s)
      gam=(e0+amp)/sqs
      bet=SQRT(1.-1./gam**2)
      ecme=(s+ame**2-amp**2)/sqs/2.
      pcm=SQRT(ecme**2-ame**2)
C
      IF(costh.LT.-1..OR.costh.GT.1.) THEN
         WRITE(6,*) ' *** Error: costh is out of range =',costh
         GO TO 999
      ENDIF
C
      sinth=DSIN(dthetx)
      sin2=(1.-costh)/2.
      cos2=SQRT(1.-sin2)
C
C---       Final state energies/momenta
C
      p=amp*e0/(e0*(1.-costh)+amp)
      e=SQRT(p**2+ame**2)
      ep=e0+amp-e
      pp=SQRT(ep**2-amp**2)
      sinp=p*sinth/pp
      thetp=ASIN(sinp)
      cosp=SQRT(1.-sinp**2)
C
C---      Q2 ...
C
      q2=-2.*(ame**2-e0*e+p0*p*costh)
C      q2=4.*p0*p*sin2
      t=q2/4./amp**2
C
C---      Formfactors 
C
      cang=0.
      IF(ABS(atar-1.).LT.0.1) THEN
C
C---      Nucleon(dipole formula)
C
         ge=1./(1.+q2/0.71)**2
         gm=ge*pmu
         cang=(ge**2+t*gm**2)/(1.+t)+t*2*gm**2*sin2/cos2
         cang=cang*ztar           ! quasielastic non-coherent scattering (z**1)
      ELSE IF(atar.GT.10.) THEN
C
C---    Heavy nucleus, spin=0 : a ball
C
         rnucl=1.07*atar**0.3333/0.197  ! raduius in GeV**-1, (1.07fm*a**3)
         qr=rnucl*SQRT(ABS(q2))
         ge=3.*ztar/qr**3*(SIN(qr)-qr*COS(qr))
         gm=0.
         cang=ge**2
      ENDIF
C
C---     Cross-section (CM) 
C
      crmott=alpha**2*cos2/4./p0**2/sin2**2
     +       /(1.+2.*p0/amp*sin2)*gevfm**2/100.
C
C      WRITE(6,*) cos2,cang
      cross=crmott*cang
      CSMOTT=cross*sinth
      IF(ABS(ALMOTT).GT.1.E-6) THEN
         CSMOTT=CSMOTT/THETX/ALMOTT
      ENDIF
C
C      WRITE(6,*) '..=',bet,gam,sqs,ecme,pcm
C     +          ,e0,e,p0,p,q2,crmott
C
 999  CONTINUE
      END
