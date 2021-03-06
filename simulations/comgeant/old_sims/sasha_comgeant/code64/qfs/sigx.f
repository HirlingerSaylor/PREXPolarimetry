
*SIGX
      REAL FUNCTION SIGX(E,TH,W,A)
C      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT NONE
      REAL E,TH,W,A
      REAL FPHENOM
      EXTERNAL FPHENOM
      
      REAL alph,pi,sig0,sig1,pimass,pm,gam0,aq,thr,qms,arg0,arg1,arg
      REAL shape,ekappa,siggam,qs,eps,flux,sigee,r,factor1

      write(6,*) 'SIGX',E,TH,W,A
      SIGEE=0.
      ALPH=1./137.03604
      PI=ACOS(-1.)
C     SIG0=111.*1.E-30
      SIG0=100.D-4
C     SIG1=60.*1.E-27
      SIG1=54.*1.D-1
      PIMASS=140.
      PM=939.
C     GAM0=550.
      GAM0=650.
C     R=0.10
      AQ=250.
      THR=TH*PI/180.
      IF(W.LT.1.E-5)GO TO 4
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      ARG0=W-QMS/2./PM-PIMASS-PIMASS**2/2./PM
      ARG1=ARG0/GAM0
      ARG=ARG1**2/2.
      IF(ARG1.GT.8.)THEN
      SHAPE=1.+SIG1/SIG0/ARG0
      ELSEIF(ARG1.LT.1.E-5)THEN
      SHAPE=0.
      ELSEIF(ARG1.LT.0.1)THEN
      SHAPE=SIG1*ARG0/2./GAM0**2/SIG0
      ELSE
      SHAPE=(1.-EXP(-ARG))*(1.+SIG1/SIG0/ARG0)
      ENDIF
      EKAPPA=W-QMS/2./PM
      SIGGAM=SIG0*SHAPE
      QS=QMS+W**2
      EPS=1./(1.+2.*QS*TAN(THR/2.)**2/QMS)
      FLUX=ALPH*EKAPPA*(E-W)/2./PI**2/QMS/E/(1.-EPS)
      IF(FLUX.LT.1.E-20)FLUX=0.
      SIGEE=FLUX*SIGGAM*FPHENOM(QMS)**2
C     SIGEE=FLUX*SIGGAM
      R=0.56*1.E6/(QMS+PM**2)
      FACTOR1=1.+EPS*R
      SIGEE=SIGEE*FACTOR1
 4    SIGX=A*SIGEE
      write(6,*) 'SIGX',SIGX
      RETURN
      END
