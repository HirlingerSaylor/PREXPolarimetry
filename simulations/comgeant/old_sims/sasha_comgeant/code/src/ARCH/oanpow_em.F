      REAL FUNCTION OANPOW_EM(THETCM)
C
C===     Analyzing power vs angle in CM THECM (degrees)
C===     COMMON/CKINTMP/:
C===     for processes KEMPROC=1 - Moller, 
C                             =2 - Bhabha, 
C                             =3 - 2-photon annihil
C===     Energy (LAB) of the beam: EEMPROC
C
      IMPLICIT NONE
      REAL THETCM
C
      COMMON/CKINTMP/ KEMPROC,EEMPROC,THCMEMPROC
      INTEGER KEMPROC
      REAL EEMPROC,THCMEMPROC
C
C      VECTOR PREAC(8)
C
      REAL    r0,em,am1(2),am2(2),pi,radeg
      DOUBLE PRECISION th,ecmp,vcmp,ecm,elab,anpw
     +                ,d1,d2,dsig,st,ct,st2,ct2,s,u,t,p1,p2,s2,t2,u2
C
      DATA pi/3.14159/
      DATA em/0.511E-3/
C
C      KEMPROC=PREAC(3)+1.1
C      EEMPROC=PREAC(1)
C
      radeg=pi/180.
      th=THETCM*radeg
C      elab=PREAC(1)
      elab=EEMPROC
C
      am1(1)=em
      am1(2)=em
      am2(1)=em
      am2(2)=em
C
      IF(KEMPROC.EQ.3) THEN
         am2(1)=0.
         am2(2)=0.
      ENDIF
C
      ecm=DSQRT(am1(1)**2+am1(2)**2+2.*am1(2)*elab)
      ecmp=ecm/2.
      p1=SQRT(ecmp**2-am1(1)**2)
      p2=SQRT(ecmp**2-am2(1)**2)
      vcmp=DSQRT(ecmp**2-am1(1)**2)/ecmp
C
C      print *,ecm,ecmp,vcmp,gam,bet
      ct=DCOS(th)
      st=DSIN(th)
      ct2=ct**2
      st2=st**2
C
      s=ecm**2
      t=am1(1)**2+am2(1)**2-2.*ecmp**2+2.*p1*p2*ct
      u=am1(1)**2+am2(2)**2-2.*ecmp**2-2.*p1*p2*ct
      s2=s**2
      t2=t**2
      u2=u**2
C
      anpw=0.
      IF(KEMPROC.EQ.1) THEN
         d1=(s2+u2)/t2+(s2+t2)/u2+2.*s2/t/u
         d2=(s2-u2)/t2+(s2-t2)/u2+2.*s2/t/u
         anpw=-d2/d1
      ELSE IF(KEMPROC.EQ.2) THEN
         d1=(s2+u2)/t2+(s2+t2)/u2+2.*s2/t/u
         d2=(s2-u2)/t2+(s2-t2)/u2+2.*s2/t/u
         anpw=d2/d1
      ELSE IF(KEMPROC.EQ.3) THEN
         d1=1.-vcmp**4+vcmp**2*(1-vcmp**2)*st2
         anpw=(vcmp**2*(1.-ct2**2)-d1)/(vcmp**2*(1.-ct2**2)+d1)
      ENDIF
C      write(6,*) ' KEMPROC=',KEMPROC,elab,anpw,THETCM,th
C
      OANPOW_EM=anpw
C
      RETURN
      END
