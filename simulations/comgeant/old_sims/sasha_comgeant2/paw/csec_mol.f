      SUBROUTINE CSEC_MOL(THELAB1,THETLAB2,PHILAB1,PHILAB2,EB)
C ---  Calculates the Moller rate into the given stereo-angle
C ---  and the particle momnets
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
      REAL     TETMOLL,HSUM
      EXTERNAL TETMOLL,HSUM
C
      INTEGER nbin,idm
      REAL    r0,x1,x2,atarg,abeam,gam,thcm1,thcm2,ctstep,phmx,radeg
     +       ,ctstep,phmx,x1,x2,ratem,rateh,t1,t2,y1,y2,crosec1
     +       ,xl1,xl2,aa,bb
C
      DATA r0/2.818/  ! electron radius in fm
C
C     ------------------------------------------------------------------
C
C---    Moller scattering
C     
      idm=999
      nbin=1200
      atarg=0.511E-5
      abeam=0.511E-5
      radeg=3.1415/180.
      gam=SQRT((ebeam+atarg)/atarg/2.)
C
      thcm1=
      x1=COS(THETMOL(2)*PI/180.)
      x2=COS(THETMOL(1)*PI/180.)
      CALL HBFUN1(idm,'Moller Theta$',nbin,x1,x2,TETMOLL)
      ebeam=EB
      csfac=r0/10./gam/2.       !    sigma=(r0/gam)**2/4.*((4-sin2)/sin2)**2 
      ctstep=(x2-x1)/nbin
      phmx=(PHIMOL(2)-PHIMOL(1))*PI/180.
      crosec=csfac**2*HSUM(IDMOLTET)*ctstep*phmx
      ratem=crosec*0.003*7.87/55.85*6.02*1.E-1*1.E-6/1.602E-19*26.
      rateh=crosec*100.0*0.00708/1.*6.02*1.E-1*1.E-6/1.602E-19*1.
      WRITE(LOUT,1000) cproc,THETMOL(1),THETMOL(2)
     +                   ,PHIMOL(1),PHIMOL(2),crosec
         WRITE(LOUT,1200) cproc,ratem,cproc,rateh
C
      ELSE IF(IOMEXKIN.EQ.2) THEN
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
C            nbin=200
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
C               CALL HCOPY(IDMOTT,102,'    ')
               CALL HDELET(IDMOTT)
               CALL HBFUN1(IDMOTT,'Mott vs exp(-aph*thet)'
     +                 ,nbin,x1,x2,CSMOTT)
               crosec1=HSUM(IDMOTT)*(x2-x1)/nbin*phmx
            ELSE
               IDRADMOT=IDHNEW(-1)
               CALL HBFUN2(IDRADMOT
     +                 ,'Rad. Mott vs exp(-alp*thet)theta vs Eloss'
     +                 ,nbin,x1,x2,nbinl,xl1,xl2,RADMOTT)
            ENDIF
            WRITE(LOUT,1050) cproc,ALMOTT,x1,x2,y1,y2
         ENDIF          
C     
         ratem=crosec*0.003*7.87/55.85*6.02*1.E-1*1.E-6/1.602E-19
         IF(CSDMOTT(2).LT.1.1) ratem=ratem*26.
         rateh=crosec*100.0*0.07080/1.*6.02*1.E-1*1.E-6/1.602E-19
         WRITE(LOUT,1100) cproc,THETMOT(1),THETMOT(2)
     +                   ,PHIMOT(1),PHIMOT(2),crosec,crosec1
         WRITE(LOUT,1200) cproc,ratem,cproc,rateh
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
      ENDIF
C
      IF(LUNDNT.NE.0) CALL ODEGNTIN
C
      RETURN
 1000 FORMAT(/' === ',A16,' cross-section '
     +     ,F7.2,'<theta CM<',F7.2,2X,F7.2,'<phi CM<',F7.2,' deg'
     +     ,' is ',E12.4,' barn per one target electron')
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
      END
C
      REAL FUNCTION TETMOLL(COSTHET)
C
C===      Moller COS(Theta) distrubution
C
      REAL COSTHET,s2
C
      s2=1.-COSTHET**2
      TETMOLL=((4.-s2)/s2)**2
C
      RETURN
      END
