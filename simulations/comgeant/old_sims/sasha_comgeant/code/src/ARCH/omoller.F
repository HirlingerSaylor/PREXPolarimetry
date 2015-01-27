* $Header:$
* $Log:$
*
      SUBROUTINE OMOLLER(NTGEN,IFL)
C
C     ******************************************************************
C     *                                                                *
C     *       Kinematics Moller/Bhabha/Annihilatin: 2-body reaction    *
C     *                                                                *
C     *    ==>Called by : OMKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER NTGEN,IFL
C
#include "geant321/gconst.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omcbeam.inc"
#include "omgeant/omceven.inc"
#include "omgeant/omcgen.inc"
#include "omgeant/cescatt.inc"
#include "moller/cntescat.inc"
C
      COMMON/CKINTMP/ KEMPROC,EEMPROC,THCMEMPROC
      INTEGER KEMPROC
      REAL EEMPROC,THCMEMPROC
C
      CHARACTER atshell*1,namshell*3
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      REAL     OANPOW_EM,HRNDM1,RNDM
      EXTERNAL OANPOW_EM,HRNDM1,RNDM
C
      INTEGER   itart,ntrymx,i1,i2,i3
     +         ,nbuf,itrtbm,itrttg,itrtut
     +         ,it,nt
     +         ,i,iptypp(2),irotat,itrt(2)
      REAL      eecm,pecm,pbt,pbm
     +         ,pcm(4,2),p(4),plab(4,2),betal(4),ubuf(20)
     +         ,chbeam,tlbeam,atarg,chtarg,tltarg
     +         ,eb2,e0
     +         ,phi,tet
     +         ,ct,st,sp,cp
     +         ,ctr,str,spr,cpr
     +         ,rm(2,2),chp(2),tlp(2),rp,re
C
      DOUBLE PRECISION debm,dpbm2,dpbm,debm2
     +                ,delab,dpb2,deb2,decm
     +                ,dm(2,2),dm2(2,2),dp,de
C
      CHARACTER nambeam*20,namtarg*20,namsec(2)*20
C
      SAVE itart,ntrymx
      DATA itart/3/
      DATA ntrymx/50/
      DATA namshell/'ZLK'/
C
C     ------------------------------------------------------------------
C
      NTGEN=0
C
      i1=IDIRDET(1)
      i2=IDIRDET(2)
      i3=IDIRDET(3)
C
      IF(IMSHELL.GE.0.AND.IMSHELL.LE.2) THEN
         atshell=namshell(IMSHELL+1:IMSHELL+1)
      ELSE
         WRITE(LOUT,*) ' *** OMOLLER: atomic shell number IMSHELL=',
     +                IMSHELL,' is out of range'
         atshell='Z'
      ENDIF
C
      ITYPRS=ITYPREST
C
C---     Types of particles for beam and target
C
      CALL GFPART(IBEAMPAR,nambeam,itrtbm,rm(1,1),chbeam,tlbeam
     +           ,ubuf,nbuf)
      CALL GFPART(   itart,namtarg,itrttg,rm(2,1),chtarg,tltarg
     +           ,ubuf,nbuf)
      IF(IFL.EQ.0) THEN  ! Bhabha/Moller
         iptypp(1)=IBEAMPAR
         iptypp(2)=itart
         IF(IBEAMPAR.EQ.2) THEN
            KEMPROC=2   ! Bhabha
         ELSE IF(IBEAMPAR.EQ.1) THEN
            KEMPROC=1   ! Moller
         ENDIF
      ELSE IF(IFL.EQ.1) THEN ! Annihilation
         iptypp(1)=1
         iptypp(2)=1
         KEMPROC=3
      ELSE
         WRITE(LOUT,*) ' *** OMOLLER error: undefined process IFL=',IFL
         IEORUN=1
         RETURN
      ENDIF
      DO it=1,2
         CALL GFPART(iptypp(it),namsec(it),itrt(it),rm(it,2),chp(it)
     +           ,tlp(it),ubuf,nbuf)
      ENDDO
C
      DO it=1,2
         DO i=1,2
            dm(i,it)=DBLE(rm(i,it))
            dm2(i,it)=dm(i,it)**2
         ENDDO
      ENDDO
C
      dpbm2=DBLE(PBEAM(1))**2+DBLE(PBEAM(2))**2+DBLE(PBEAM(3))**2
      dpbm=DSQRT(dpbm2)
      pbm=REAL(dpbm)
      debm2=dpbm2+dm2(1,1)
      eb2=REAL(debm2)
C
      EEMPROC=SQRT(eb2)
C
      delab=DSQRT(debm2)+dm(2,1)
      decm=DSQRT(delab**2-dpbm2)
      de=(decm**2+dm(1,2)**2-dm(2,2)**2)/2./decm ! secondary particle energy
      dp=DSQRT(de**2-dm2(1,2)) ! momentum
      re=REAL(de)
      rp=REAL(dp)
C
C      pb2=PBEAM(1)**2+PBEAM(2)**2+PBEAM(3)**2
C      eb2=pb2+ab2
C      elab=SQRT(eb2)+atarg
C      elab2=elab**2
C      tecm2=elab2-pb2
C      tecm=SQRT(tecm2)
C
C---      Lorentz boost vector
C
      betal(i2)=0.
      betal(i3)=0.
      betal(i1)=REAL(-dpbm/delab)
      betal(4)=REAL(delab/decm)
C
C---      Rotation matrix
C
      irotat=0
      pbt=SQRT(PBEAM(i2)**2+PBEAM(i3)**2)
      IF(pbt.GT.pbm*1.E-6) THEN
        irotat=1
        str=pbt/pbm
        ctr=SQRT(1.-str**2)
        cpr=PBEAM(i2)/pbt
        spr=PBEAM(i3)/pbt
      ENDIF

C
C---    2-body production in CM, Z - along the beam
C
      phi=PI/180.*(RNDM(eecm)*(PHIMOL(2)-PHIMOL(1))+PHIMOL(1))
      sp=SIN(phi)
      cp=COS(phi)
      ct=HRNDM1(IDMOLTET)
      st=SQRT(1.-ct**2)
      tet=ACOS(ct)
C
C---      NTUPLE variables
C
      IETYPE=1
      THETCM=tet
      PHICM=phi
      WGCROS=1.
C      ANPOWER=st**2*(8.-st**2)/(4.-st**2)**2
      ANPOWER=OANPOW_EM(ACOS(ct)*180./3.14159)
C
      pcm(i2,1)=rp*st*cp
      pcm(i3,1)=rp*st*sp
      pcm(i1,1)=rp*ct
      pcm(4,1)=re
C
      DO i=1,3
         pcm(i,2)=-pcm(i,1)
      END DO
      pcm(4,2)=REAL(decm-de)
C
C---     Rotation
C
      IF(irotat.NE.0) THEN
         DO it=1,2
            DO i=1,3
               p(i)=pcm(i,it)
            ENDDO
            pcm(i1,it)=p(i1)*ctr    -p(i2)*str
            pcm(i2,it)=p(i1)*str*cpr+p(i2)*ctr*cpr-p(i3)*spr
            pcm(i3,it)=p(i1)*str*spr+p(i2)*ctr*spr+p(i3)*cpr
         ENDDO
      ENDIF

         
      IF(atshell.NE.'K'.AND.atshell.NE.'L') THEN
         DO it=1,2
            CALL GLOREN(betal,pcm(1,it),plab(1,it))
         END DO
      ELSE
         e0=SQRT(eb2)
         CALL OATOMIC(atshell,e0,ct,st,cp,sp,plab)
      ENDIF
      DO it=1,2
         CALL GSKINE(plab(1,it),iptypp(it),IVRTPRIM,0,0,nt)
      END DO
C
      NTGEN=2
C
C      WRITE(6,*) 'pcm=',pcm
C
      RETURN
      END
