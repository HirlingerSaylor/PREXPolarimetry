* $Header:$
* $Log:$
*
      SUBROUTINE OMOTT(NTGEN)
C
C     ******************************************************************
C     *                                                                *
C     *       Kinematics for Mott scattering                           *
C     *                                                                *
C     *    ==>Called by : OMKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER NTGEN
C
#include "geant321/gconst.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcbeam.inc"
#include "omgeant/omceven.inc"
#include "omgeant/omcgen.inc"
#include "moller/cescatt.inc"
#include "moller/cntescat.inc"
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      COMMON/CCSMOTT/ CSDMOTT(4),ALMOTT
      REAL            CSDMOTT,ALMOTT
C---                  CSMOTT(1:3) = E0 electron (LAB), Z target, A target.
C---                  ALMOTT      - exponential slope of the variable
C
      REAL     HRNDM1,RNDM,CSMOTT
      EXTERNAL HRNDM1,RNDM,CSMOTT
C
      INTEGER   nbuf,itrtbm,itrttg,itrtut
     +         ,it,nt,kt
     +         ,i,j,iptypp(3),irmot,ifermi
     +         ,ntry
      REAL      plab0(4,3)           ! initial, LAB frame
     +         ,plab(4,2)            ! final  , LAB     
     +         ,pref0(4,2)           ! initial, in CM of the target
     +         ,pref01(4)            ! initial, in CM of the target, rotated that Pbeam is along Z'
     +         ,pref(4,3)            ! final,   in CM of the target, rotated that Pbeam is along Z'
     +         ,pref1(4,3)           ! final,   in CM of the target
     +         ,chbeam,tlbeam,atarg,chtarg,tltarg
     +         ,ab2,at2,pb2,eb2
     +         ,abeam
     +         ,phi,tet
     +         ,ct,st,sp,cp
     +         ,ebeam,xt,t1,t2,ubuf(10),elos
     +         ,bet(4),bet1(4),efermi,pfermi,rot(3,3),cf,sf,pf
     +         ,ptot,pxz,pp
      CHARACTER nambeam*20,namtarg*20
C
C     ------------------------------------------------------------------
C
      NTGEN=0
C
C
      ITYPRS=ITYPREST
C
C---     Types of particles for beam and target
C
      CALL GFPART(IBEAMPAR,nambeam,itrtbm,abeam,chbeam,tlbeam,ubuf,nbuf)
      CALL GFPART(MOTTARG,namtarg,itrttg,atarg,chtarg,tltarg,ubuf,nbuf)
      iptypp(1)=IBEAMPAR
      iptypp(2)=MOTTARG
      iptypp(3)=1    !  photon - in case of Rad. Mott.
      ab2=abeam**2
      at2=atarg**2
C
      DO i=1,3
        plab0(i,1)=PBEAM(i)
        plab0(i,2)=0.
      ENDDO
      plab0(4,1)=SQRT(abeam**2+PBEAM(1)**2+PBEAM(2)**2+PBEAM(3)**2)
      plab0(4,2)=atarg
      DO i=1,4
         pref01(i)=plab0(i,1)
      ENDDO
C
C---     Fermi motion (quasielastic scattering): change the reference frame that
C        the target proton is at rest 
C
      ifermi=0
      IF(ABS(atarg-0.938).LT.0.5.AND.AMATTAR.GT.1.5.AND.
     +   FERMIMOM.GT.0.00001) THEN
         IF(AMATTAR.GT.1.5) THEN
C
C
C
C
C---          Generate the target momentum
C---          Assumed Fermi distribution: dn/(dpx dpy dpz) = const(p), or dn/dp = p**2*const 
C
            pfermi=RNDM(pfermi)**0.333333*FERMIMOM
            cf=2.*RNDM(cf)-1.
            sf=SQRT(1.-cf**2)
            pf=TWOPI*RNDM(sf)
C!!
C            pfermi=FERMIMOM
C            pf=0.
C            cf=1.
C            sf=0.
C!!
            plab0(1,2)=pfermi*sf*COS(pf)
            plab0(2,2)=pfermi*sf*SIN(pf)
            plab0(3,2)=pfermi*cf
            plab0(4,2)=SQRT(pfermi**2+at2)
C
C---          Define the frame where the target is at rest
C
            DO i=1,3
               bet(i)=plab0(i,2)/plab0(4,2)
               bet1(i)=-bet(i)
            ENDDO
            bet(4)=plab0(4,2)/atarg
            bet1(4)=bet(4)
C
C---          Boost the beam into new frame
C
            CALL GLOREN(bet(1),plab0(1,1),pref0(1,1))
C
C---          Rotation matrix (from the rotated frame (momentum=Z) to LAB)
C---           Rotation: turn around Y, then turn around X' 
C
            ptot=SQRT(pref0(1,1)**2+pref0(2,1)**2+pref0(3,1)**2)
            pxz =SQRT(pref0(1,1)**2+pref0(3,1)**2)
            IF(ptot.GT.0.AND.pxz.GT.0.) THEN
               ifermi=1
               rot(1,1)=pref0(3,1)/pxz
               rot(1,2)=0.
               rot(1,3)=-pref0(1,1)/pxz
               rot(2,1)=-pref0(2,1)*pref0(1,1)/pxz/ptot
               rot(2,2)=pxz/ptot
               rot(2,3)=-pref0(2,1)*pref0(3,1)/pxz/ptot
               rot(3,1)=pref0(1,1)/ptot
               rot(3,2)=pref0(2,1)/ptot
               rot(3,3)=pref0(3,1)/ptot
            ENDIF
            pref01(1)=0.
            pref01(2)=0.
            pref01(3)=ptot
            pref01(4)=SQRT(ptot**2+ab2)
C
         ENDIF
      ENDIF
C!!
C      WRITE(6,*) 'event'
C      WRITE(6,FMT='(A8,12F10.5)') ' plab0= ',plab0
C      WRITE(6,FMT='(A8,12F10.5)') ' pref0= ',pref0
C      WRITE(6,FMT='(A8,12F10.5)') ' pref01=',pref01
C      WRITE(6,FMT='(A8,12F10.5)') ' bet=   ',bet
C      WRITE(6,FMT='(A8,12F10.5)') ' bet1=  ',bet1
C      WRITE(6,FMT='(A8,12F10.5)') ' rot=   ',rot
C
      pb2=pref01(1)**2+pref01(2)**2+pref01(3)**2
      eb2=pb2+ab2
      ebeam=SQRT(eb2)
C
      CSDMOTT(1)=ebeam
      CSDMOTT(2)=chtarg
      CSDMOTT(3)=atarg
      ALMOTT=0.
C
C---    Mott production in LAB
C
      phi=PI/180.*(RNDM(ebeam)*(PHIMOT(2)-PHIMOT(1))+PHIMOT(1))
      sp=SIN(phi)
      cp=COS(phi)
C
      t1=PI/180.*THETMOT(1)
      t2=PI/180.*THETMOT(2)
C
      WGCROS=1.
      irmot=0
      IF(IMOTTWG.EQ.2) THEN
         tet=RNDM(cp)*(t2-t1)+t1
         WGCROS=CSMOTT(tet)
      ELSE
C 
         IF(IMOTTWG.EQ.1) THEN
            ALMOTT=ALPHMOTT
            xt=RNDM(cp)*(EXP(-ALMOTT*t1)-EXP(-ALMOTT*t2))
     +                 +EXP(-ALMOTT*t2)
            WGCROS=CSMOTT(xt)
         ELSE
            WGCROS=((BEAMOLIM(1)+BEAMOLIM(2))/2./ebeam)**2
            IF(IRADMOT.EQ.0) THEN
               xt=HRNDM1(IDMOTT)
               elos=1.
            ELSE
               CALL HRNDM2(IDRADMOT,xt,elos)
               irmot=1                         ! radiation after scattering
            ENDIF
         ENDIF
C
         IF(ALPHMOTT.LE.0) THEN
            WRITE(LOUT,*) '*** OMOTT error: ALPHMOTT=',ALPHMOTT,' abort'
            IEORUN=1
            ALMOTT=1.
         ENDIF
         tet=-ALOG(xt)/ALPHMOTT
      ENDIF
C
C!!
C      tet=0.
C
      st=SIN(tet)
      ct=COS(tet)
C
C---   Find the energies/momenta in the Reference frame frame
C
C---       e-
C
      ptot=ebeam*atarg/(atarg+ebeam*(1.-ct))
C
C---       Radiation
C
      IF(irmot.GT.0) THEN
         ptot=ptot*elos
      ENDIF
      pref(4,1)=SQRT(ptot**2+ab2)
      pref(1,1)=ptot*st*cp
      pref(2,1)=ptot*st*sp
      pref(3,1)=ptot*ct
C!!
C      WRITE(6,FMT='(A8,12F10.5)') ' pref=  ',pref
C      WRITE(6,FMT='(A8,12F10.5)') ' ptot..=',ptot,ebeam,ct
C
C---       Radiation: define the photon (along incoming e- - Z' or the outgoing e-)
C
      IF(irmot.GT.0) THEN
         IF(RNDM(ptot).LT.0.4) THEN
            pref(1,3)=0.
            pref(2,3)=0.
            pref(3,3)=ptot*(1.-elos)          ! along Z'
         ELSE
            DO i=1,3
               pref(i,3)=pref(i,1)*(1.-elos)  ! along e-
            ENDDO
         ENDIF
         pref(4,3)=SQRT(pref(1,3)**2+pref(2,3)**2+pref(3,3)**2)
      ENDIF
C
C---        Recoil
C
      DO i=1,3
         pref(i,2)=pref01(i)-pref(i,1)
         IF(irmot.GT.0) pref(i,2)=pref(i,2)-pref(i,3)
      ENDDO
      pref(4,2)=SQRT(pref(1,2)**2+pref(2,2)**2+pref(3,2)**2+at2)
C
C---        Turn back to X,Y,Z in the frame with the initial target at resr 
C
      kt=2
      IF(irmot.NE.0) kt=3
      DO it=1,kt
         IF(ifermi.NE.0) THEN                   ! Fermi motion: rotate the momenta and boost to the LAB
            DO i=1,3
               pp=0.
               DO j=1,3
                  pp=pp+rot(j,i)*pref(j,it)
               ENDDO
               pref1(i,it)=pp
            ENDDO
            pref1(4,it)=pref(4,it)
            CALL GLOREN(bet1(1),pref1(1,it),plab(1,it))
         ELSE                                   ! No Fermi motion - use the momenta as they are 
            DO i=1,4
               plab(i,it)=pref(i,it)
            ENDDO
         ENDIF
      ENDDO
C!!
C      WRITE(6,FMT='(A8,12F10.5)') ' pref=  ',pref
C      WRITE(6,FMT='(A8,12F10.5)') ' pref1= ',pref1
C      WRITE(6,FMT='(A8,12F10.5)') ' plab=  ',plab
C
C---      NTUPLE variables
C
      IETYPE=2
      THETCM=tet
      PHICM=phi
      ANPOWER=0.
C
      DO it=1,kt
         CALL GSKINE(plab(1,it),iptypp(it),IVRTPRIM,0,0,nt)
      END DO
C
      NTGEN=2
C
C      WRITE(6,*) 'pcm=',pcm
C
      RETURN
      END




