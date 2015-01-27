* $Header:$
* $Log:$
*
      SUBROUTINE OPHOTBRM(NTGEN)
C
C     ******************************************************************
C     *                                                                *
C     *       Kinematics for pion, kaon... photoproduction             *
C     *       from Bremstralung of an electron beam in the target      *
C     *       using Bosted - Wiser program/distributions               *
C     *                                                                *
C     *    ==>Called by : OEXCLGEN                                     *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER NTGEN
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
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      REAL     RNDM
      EXTERNAL RNDM
C
      INTEGER   i1,i2,i3,nt,ip
      REAL      plab(4,2)
     +         ,rmom,pbm,ppar
     +         ,phi,tet

C     ------------------------------------------------------------------
C
      NTGEN=0
C
      i1=IDIRDET(1)
      i2=IDIRDET(2)
      i3=IDIRDET(3)
C
      ITYPRS=ITYPREST
C
C---     Types of particles for beam and target
C
      pbm=SQRT(PBEAM(1)**2+PBEAM(2)**2+PBEAM(3)**2)
C
      CALL HRNDM2(IDPHOTPR,tet,rmom)
      ppar=rmom*pbm
      phi=PI*2.*RNDM(tet)
      plab(i1,1)=ppar*COS(tet)
      plab(i2,1)=ppar*SIN(tet)*COS(phi)
      plab(i2,1)=ppar*SIN(tet)*SIN(phi)
C
      ip=5
      IF(IPHOTPR.GT.0.AND.IPHOTPR.LE.6) THEN
         ip=IPHPRGEA(IPHOTPR)
      ENDIF
C
C---      NTUPLE variables
C
      IETYPE=3
      THETCM=tet
      PHICM=phi
      ANPOWER=0.
C
      CALL GSKINE(plab(1,1),ip,IVRTPRIM,0,0,nt)
C
      NTGEN=1
C
C      WRITE(6,*) 'pcm=',pcm
C
      RETURN
      END


