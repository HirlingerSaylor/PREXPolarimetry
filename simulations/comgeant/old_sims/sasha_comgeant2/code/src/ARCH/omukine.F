* $Header:$
* $Log:$
*
      SUBROUTINE OMUKINE(NTGEN)
C
C     ******************************************************************
C     *                                                                *
C     *   Users's program to simulate the kinematics                   *
C     *   Template - 1 particle with a given type, momentum and angles * 
C     *                                                                *
C     *    ==>Called by : OEXCLGEN                                     *
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
#include "omgeant/omceven.inc"
#include "omgeant/omcgen.inc"
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
      INTEGER   i1,i2,i3,i,ip,itrp,nbuf,nt,ifirst
      REAL      amp,chp,tlp,ubuf(20),p
     +         ,ct,ct1,ct2,st,phi,plab(4)
C
      CHARACTER namp*20
      DATA ifirst/1/
      SAVE ifirst
C
C     ------------------------------------------------------------------
C
      IF(ifirst.NE.0) THEN
         ifirst=0
C
C---        Initialization - 1-st call 
C
      ENDIF
C
      NTGEN=0
C
      i1=IDIRDET(1)
      i2=IDIRDET(2)
      i3=IDIRDET(3)
C
C---     Types of particles for beam and target
C
      ip=9  ! GEANT number for pi- 
      itrp=0
      CALL GFPART(ip,namp,itrp,amp,chp,tlp,ubuf,nbuf)  ! get the particle parameters
      IF(itrp.EQ.0) THEN
         WRITE(LOUT,1000) itrp
 1000    FORMAT(' *** OMUKINE error: no particle ',I6)
         GO TO 999
      ENDIF
C     
      p=1.+RNDM(p)*2.     !  momentum 
      ct1=COS(21.*DEGRAD) !  max polar angle: cos  
      ct2=COS(20.*DEGRAD) !  min polar angle: cos  
      ct=ct1+(ct2-ct1)*RNDM(ct)            ! cos of the polar angle
      st=SQRT(1.-ct**2)
      phi=(180.+(1-2.*RNDM(st))*0.1)*DEGRAD  ! azimuthal angle
      plab(i2)=p*st*COS(phi)               ! output momentum, X-projection
      plab(i3)=p*st*SIN(phi)               ! output momentum, Y-projection
      plab(i1)=p*ct                        ! output momentum, Z-projection
      plab(4)=SQRT(p**2+amp**2)
C
      CALL GSKINE(plab(1),ip,IVRTPRIM,0,0,nt)
C
      NTGEN=1
C
 999  RETURN
      END
