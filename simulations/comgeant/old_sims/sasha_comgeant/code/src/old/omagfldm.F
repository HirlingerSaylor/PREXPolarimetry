* $Header:$
* $Log:$
*
      SUBROUTINE OMAGFLDM(MAG,XX,BX)
C
C     *****************************************************************
C     *                                                               *
C     *   Define the magnetic field for simple dipoles and quadrupoles*
C     *   (Used for the Moller polarimeter)                           *
C     *                                                               *
C     *   called from: GUFLD                                          *
C     *                                                               *
C     *****************************************************************
C
      IMPLICIT NONE
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
      INTEGER MAG
      REAL XX(3),BX(3)
C
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omcsvol.inc"
#include "moller/cmolmag.inc"
C
      INTEGER i,molm,i1,i2,i3,ifl,j2,j3
      REAL    rquad,r,z,bzr(2),cphi,sphi,b0,sca,gap,gaph
C
      SAVE b0
      DATA b0/0./
C
C     ------------------------------------------------------------------
C
      DO i=1,3
         BX(i)=0.
      END DO
      IF(IDEBUG.NE.0.AND.ISWIT(3).EQ.13) THEN
         WRITE(LOUT,1200) MAG,XX
      ENDIF
      IF(MAG.LE.0.OR.MAG.GT.NMOLMAG) GO TO 999
C
      molm=MAG
C
      i1=IDIRDET(1)  ! direction along the beam
      i2=IDIRDET(2)  ! the 1-st transverse projection
      i3=IDIRDET(3)
C
      IF(MOLMAGP(1,molm).EQ.1) THEN
C
C---      A quadrupole
C
         DO i=1,3
            IF(ABS(XX(i)).GT.PMOLMAGP(1+i,molm)/2.) GO TO 999
         END DO
C
         rquad=PMOLMAGP(2,molm)/2.
         IF(rquad.GT.0) THEN
            IF(XX(i2)**2+XX(i3)**2.LT.rquad**2) THEN
               BX(i2)=-XX(i3)*AMAGPAR(5,MAG)/rquad
               BX(i3)=-XX(i2)*AMAGPAR(5,MAG)/rquad
               BX(i1)= 0.
            ENDIF
         ENDIF
C
      ELSE IF(IABS(MOLMAGP(1,molm)).EQ.2) THEN
C
C---      A dipole, flag= 2, aprox Bx=B, By=Bz=0
C---                    =-2        By=B, Bz=Bx=0
C
         IF(IDEBUG.NE.0.AND.ISWIT(3).EQ.13) THEN
            WRITE(LOUT,1300) MOLMAGP(1,molm)
     +               ,(PMOLMAGP(i+1,molm),i=1,3)
         ENDIF
         IF(ABS(XX(i2)).GT.PMOLMAGP(1+i2,molm)/2.) GO TO 999
         IF(ABS(XX(i3)).GT.PMOLMAGP(1+i3,molm)/2.) GO TO 999
C
         j2=i2 ! - field direction
         j3=i3 ! - perpend to the beam and to the field
         IF(MOLMAGP(1,molm).EQ.-2) THEN
            j2=i3
            j3=i2
         ENDIF
C
         ifl=0
         IF(ABS(AMAGPAR(6,MAG)).GT.1.E-3) ifl=1
         IF(ifl.EQ.0) THEN
C
C---        A simple rectangular field
C            
            IF(ABS(XX(i1)).GT.PMOLMAGP(1+i1,molm)/2.) GO TO 999
            BX(j2)=AMAGPAR(5,MAG)
            BX(j3)=0.
            BX(i1)=0.
         ELSE IF(ifl.EQ.1) THEN
C
C---        Stray field in the longit. direction at the entrance and exit,
C---        Gradual reduction of the main field.
C---        The length of the area with variable field: -GAP: +GAP from the edge
C---        Longitudinal field Bz=coeff*x/(gap/2)*Bx, at the enterance and -Bz at the exit
C
            gap=PMOLMAGP(1+j2,molm)
            gaph=gap/2.
            IF(ABS(XX(i1)).GT.PMOLMAGP(1+i1,molm)/2.+gaph) THEN
               GO TO 999
            ELSE IF(ABS(XX(i1)).GT.PMOLMAGP(1+i1,molm)/2.-gaph) THEN
               BX(j2)=AMAGPAR(5,MAG)
     +               *(PMOLMAGP(1+i1,molm)/2.+gaph-ABS(XX(i1)))/gap   ! B linearly goes to 0
               BX(j3)=0.
C               BX(i1)=-AMAGPAR(5,MAG)*XX(i2)/gap*AMAGPAR(6,MAG)        ! Longit. component
               BX(i1)=-AMAGPAR(5,MAG)*AMAGPAR(6,MAG)*XX(j2)/gaph        ! Longit. component X-dependence
     +               *(1.-(ABS(XX(i1))-PMOLMAGP(1+i1,molm)/2.)/gaph)**2 !      longit dependence
     +                  *SIGN(1.,XX(i1))                                !      sign
            ELSE
               BX(j2)=AMAGPAR(5,MAG)
               BX(j3)=0.
               BX(i1)=0.
            ENDIF
C           
         ENDIF
C
      ELSE IF(MOLMAGP(1,molm).EQ.3) THEN
C
C---      A solenoid
C
         IF(ABS(XX(i1)).GT.PMOLMAGP(1+i1,molm)/2.) GO TO 999
C
         IF(SQRT(XX(i2)**2+XX(i3)**2).LT.PMOLMAGP(1+i2,molm)/2.) THEN
            BX(i1)=AMAGPAR(5,MAG)
            BX(i2)=0.
            BX(i3)=0.
         ENDIF
C
      ELSE IF(MOLMAGP(1,molm).EQ.4) THEN
C
C---      A solenoid..(axially symm field, the axis along Z) 
C---           with a field map provided
C
         z=ABS(XX(i1))
         r=SQRT(XX(i2)**2+XX(i3)**2)
         IF(r.GT.0.) THEN
            cphi=XX(i2)/r
            sphi=XX(i3)/r
         ELSE
            cphi=1.
            sphi=0.
         ENDIF
C---        Get the value from the map
         IF(b0.EQ.0.) THEN
            CALL OMSOLMAP(molm,0.,0.,bzr(1))
            b0=SQRT(bzr(1)**2+bzr(2)**2)
         ENDIF
         CALL OMSOLMAP(molm,z,r,bzr(1))
C
         sca=1.
         IF(b0.GT.0.) sca=AMAGPAR(5,MAG)/b0
C         
         BX(i1)=bzr(1)*sca
         BX(i2)=bzr(2)*sca*cphi
         BX(i3)=bzr(2)*sca*sphi
C
      ENDIF
C
 999  RETURN
 1200 FORMAT(' --- OMAGFLDM called for MAG=',I3,' at X=',3F10.2)
 1300 FORMAT(' --- OMAGFLDM dipole type=',I3,' sizes',3F10.2)
      END





