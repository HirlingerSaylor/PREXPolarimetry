* $Header:$
* $Log:$
*
      SUBROUTINE OMROTS(V,ROT)
C
C     ******************************************************************
C     *                                                                *
C     *       Fill a rotation matrix V=ROT*V1                          *
C     *  INPUT: V - 3 vector                                           *
C     * OUTPUT: ROT  - rotation matrix (V1 looks along Z)              *
C     *                                                                *
C     *    ==>Called by : kinematics programs                          *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      REAL     V(3),ROT(3,3)
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      INTEGER   i,j
      REAL ptot,pxz
C
C     ------------------------------------------------------------------
C
      ptot=SQRT(V(1)**2+V(2)**2+V(3)**2)
      pxz =SQRT(V(1)**2+V(3)**2)
C!!
C      WRITE(6,FMT='(A8,12F10.5)') ' ptot=  ',ptot,pxz
      IF(ptot.GT.0.AND.pxz.GT.0.) THEN
         ROT(1,1)=V(3)/pxz
         ROT(1,2)=0.
         ROT(1,3)=-V(1)/pxz
         ROT(2,1)=-V(2)*V(1)/pxz/ptot
         ROT(2,2)=pxz/ptot
         ROT(2,3)=-V(2)*V(3)/pxz/ptot
         ROT(3,1)=V(1)/ptot
         ROT(3,2)=V(2)/ptot
         ROT(3,3)=V(3)/ptot
      ELSE
         DO i=1,3
            DO j=1,3
               IF(j.EQ.i) THEN
                  ROT(j,i)=1.
               ELSE
                  ROT(j,i)=0.
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
