      REAL FUNCTION COIL_F(XX)
*
* ===   The field B (Gs) from a coil (1 turn) of a thin wire, in the coil plane (XY)
*        XX - the radius of the point, R coil=1cm, I=1A
*
      IMPLICIT NONE
      REAL XX
C
      INTEGER nphi,iphi,i
C
      REAL x0(3),x(3),d(3),dl(3),b(3),dphi,cur,r1,dm
     +    ,rr,phi,bb,b1,db(3),pi
C
C     -----------------------------------------------------------------
C
      nphi=100
      pi=ACOS(0.)*2.
C
      DO i=1,3
         b(i)=0.
      ENDDO
      x0(1)=XX
      x0(2)=0.
      x0(3)=0.
C      WRITE(6,*) x0
C
      r1=1.
      cur=1.
C
      dphi=2.*pi/nphi
C            WRITE(6,*) x0,r1,r2,zl,curdens
C
C--- Formula:         dB = mu0/4pi * d x dl /d**3  (muo=4pi*1.e-7)
C
      b1=1.E-7*cur*dphi         
      b1=b1*10.*100.*1000.            ! convert to Gs
C
      x(3)=0.
C
      rr=r1
C
      DO iphi=1,nphi
         phi=(iphi-0.5)*dphi
C
         dl(1)=-SIN(phi)
         dl(2)=COS(phi)
         dl(3)=0.
C
         x(1)=rr*dl(2)
         x(2)=rr*(-dl(1))
         DO i=1,3
            d(i)=x0(i)-x(i)
         ENDDO
         dm=SQRT(d(1)**2+d(2)**2+d(3)**2)
C
         bb=b1/dm**2/dm*rr
         db(1)=bb*(d(2)*dl(3)-d(3)*dl(2))
         db(2)=bb*(d(3)*dl(1)-d(1)*dl(3))
         db(3)=bb*(d(1)*dl(2)-d(2)*dl(1))
C
         DO i=1,3
            b(i)=b(i)+db(i)
         ENDDO
      ENDDO
C
C      COIL_F=SQRT(b(1)**2+b(2)**2+b(3)**2)
      COIL_F=b(3)
C
      RETURN
      END






