      REAL FUNCTION SOLEN_F(XX,YY,ZZ)
*
* ===   Solenoid field
*
      IMPLICIT NONE
      REAL XX,YY,ZZ
C
      COMMON/CSOLEN/ PARSOL(12,3),BSOL(3)
      REAL           PARSOL,BSOL
C
      INTEGER nz,nr,nphi,iz,ir,iphi,i,nsol,isol
C
      REAL x0(3),x(3),d(3),dl(3),b(3),dm,dz,dr,dphi,curdens,zl,r1,r2
     +    ,rr,phi,bb,b1,db(3),pi,xc(3)
C
C     -----------------------------------------------------------------
C
      nz=50
      nr=4
      nphi=100
      nsol=3
      pi=ACOS(0.)*2.
C
      DO i=1,3
         b(i)=0.
      ENDDO
      x0(1)=XX
      x0(2)=YY
      x0(3)=ZZ
C      WRITE(6,*) x0
C
      DO isol=1,nsol
C
         IF(PARSOL(1,isol).GT.0.1) THEN
C
C---            Get the soledoind parameters
C
            DO i=1,3
               xc(i)=PARSOL(i+1,isol)
            ENDDO
            r1=PARSOL(5,isol)
            r2=PARSOL(6,isol)
            zl=PARSOL(7,isol)
            curdens=PARSOL(8,isol)
C
            dz=zl/nz*2.
            dr=(r2-r1)/nr
            dphi=2.*pi/nphi
C            WRITE(6,*) x0,r1,r2,zl,curdens
C
C--- Formula:         dB = mu0/4pi * d x dl /d**3  (muo=4pi*1.e-7)
C
            b1=1.E-7*curdens*dz*dr*dphi         
            b1=b1*10.*100.        ! convert to kGs
C
            DO iz=1,nz
               x(3)=-zl+(iz-0.5)*dz
C
               DO ir=1,nr
                  rr=r1+(ir-0.5)*dr
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
                        d(i)=x0(i)-xc(i)-x(i)
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
               ENDDO
            ENDDO
C
         ENDIF
      ENDDO
C
      SOLEN_F=SQRT(b(1)**2+b(2)**2+b(3)**2)
      DO i=1,3
         BSOL(i)=b(i)
      ENDDO
C
      RETURN
      END
