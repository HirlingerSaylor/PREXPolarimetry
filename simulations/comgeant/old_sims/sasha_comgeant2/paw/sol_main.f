      PROGRAM SOL_MAIN
      IMPLICIT NONE
C
C---     Calculate a solenoid field and write out the results
C
      COMMON/CSOLEN/ PARSOL(12,3),BSOL(3)
      REAL           PARSOL,BSOL
C
      REAL SOLEN_F
C
      INTEGER    mxpnt
      PARAMETER (mxpnt=100000)
      COMMON/CSOLRES/ NPNT,PZX(2,mxpnt),BZX(2,mxpnt)
      INTEGER npnt
      REAL    PZX,BZX
C
      INTEGER nx,nz,ix,iz,i,j
      REAL zx(2,2),dzx(2),x,z,b
C
      DATA PARSOL/1.,0.,0.,0.  ,92. ,102., 70.,-400.,4*0.
     +           ,0.,0.,0.,200.,14. , 17., 70.,-400.,4*0.
     +           ,12*0./
C
C     ------------------------------------------------------------------
C
      dzx(1)=1.
      dzx(2)=1.
      zx(1,1)=-2.*dzx(1)
      zx(2,1)=-2.*dzx(2)
      zx(1,2)=200.
      zx(2,2)=120.
      nz=(zx(1,2)-zx(1,1))/dzx(1)+1
      nx=(zx(2,2)-zx(2,1))/dzx(2)+1
C
      NPNT=0
C
      DO ix=1,nx
         x=zx(2,1)+dzx(2)*(ix-1)
         DO iz=1,nz
            z=zx(1,1)+dzx(1)*(iz-1)
            IF(NPNT.LT.mxpnt) THEN
               NPNT=NPNT+1
            ELSE
               WRITE(6,*) ' *** Error: too many poins: abort'
               GO TO 999
            ENDIF
C
            b=SOLEN_F(x,0.,z)
            PZX(1,NPNT)=z
            PZX(2,NPNT)=x
            BZX(1,NPNT)=BSOL(3)
            BZX(2,NPNT)=BSOL(1)
C
         ENDDO
         IF(MOD(ix,10).EQ.0) WRITE(6,*) ' ix,npnt=',ix,npnt 
      ENDDO
C
      WRITE(45) NPNT,nz,nx,zx,dzx
     +         ,((PZX(j,i),j=1,2),(BZX(j,i),j=1,2),i=1,NPNT)
C
 999  CONTINUE
      END
C
      INCLUDE 'solen_f.f'



