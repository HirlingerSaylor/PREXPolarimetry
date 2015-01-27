      PROGRAM RD_MAIN
      IMPLICIT NONE
C
C---     Read a field map
C
      INTEGER    mxpnt
      PARAMETER (mxpnt=100000)
      COMMON/CSOLRES/ NPNT,PZX(2,mxpnt),BZX(2,mxpnt)
      INTEGER npnt
      REAL    PZX,BZX
C
      INTEGER nx,nz,ix,iz,i,j
      REAL zx(2,2),dzx(2)
C
C     ------------------------------------------------------------------
C
      NPNT=0
C
      READ(45) NPNT,nz,nx,zx,dzx
     +        ,((PZX(j,i),j=1,2),(BZX(j,i),j=1,2),i=1,NPNT)
      WRITE(6,1000) NPNT,nz,nx,zx,dzx
     +      ,(i,(PZX(j,i),j=1,2),(BZX(j,i),j=1,2),i=1,NPNT)
 1000 FORMAT(' Number of points=',I7
     +      /' nz,nx,zmin,zmax,xmin,xmax=',2I5,3X,4F8.2,3X,2F7.2
     +      /(I7,2F8.2,3X,2F11.5))
C
 999  CONTINUE
      END



