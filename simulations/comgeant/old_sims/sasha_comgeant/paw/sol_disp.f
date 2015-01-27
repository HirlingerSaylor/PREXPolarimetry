* $Header:$
* $Log:$
*
      SUBROUTINE SOL_DISP(MAG,Z,R)
C
C     *****************************************************************
C     *                                                               *
C     *   Read the map for an axially symmetric field (solenoid)      * 
C     *   And return the value of the field                           *
C     *   INPUT: MAG - the magnet number, map file: magnet_MAG.map    * 
C     *          Z,R - the coordinates                                *
C     *   OUTPUT: BZR(1-2) Z,R components of the field                *
C     *                                                               *
C     *   called from: OMAGFLDM                                       *
C     *                                                               *
C     *****************************************************************
C
      IMPLICIT NONE
C #if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
C #endif
      INTEGER MAG
      REAL Z,R
C
C #include "geant321/gcunit.inc"
C #include "omgeant/omcunit.inc"
C
C---     Read a field map
C
      INTEGER    mxpnt
      PARAMETER (mxpnt=100000)
      INTEGER    npnt
C      REAL       zrp(2,mxpnt),bzrp(2,mxpnt)
C
      VECTOR zrp(2,100000),bzrp(2,10000),BZR(2)
C
      INTEGER i,j,nz,nr,lcw,nzr(2),izr(2),ip,iz,ir
      REAL    zr0(2,2),dzr(2),bcen,dis,dismin,zr(2),ds(2,2)
      CHARACTER fnam*16,cnm*2
C
      INTEGER LUNTMP,LOUT
      SAVE npnt,zrp,bzrp,zr0,dzr,nz,nr,nzr
C
      DATA npnt/0/
C
C     ------------------------------------------------------------------
C
C---    First call - read the map
C
      LUNTMP=4
      LOUT=6
C
      BZR(1)=0.
      BZR(2)=0.
C
      IF(npnt.LT.0) GO TO 999
      IF(npnt.EQ.0) THEN
         IF(MAG.LT.10) THEN
            WRITE(cnm,FMT='(I1)') MAG
            lcw=1
         ELSE
            WRITE(cnm,FMT='(I2)') MAG
            lcw=2
         ENDIF
         fnam='magnet_'//cnm(1:lcw)//'.map'
         lcw=INDEX(fnam,'map')+2
         OPEN(LUNTMP,FILE=fnam(1:lcw),STATUS='OLD',ERR=991
     +       ,FORM='UNFORMATTED')
         READ(LUNTMP,END=992,ERR=993) npnt,nz,nr,zr0,dzr
     +        ,((zrp(j,i),j=1,2),(bzrp(j,i),j=1,2),i=1,npnt)
         CLOSE(LUNTMP)
C
C---      Find the field at the point closest to the center
C
         dismin=999999.
         DO i=1,npnt
            dis=SQRT(zrp(1,i)+zrp(2,i)**2)
            IF(dis.LT.dismin) THEN
               dismin=dis
               bcen=SQRT(bzrp(1,i)**2+bzrp(2,i)**2)
            ENDIF
         ENDDO
C
         nzr(1)=nz
         nzr(2)=nr
C
         WRITE(LOUT,1100) MAG,nz,nr,npnt,zr0,dzr,dismin,bcen
 1100    FORMAT(' === Read the field map of mag=',I2
     +         /6X,'nz,nr,npnt=',3I6,3X,' z,r limits=',4F8.2,3X
     +         ,'z,r steps',2F8.2,3X,'B(',F5.1,')=',F9.4,' kGs')
         IF(nz.LE.0.OR.nr.LE.0.OR.dzr(1).LE.0..OR.dzr(2).LE.0..OR.
     +      nz*nr.NE.npnt.OR.
     +      zr0(1,1)+dzr(1)*(nz-2).GT.zr0(1,2).OR.
     +      zr0(2,1)+dzr(2)*(nr-2).GT.zr0(2,2)) THEN
            WRITE(LOUT,1200)
 1200       FORMAT(' *** Error in the map parameters')
            npnt=-2
         ENDIF
      ENDIF
C
      zr(1)=Z
      zr(2)=R
C
C---    Take a linearly weighed average of 4 points
C
      DO i=1,2
         IF(zr(i).LT.zr0(i,1).OR.zr(i).GT.zr0(i,2)) GO TO 999
         izr(i)=(zr(i)-zr0(i,1))/dzr(i)+1
         IF(izr(i).GE.nzr(i)) GO TO 999
         ds(2,i)=(zr(i)-(zr0(i,1)+dzr(i)*(izr(i)-1)))/dzr(i)
         ds(1,i)=1.-ds(2,i)
      ENDDO
C
      DO iz=1,2
         DO ir=1,2
            ip=iz+izr(1)-1+nzr(1)*(ir-2+izr(2))
C            WRITE(6,*) ip,izr,zrp(1,ip),zrp(2,ip),bzrp(1,ip)
C     +                ,ds(iz,1),ds(ir,2)
            DO i=1,2
               BZR(i)=BZR(i)+bzrp(i,ip)*ds(iz,1)*ds(ir,2)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
C
 991  WRITE(LOUT,2000) MAG,fnam(1:lcw)
 2000 FORMAT(' *** Error in OMSOLMAP, magnet=',I2,' no file ',A16)
      npnt=-1
      RETURN
C
 992  WRITE(LOUT,2100) MAG
 2100 FORMAT(' *** Error in OMSOLMAP, magnet=',I2,' unexpected EOF')
      npnt=-1
      RETURN
C
 993  WRITE(LOUT,2200) MAG,npnt,nz,nr,zr0,dzr
 2200 FORMAT(' *** Error in OMSOLMAP, magnet=',I2,' read error'
     ,      ,3X,3I8,6F9.2)
      npnt=-1
C
 999  RETURN
      END
