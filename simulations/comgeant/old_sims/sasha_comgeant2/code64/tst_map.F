      PROGRAM TST_MAP
      IMPLICIT NONE
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcunit.inc"
      REAL b(3),x(3),rr,phideg,phi,pi,bpol(3)
      INTEGER mtyp,ierr
      LOUT=6
      ISWIT(9)=0
      LUNTMP=9
      pi=2.*ACOS(0.)
C
 10   WRITE(6,*) 'Enter R, phi(deg), Z '
      READ(5,*) rr,phideg,x(3)
C
      phi=phideg*pi/180.
      x(1)=rr*COS(phi)
      x(2)=rr*SIN(phi)
C
      mtyp=1
      CALL OMAGMAP(1,mtyp,x,b,ierr)
      bpol(3)=b(3)
      bpol(1)= b(1)*COS(phi)+b(2)*SIN(phi)
      bpol(2)=-b(1)*SIN(phi)+b(2)*COS(phi)
      WRITE(6,FMT='(3F7.1,4X,3E14.5,4X,2F7.1,4X,3F14.5)') 
     +            x,b,rr,phideg,bpol
      GO TO 10
      END
C
* $Header:$
* $Log:$
*
      SUBROUTINE OMAGMAP(MAG,MTYP,XX,BX,IERR)
C
C     *****************************************************************
C     *                                                               *
C     *   Returns the mag. field vector for the point X,              *
C     *    For the magnet #MAG                                        *
C     *    If MAG is not defined: tries to read in the map,           *
C     *    accordingly to the map type MTYP                           *
C     *       MTYP=0 - SNAKE type maps (from formatted files)         *
C     *            1 - polar coordinates (from unformatted files)     *
C     *   INPUT: MAG - the magnet number, map file: magnet_MAG.map    * 
C     *          XX  - the coordinates                                *
C     *   OUTPUT: BX(1-3)  the components of the field                *
C     *           IERR>0 - error                                      *
C     *                                                               *
C     *   Only a part of the field is often set. For the rest of      *
C     *  the volume an extrapolation done, for various cases:         *
C     *   An extrapolation flag is assigned to each of 6 boundaries   *
C     *   kextrap(2,3,map). Let us consider the boundary (j,i):       *
C     *     =0 - no extrapolation, the end                            *
C     *     =1 - mirror across the boundary:                          *
C     *                            x(i)->2*xlim(j,i,map)-x(i)         *
C     *                     the other coordinates are unchanged       *
C     *     =2 - periodic structure, there are cases:                 *
C     *         other side =2 - the period is    xlim(2..)-xlim(1..)  *
C     *         other side =1 - the period is 2*(xlim(2..)-xlim(1..)) *
C     *               x(i)->MOD(x(i)-xlim(1..),period)                *
C     *                                                               * 
C     *   Whenever extrapolation is done across a mirror boundary,    *
C     *         the B directions may change:                          * 
C     *     possible permutations are defined in array SYMMF(3,isym)  * 
C     *      as B(i1)=B(i1)*SYMMF(j,isym), where j=1,3, i1=i+j-1, cycle
C     *      isym          SYMMF                                      * 
C     *       1          -1.  1.  1.     (asymmetrical)               * 
C     *       2           1. -1. -1.      (symmetrical)               * 
C     *       3           1. -1.  1.                                  * 
C     *       4           1.  1. -1.                                  * 
C     *                                                               * 
C     *   Read the #MAG map from SNAKE, or similar programs           * 
C     *   SNAKE frame: XY-bend frame, Y - along the beam (left handed)*
C     *   We convert it to: XYZ(SNAKE)==>XZY(ours) (right handed)     *
C     *   A SNAKE-emulation format can be read (EMUL key) -           *
C     *            in this case no conversion is used                 *
C     *                                                               *
C     *   file FORMAT:                                                *
C     *     1-st line: if contains "EMUL" - SNAKE emulation, normal   *
C     *                                     coordinate system         *
C     *                otherwise - SNAKE map (needs Y->Z, Z->Y)       *
C     *     2-nd line: Symmetry flag - 3 characters, 1 character per  *
C     *         projection, say "NSA". This means:                    *
C     *         At reflection of coordinate                           *
C     *         X=-X - "N" - no symmetry defined                      *
C     *         Y=-Y - "S" -  "symmetric"  By= By, Bz=-Bz, Bx=-Bx     *
C     *         Z=-Z - "A" - "asymmetric"  Bz=-Bz, Bz= Bz, Bx= Bx     *
C     *              other keys: Q,R, example QRN:                    *
C     *         X=-X - "Q" -               Bx= Bx, By=-By, Bz= Bz     *
C     *         Y=-Y - "R" -               By= By, Bz= Bz, Bx=-Bx     *
C     *     3-rd line: X: number of points and range                  *
C     *     lines 4-5 - the same for Y,Z                              *
C     *     line 6 - format                                           *
C     *     line 7-... - data                                         *
C     * Example:                                                      *
C     *  SNAKE EMULATION of Big Bite Magnet - from MAFIA
C     *SNN
C     *   8     -10.000   130.000
C     *  76    -750.000   750.000
C     *  86    -850.000   850.000
C     * (3G18.8,3X,3G18.8)              
C     *    0.40961332E-02   -0.38677853E-03   -0.45087747E-03       0.40961332E-02    0.38677853E-03    0.45087747E-03
C     * ... etc                                                       *
C     *                                                               *
C     *   called from: OMAGFLDN                                       *
C     *                                                               *
C     *****************************************************************
C
      IMPLICIT NONE
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
      INTEGER MAG,MTYP,IERR
      REAL XX(3),BX(3)
C
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "omgeant/omcunit.inc"
#include "omgeant/cmagmap.inc"
C
      INTEGER map,maga,iost,nline,i,j,m,lcw,ip,jp,ixyz(3),ipri,ix,iy,iz
     +       ,ipz,ipy,n12,isym(2),k,isp,nper
      CHARACTER fnam*16,cnm*2
C
      REAL    bmax,bb,xyz(3),qq,sxyz(2,3),btmp(3),s3,s32,s321,a1,a2
     +       ,refl(3),rr,phi,degrad,pi,xyz1(3),dsiz,xlima(2,3),dbin(3)
      REAL symm(3,MXSYM)      ! (j,i) symmetry in projection Bj for case i, i=1 - "asymm", =2 -"symm" etc
C
      INTEGER nhead,icoilt,ncoils  ! parameters of the toroid's coils
      REAL    curwire
C
      INTEGER ifirst
      DATA ifirst/1/
      DATA ipri/0/
C   
C---    "Symmetry" flags for reflection of 3 map coordinates: (ex: NAS, N - no symmetry,                                   
C    X(j) --> -X(j), where j=1,2,3  B(j)     B(j+1)    B(j+2)     B(j) - j-th projection of B, if j+1=4 then j=1 ...(a cyclic index)
      DATA symm/-1., 1., 1.     !  -B(j)     B(j+1)    B(j+2)  - "asymmetric", key "A"
     +         , 1.,-1.,-1.     !   B(j)    -B(j+1)   -B(j+2)  - "symmetric",  key "S"
     +         , 1.,-1., 1.     !   B(j)    -B(j+1)    B(j+2)  -               key "Q"
     +         , 1., 1.,-1.     !   B(j)     B(j+1)   -B(j+2)  -               key "R"
     +         / 
C
C     ------------------------------------------------------------------
C
C---    First call - read the map
C
      ipri=0
      IF(ISWIT(9).GE.10) ipri=ISWIT(9)
      IF(ipri.EQ.20) WRITE(LOUT,*) ' OMAGMAP for MAG=',MAG
     +       ,' called for ',XX
      BX(1)=0.
      BX(2)=0.
      BX(3)=0.
      IERR=0
C
C---     First call - initialize the COMMON
C
      IF(ifirst.EQ.1) THEN
         ifirst=0
         NMAP=0
         DO maga=1,MXMAG
            MAPMAG(maga)=0  ! map counter
         ENDDO
         DO j=1,MXSYM
            DO i=1,3
               SYMMF(i,j)=symm(i,j)
            ENDDO
         ENDDO
      ENDIF
C
      IF(NMAP.LT.0) THEN
         IERR=1
         GO TO 999
      ENDIF
C
      IF(MAG.LE.0.OR.MAG.GT.MXMAG) THEN
         WRITE(LOUT,1500) MAG
         GO TO 998
      ENDIF
C
C---    Try to find the map
C
      IF(MAPMAG(MAG).EQ.0) THEN
C
C---      Try to read this map
C
         IF(NMAP.GE.MXMAP) THEN
            WRITE(LOUT,1600) MAG,NMAP
            GO TO 998
         ENDIF
C
         IF(MAG.LT.10) THEN
            WRITE(cnm,FMT='(I1)') MAG
            lcw=1
         ELSE
            WRITE(cnm,FMT='(I2)') MAG
            lcw=2
         ENDIF
         fnam='magnet_'//cnm(1:lcw)//'.map'
         lcw=INDEX(fnam,'map')+2
C
         IF(MTYP.EQ.0) THEN  ! SNAKE
C
            OPEN(LUNTMP,FILE=fnam(1:lcw),STATUS='OLD',IOSTAT=iost
     +           ,FORM='FORMATTED')
            IF(iost.NE.0) THEN
               WRITE(LOUT,2000) MAG,fnam(1:lcw),iost
               GO TO 998
            ENDIF
C
C     
         ELSE IF(MTYP.EQ.1) THEN  ! r,phi,z  periodic in phi
C
            OPEN(LUNTMP,FILE=fnam(1:lcw),STATUS='OLD',IOSTAT=iost
     +           ,FORM='UNFORMATTED')
            IF(iost.NE.0) THEN
               WRITE(LOUT,2000) MAG,fnam(1:lcw),iost
               GO TO 998
            ENDIF
C     
            NMAP=NMAP+1
            MAGMAP(NMAP)=MAG
            MAPMAG(MAG)=NMAP
            ISPNT(NMAP)=0
            IF(NMAP.GT.1) ISPNT(NMAP)=ISPNT(NMAP-1)+NPNT(NMAP-1) ! address of the map
C
            nline=1
            READ(LUNTMP,IOSTAT=iost) nhead,(NXYZ(i,NMAP),i=1,3)
     +                                    ,icoilt,ncoils
            IF(iost.NE.0) GO TO 997
            nline=nline+1
            READ(LUNTMP,IOSTAT=iost) nhead
     +                              ,((XLIM(j,i,NMAP),j=1,2),i=1,3)
     +                              ,curwire
            IF(iost.NE.0) GO TO 997
C
            NPNT(NMAP)=NXYZ(1,NMAP)*NXYZ(2,NMAP)*NXYZ(3,NMAP)
            IF(NPNT(NMAP).LT.1.OR.
     +         XLIM(1,1,NMAP).GE.XLIM(2,1,NMAP).OR.
     +         XLIM(1,2,NMAP).GE.XLIM(2,2,NMAP).OR.
     +         XLIM(1,3,NMAP).GE.XLIM(2,3,NMAP)) THEN
               WRITE(LOUT,2010) MAG,NMAP
     +                 ,(i,NXYZ(i,NMAP),(XLIM(j,i,NMAP),j=1,2),i=1,3)
               GO TO 998
            ENDIF
            IF(NPNT(NMAP)+ISPNT(NMAP).GT.MXPNT) THEN
               WRITE(LOUT,2020) MAG,NMAP,ISPNT(NMAP),NPNT(NMAP)
               GO TO 998
            ENDIF
C
            k=0
            DO j=1,NXYZ(3,NMAP)*NXYZ(2,NMAP)
               READ(UNIT=LUNTMP,IOSTAT=iost)
     +             ((BXYZ(i,ISPNT(NMAP)+k+ix),i=1,3),ix=1,NXYZ(1,NMAP))
               IF(iost.NE.0) THEN
                  WRITE(LOUT,2030) MAG,NMAP,j
                  GO TO 998
               ENDIF
               k=k+NXYZ(1,NMAP)
            ENDDO
C
            MAPAXIS(1,NMAP)=4  ! R
            MAPAXIS(2,NMAP)=5  ! phi
            MAPAXIS(3,NMAP)=3  ! Z
C
            DO i=1,3
               DO j=1,2
                  IF(MAPAXIS(i,NMAP).EQ.5) THEN
                     MAPLIMX(j,i,NMAP)=2  ! symmetric extrapolations (B_phi - the same, B_R,B_Z - reflected)
                  ELSE
                     MAPLIMX(j,i,NMAP)=0
                  ENDIF
               ENDDO
            ENDDO
C
         ENDIF
C
         MAPTYP(NMAP)=MTYP
C
         CLOSE(LUNTMP)
C
         bmax=0.
         DO ip=1,NPNT(NMAP)
            isp=ISPNT(NMAP)+ip
            bb=0.
            DO i=1,3
               BXYZ(i,isp)=BXYZ(i,isp)*10.   !  T --> kGs
               bb=bb+BXYZ(i,isp)**2
            ENDDO
            bmax=MAX(bmax,SQRT(bb))
         ENDDO
C
         DO i=1,3
            DXYZ(i,NMAP)=(XLIM(2,i,NMAP)-XLIM(1,i,NMAP))
     +                  /(NXYZ(i,NMAP)-1)
         ENDDO
C
         WRITE(LOUT,1100) MAG,NMAP,(NXYZ(i,NMAP),i=1,3),NPNT(NMAP)
     +                 ,((XLIM(j,i,NMAP),j=1,2),i=1,3)
     +                 ,(DXYZ(i,NMAP),i=1,3),bmax
 1100    FORMAT(' === Read the field map of mag=',I2,' into map=',I2
     +         /6X,'nx,ny,nz,npnt=',3I5,I8,1X
     +         /,' x,y,z limits=',3(2F8.2,1X)
     +        /,' x,y,z steps=' ,3F7.2,2X,'Bmax=',F10.3)
C
         IF(ipri.EQ.20) THEN
            DO ip=1,MIN(NPNT(NMAP),1000)
               k=ip
               DO i=3,1,-1
                  ixyz(i)=MOD(k-1,NXYZ(i,NMAP))+1
                  k=(k-1)/NXYZ(i,NMAP)
                  xyz(i)=XLIM(1,i,NMAP)+(ixyz(i)-1)*DXYZ(i,NMAP)
               ENDDO
               WRITE(LOUT,FMT='(I8,2X,3I4,3X,3F8.2,3X,3F10.3)') 
     +              ip,(ixyz(i),i=1,3)
     +                ,(xyz(i),i=1,3)
     +                ,(BXYZ(i,ip+ISPNT(NMAP)),i=1,3)
            ENDDO
         ENDIF
      ENDIF
C
      map=MAPMAG(MAG)
      IF(map.LT.0.OR.map.GT.MXMAP) THEN
         WRITE(LOUT,1900) MAG,map
         GO TO 998
      ENDIF
C
C---    
C
      pi=ACOS(0.)*2.
      rr=SQRT(XX(1)**2+XX(2)**2)
      phi=ATAN2(XX(2),XX(1))
      IF(phi.LT.0.) phi=phi+2.*pi
      degrad=pi/180.
C
      DO i=1,3
         refl(i)=1.
         dbin(i)=DXYZ(i,map)
         IF(MAPAXIS(i,map).EQ.5) THEN
            dbin(i)=dbin(i)*degrad
         ENDIF
         DO j=1,2
            xlima(j,i)=XLIM(j,i,map)
            IF(MAPAXIS(i,map).EQ.5) THEN
               xlima(j,i)=xlima(j,i)*degrad
            ENDIF
         ENDDO
      ENDDO
C
C---     Find the location of the point on the grid
C
      DO i=1,3
         dsiz=xlima(2,i)-xlima(1,i)
         IF(MAPAXIS(i,map).GE.1.AND.MAPAXIS(i,map).LE.3) THEN
            xyz(i)=XX(i)
         ELSE IF(MAPAXIS(i,map).EQ.4) THEN
            xyz(i)=rr
         ELSE IF(MAPAXIS(i,map).EQ.5) THEN
            xyz(i)=phi
         ELSE
            WRITE(LOUT,2300) MAG,map,i,MAPAXIS(i,map)
            GO TO 998
         ENDIF
C
         xyz1(i)=xyz(i)
         DO j=1,2
            isym(j)=MAPLIMX(j,i,map)
            IF(isym(j).LT.0.OR.isym(j).GT.MXSYM) THEN  ! no symmetry defined
               WRITE(LOUT,1950) MAG,map,i,j,isym(j)
               GO TO 998
            ENDIF
         ENDDO
C
         IF(xyz(i).LT.xlima(1,i)) THEN
            k=1                 ! below the range
         ELSE IF(xyz(i).GT.xlima(2,i)) THEN
            k=2                 ! above the range
         ELSE
            k=0                 ! within the range
         ENDIF
C
C         write(6,FMT='(4I4,3F10.3)') i,k,isym,XX(i),xyz(i),xyz1(i)

         IF(k.NE.0) THEN
            IF(isym(k).EQ.0) GO TO 999 ! outside of the volume
            nper=INT((xyz(i)-xlima(3-k,i))/(2.*dsiz))      ! extended period=2*period 
            IF(isym(3-k).EQ.0.AND.ABS(nper).GT.0) GO TO 999   ! no extrapolation on the other side - no period 
            xyz1(i)=xyz(i)-2.*dsiz*nper                       ! translation to the main double range
C            write(6,*) 'nper=',nper,xyz1(i)
            IF( xyz1(i).LT.xlima(1,i).OR.
     +          xyz1(i).GT.xlima(2,i)) THEN
               xyz1(i)=2.*xlima(k,i)-xyz1(i)                  ! reflection
C
               DO j=1,3
                  m=i+j-1
                  IF(m.GT.3) m=m-3
C                  write(6,*) 'refl',j,m,k,isym(k),refl(m)
C     +                           ,SYMMF(j,isym(k))
                  refl(m)=refl(m)*SYMMF(j,isym(k))
               ENDDO
C                  
            ENDIF
         ENDIF
C
         qq=xyz1(i)-xlima(1,i)
         ixyz(i)=INT(qq/dbin(i))+1
         sxyz(2,i)=(qq-(ixyz(i)-1)*dbin(i))/dbin(i)
         sxyz(1,i)=1.-sxyz(2,i)
         btmp(i)=0.
C
      ENDDO
C
C---       Take an average of 8 points using weights depending on the distances. 
C---           On a point the result coincides with the value at this point
C
C      WRITE(LOUT,FMT='(3I4,2X,6F6.3)') ixyz,sxyz 
      n12=NXYZ(1,map)*NXYZ(2,map)
      DO iz=0,1
         ipz=(ixyz(3)+iz-1)*n12+ISPNT(map)
         s3=sxyz(iz+1,3)
         DO iy=0,1
            ipy=ipz+(ixyz(2)+iy-1)*NXYZ(1,map)
            s32=sxyz(iy+1,2)*s3
            DO ix=0,1
               ip=ipy+ixyz(1)+ix
               s321=sxyz(ix+1,1)*s32
C               WRITE(LOUT,FMT='(3I4,2X,3F9.2,2X,3F9.3)') 
               DO i=1,3
                  btmp(i)=btmp(i)+BXYZ(i,ip)*s321*refl(i)
               ENDDO
               IF(ipri.EQ.12) WRITE(LOUT,
     +          FMT='(3I3,2X,F8.3,I8,3(2X,2F9.3),2X,3F9.3)')
     +          ix,iy,iz,s321,ip,(refl(j),BXYZ(j,ip),j=1,3),btmp
            ENDDO
         ENDDO
      ENDDO
C
      IF(MAPTYP(map).EQ.0) THEN
         DO i=1,3
            BX(i)=btmp(i)
         ENDDO
      ELSE IF(MAPTYP(map).EQ.1) THEN
         BX(1)=btmp(1)*COS(phi)-btmp(2)*SIN(phi)
         BX(2)=btmp(1)*SIN(phi)+btmp(2)*COS(phi)
         BX(3)=btmp(3)
      ENDIF
C     
      RETURN
C
 997  IERR=1
      IF(iost.LT.0) WRITE(LOUT,2100) MAG,nline
      IF(iost.GT.0) WRITE(LOUT,2200) MAG,nline
C
 998  IERR=1
C
 999  RETURN
 1500 FORMAT(' *** Error in OMAGMAP, magnet=',I4,' is out of range')
 1600 FORMAT(' *** Error in OMAGMAP, magnet=',I2
     +      ,' no space for a new map, already define ',I3)
 2010 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +      ,' wrong steps/limits ',/(3X,I2,I6,2E13.4))
 2020 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +      ,' no space for this map ',2I8)
 2030 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +      ,' error reading the column iz*iy=',I8)
 1900 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +     ,' map reference is out of range ')
 1950 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +     ,' invalid symmetry ',3I6)
 2000 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' no file opened',A16
     +      ,'  iost=',I9)
 2100 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' unexpected EOF'
     +      ,'  line=',I8)
 2200 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' read error '
     +      ,',  line=',I8)
 2300 FORMAT(' *** Error in OMAGMAP, magnet=',I2,' map=',I2
     +     ,' MAPAXIS not defined ',2I6)
      END
