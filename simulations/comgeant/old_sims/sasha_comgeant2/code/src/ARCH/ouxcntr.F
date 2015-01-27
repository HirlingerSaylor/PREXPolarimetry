* $Header:$
* $Log:$
*
      SUBROUTINE OUXCNTR
C
C     ******************************************************************
C     *                                                                *
C     *       USER's routine                                           *
C     *    KUIP control commands                                       *
C     *                                                                *
C     *    ==>Called by :                                              *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "omgeant/omcsvol.inc"
#include "omgeant/omcbeam.inc"
#include "omgeant/omcunit.inc"
#include "omgeant/cescatt.inc"
#include "moller/comnfla.inc"
#include "omgeant/ranma2.inc"
#include "moller/csnaprnd.inc"
C
      INTEGER npar,i,n,numed1,lenc,iseq,iev,iost,ipr,lrec
      REAL    ecut,a,x1,x2,xyz(3),bxyz(3)
      CHARACTER chpatl*32,cnam*80
C
C     ------------------------------------------------------------------
C
      CALL KUPATL(chpatl,npar)
C      write(lout,*) npar,'   ',chpatl
C
C      write(6,*) ' ** command,npar: ',chpatl,npar
C
      IF(chpatl(1:8).EQ.'MAGFIELD'.AND.npar.GE.2) THEN
        CALL KUGETI(n)
        IF(n.GE.1.AND.n.LE.MXMAG) THEN
           DO i=1,npar-1
              CALL KUGETR(a)
              IF(ABS(a+9999.).GT.0.5) THEN
                 AMAGPAR(9+i,n)=a
              ENDIF
           ENDDO
        ELSE
          write(6,*) ' *** Error: IPAR=',i,' is out of range'
        ENDIF
C
      ELSE IF(chpatl(1:8).EQ.'MOLLANGL'.AND.npar.EQ.4) THEN
         CALL KUGETR(THETMOL(1))
         CALL KUGETR(THETMOL(2))
         CALL KUGETR(PHIMOL(1))
         CALL KUGETR(PHIMOL(2))
         IF(IDMOLTET.NE.0) THEN
            WRITE(LOUT,1000) THETMOL,PHIMOL
 1000       FORMAT(' == Change the Moller scattering angles: '
     +            ,(2X,2F7.2)) 
            CALL OMOLINI
         ELSE
            WRITE(LOUT,*) '*** Error: Moller scattering is not set up'
         ENDIF
C
      ELSE IF(chpatl(1:8).EQ.'MOTTANGL'.AND.npar.EQ.4) THEN
         CALL KUGETR(THETMOT(1))
         CALL KUGETR(THETMOT(2))
         CALL KUGETR(PHIMOT(1))
         CALL KUGETR(PHIMOT(2))
         IF(IDMOTT.NE.0) THEN
            WRITE(LOUT,1050) THETMOT,PHIMOT
 1050       FORMAT(' == Change the Mott scattering angles: '
     +            ,(2X,2F7.2)) 
            CALL OMOTTINI
         ELSE
            WRITE(LOUT,*) '*** Error: Mott scattering is not set up'
         ENDIF
C
      ELSE IF(chpatl(1:7).EQ.'MAGSIZE'.AND.npar.GT.1) THEN
         CALL KUGETI(n)
         IF(n.GT.0.AND.n.LE.MXMAG) THEN
            DO i=2,npar
               CALL KUGETR(a)
               IF(a.GT.0.0001) AMAGPAR(6+i,n)=a
            END DO
         ELSE
            WRITE(LOUT,*) ' OUXCNTR error: mag=',n,' is out of range' 
         ENDIF
C
      ELSE IF(chpatl(1:7).EQ.'BEAMMOM'.AND.npar.EQ.3) THEN
         CALL KUGETR(BEAMOLIM(1))
         CALL KUGETR(BEAMOLIM(2))
         CALL KUGETR(BEAMDSP)
C
      ELSE IF(chpatl(1:7).EQ.'CHERFAC'.AND.npar.EQ.1) THEN
         CALL KUGETR(CHERFAC)
      ELSE IF(chpatl(1:6).EQ.'PRIMAG'.AND.npar.EQ.3) THEN
         CALL KUGETR(xyz(1))
         CALL KUGETR(xyz(2))
         CALL KUGETR(xyz(3))
         NLEVEL=0
         CALL GUFLD(xyz,bxyz)
         WRITE(LOUT,FMT='(3X,3F9.2,4X,3F9.3)') xyz,bxyz
C
      ELSE IF(chpatl(1:6).EQ.'PRIMED'.AND.npar.EQ.3) THEN
         CALL KUGETR(xyz(1))
         CALL KUGETR(xyz(2))
         CALL KUGETR(xyz(3))
         NLEVEL=0
         CALL GMEDIA(xyz,numed1)
         WRITE(LOUT,1100) xyz,numed1
 1100    FORMAT(' == Point ',3F9.2,3X,' medium=',I6/' == Volume tree:')
         WRITE(LOUT,FMT='(6X,A4,3X,I3)') (NAMES(i),NUMBER(i),i=1,NLEVEL)
C
      ELSE IF(chpatl(1:6).EQ.'KINONE'.AND.npar.EQ.7) THEN
         CALL KUGETI(IPARONE)
         DO i=1,3
            CALL KUGETR(PAVONE(i))
         ENDDO
         DO i=1,3
            CALL KUGETR(DPAVONE(i))
         ENDDO
C
      ELSE IF(chpatl(1:7).EQ.'GETSEED'.AND.npar.EQ.3) THEN
         CALL KUGETI(iev)
         CALL KUGETI(ipr)
         CALL KUGETS(cnam,lenc)
         lrec=(103*MXRMSEQ+1)*4
         OPEN(LUNTMP,FILE=cnam(1:lenc),ACCESS='DIRECT'
     +       ,FORM='UNFORMATTED',RECL=lrec,STATUS='OLD',IOSTAT=iost)
         IF(iost.NE.0) THEN
            WRITE(LOUT,*) ' *** no file with name ',cnam,' length=',lenc
         ELSE
            READ(LUNTMP,REC=iev,IOSTAT=iost) IEVSEED,ISEEDEV
            IF(iost.NE.0) THEN
               WRITE(LOUT,*) ' *** error reading file ',cnam(1:lenc)
     +                    ,' record=',iev
            ELSE
               IF(ipr.NE.0) THEN
                  WRITE(LOUT,FMT='('' event '',I10)') IEVSEED
                  DO iseq=1,MXRMSEQ
                     WRITE(LOUT,FMT='(10I12)') (ISEEDEV(i,iseq),i=1,103)
                  ENDDO
               ENDIF
               DO iseq=1,MXRMSEQ
                  CALL RMMAQ(ISEEDEV(1,iseq),iseq,'SV')
               ENDDO
            ENDIF
            CLOSE(LUNTMP)
         ENDIF
      ENDIF
C
      RETURN
      END

