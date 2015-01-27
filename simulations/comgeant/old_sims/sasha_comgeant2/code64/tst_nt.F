      PROGRAM TST_NTUP
C
C===  Reads ntuple
C
      IMPLICIT NONE
C
      INTEGER    mxlenb,mxeve,mxhbook
      PARAMETER (mxlenb=8192,mxeve=10000000,mxhbook=2000000)
C
      INTEGER mhbook
      COMMON/PAWC/ mhbook(mxhbook)
C
#include "omgeant/omcntup.inc"
#include "omgeant/omcnthea.inc"
#include "omgeant/bg_evec.inc"
#include "omgeant/cquest.inc"
C
      CHARACTER fnam*64,cht*4,chtit*80,ctyp(2)*6,cdir*32
      INTEGER lrec,istat,lenp,id,iev,ierr,idh,io,lout,icycle
      DATA ctyp/'output','input '/
C
      CALL HLIMIT(mxhbook)
C
      lout=6
      NTUPNAM(1)='OMNTUPO'
      NTUPNAM(2)='OMNTUPI'
      LUNNTUP(1)=21
      LUNNTUP(2)=22
      IQUEST(10)= 128000        !  max number of records for NTUPLE
C
      fnam='csimout.nt'
      lrec=1024
      io=1
      CALL HROPEN(LUNNTUP(IO),NTUPNAM(IO),fnam,'NQ',lrec,istat)
      CALL HCDIR('//'//NTUPNAM(IO),' ')
C
C---          OUTPUT
C
      id=2
      CALL HBNT(id,'OMGEANT every event',' ')
      CALL HBNAME(id,'run',IEVE   ,'ieve')
      CALL HBNAME(id,'run',IRUN   ,'irun')
      CALL HBNAME(id,'run',IEND   ,'iend[-120,120]')
C
      lenp=100
      id=2
      write(6,*) ' output call ..xini'
      CALL OMNTXINI(0,lenp,6,id)
      CALL HLDIR(' ','R')
C      
      
      id=9
C      CALL HCDIR('//PAWC',' ')
C      CALL HLDIR(' ',' ')
C      CALL HDELET(id)
C
      io=2
      fnam='csimin.nt'
      lrec=0
      CALL HROPEN(LUNNTUP(IO),NTUPNAM(io),fnam,' ',lrec,istat)
      WRITE(6,*) 'istat=',istat
      WRITE(LOUT,1300) ctyp(IO),IO,NTUPNAM(IO),LUNNTUP(IO),fnam,lrec
        CALL HLDIR(' ',' ')
C        CALL HLDIR(' ','A')
C        CALL HLDIR(' ','I')
C        CALL HLDIR(' ','N')
C        CALL HLDIR(' ','R')
C        CALL HLDIR(' ','T')
        WRITE(6,1120) 
        idh=0
 10     CONTINUE
        cht='N'
        CALL HLNEXT(idh,cht,chtit,' ')    ! scan the contents of the NTUPLE
        IF(idh.NE.0) THEN
          WRITE(6,1130) idh,cht,chtit
          IF(NIDNTUP(IO).LT.MXIDNT) THEN
            NIDNTUP(IO)=NIDNTUP(IO)+1
            IDNTUP(NIDNTUP(IO),IO)=idh
          ELSE
            WRITE(6,*) ' *** Error: too many ID-s'
          ENDIF
          GO TO 10
        ENDIF
C

C      CALL HCDIR('//'//NTUPNAM(io),' ')
      CALL HCDIR(cdir,'R')
      write(6,*) ' Directory=',cdir
C      CALL HPRNT(1)
C      CALL HPRNT(3)
      CALL HRIN(id,9999,0)
      CALL HPRNT(id)
      CALL HBNAME(id,' ',0,'$CLEAR')
C
      CALL HBNAME(id,'run',IEVE   ,'$SET:ieve')
      CALL HBNAME(id,'run',IRUN   ,'$SET:irun')
      CALL HBNAME(id,'run',IEND   ,'$SET:iend')
C
      lenp=-1
C      id=2
      write(6,*) ' call ..xini'
      CALL OMNTXINI(0,lenp,6,id)
C
      DO iev=1,10000
C
         IEVE=0
         IRUN=0
         IEND=0
         IEVPROC=-1
C
C         write(6,*) ' try to read iev=',iev
         CALL OMNTUPRD(iev,ierr)

         id=2
C         CALL HCDIR('//PAWC',' ')
C         CALL HDELET(id)
         CALL HCDIR('//'//NTUPNAM(1),' ')
C         CALL HCDIR(cdir,'R')
C         write(6,*) ' Directory=',cdir
         IEVPROC=IEVPROC+100
         CALL HFNT(id)
      ENDDO
C
      id=2
      CALL HCDIR('//'//NTUPNAM(1),' ')
      CALL HROUT(id,icycle,' ')
      CALL HREND(NTUPNAM(1))
C
 1120 FORMAT(5X,'===>Directory:'
     +     /20X,'ID',5X,'type',5X,'name')
 1130 FORMAT(10X,I12,5X,A4,5X,A80)
 1300 FORMAT(' --- OMNTUPIN: ',A6,' ntuple ',I3,1X,A8,'  ==> LUN=',I3
     +      ,'  file=',A16,'  lrec=',I8)
      END
      

*-- Author :    E.Chudakov   15/02/2008
C
      SUBROUTINE OMNTXINI(NTFLA,LENMX,LOUT,ID)
C
C     ******************************************************************
C     *                                                                *
C     *   Books a CWN ntuple for the BGGEN parameters                  *
C     *   INPUT :  NTFLA - flag                                        *
C     *            LENMX - max length                                  *
C     *            LENMX >=0 - output NTUPLE                           *
C     *                  =-1 - input  NTUPLE                           *
C     *            LOUT  - logical unit for printouts                  *
C     *   OUTPUT:  ID - id of the ntuple                               *
C     *                                                                *
C     *    ==>Called by : OMNTUPIN                                     *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER  NTFLA,LENMX,LOUT,ID
#include "omgeant/omcfla.inc"
#include "omgeant/bg_evec.inc"
CC
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id: $'/
#endif
C
      INTEGER   nmxt
      CHARACTER lenv*3,lent*3
C
C     ------------------------------------------------------------------
C
C      ID=0
C
      IEVPROC=-1
      NTRABG=0
C
      IF(LENMX.GT.2) THEN
C
C---       INPUT ntuple
C
        nmxt=MXTRABG
C
        WRITE(lent,FMT='(I3)') nmxt
C
C        ID=1
C
        CALL HBNAME(ID,'bgki',IEVPROC,'iproc')
        CALL HBNAME(ID,'bgki',ITPIN(1,1) ,'itypin(2,2)')
        CALL HBNAME(ID,'bgki',AMIN(1) ,'amin(2)')
        CALL HBNAME(ID,'bgki',PIN(1,1),'pin(3,2)')
        CALL HBNAME(ID,'bgki',NTRABG,'np[0,'//lent//']')
        CALL HBNAME(ID,'bgki',ITPTRA(1,1) ,'ityp(6,np)')
        CALL HBNAME(ID,'bgki',AMTRA(1) ,'am(np)')
        CALL HBNAME(ID,'bgki',PTRA(1,1),'pout(3,np)')
C
      ELSE IF(LENMX.EQ.-1) THEN
C
C---       INPUT ntuple
C
C        ID=1
C
        CALL HBNAME(ID,'bgki',IEVPROC    ,'$SET:iproc')
        CALL HBNAME(ID,'bgki',ITPIN(1,1) ,'$SET:itypin')
        CALL HBNAME(ID,'bgki',AMIN(1)    ,'$SET:amin')
        CALL HBNAME(ID,'bgki',PIN(1,1)   ,'$SET:pin')
        CALL HBNAME(ID,'bgki',NTRABG     ,'$SET:np')
        CALL HBNAME(ID,'bgki',ITPTRA(1,1),'$SET:ityp')
        CALL HBNAME(ID,'bgki',AMTRA(1)   ,'$SET:am')
        CALL HBNAME(ID,'bgki',PTRA(1,1)  ,'$SET:pout')
C
      ELSE
        WRITE(LOUT,1000) LENMX
        GO TO 999
      ENDIF
C
 999  RETURN
 1000 FORMAT(' *** OMNTXINI error: the length specified in NTUPLEN ='
     +      ,I5,' is too short')
      END
*-- Author :    E.Chudakov   15/02/2008
C
      SUBROUTINE OMNTUPRD(IEV,IRET)
C
C     ******************************************************************
C     *                                                                *
C     *   Reads a CWN ntuple (BGGEN)                                   *
C     *   INPUT :  IEV - event number                                  *
C     *   OUTPUT:  IRET=0 - OK,                                        *
C     *                =1 - IEV > number of enets in NTUPLE            *
C     *                =2 - error                                      *
C     *                                                                *
C     *    ==>Called by : GUKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER  IEV,IRET
#include "geant321/gcunit.inc"
#include "omgeant/omcfla.inc"
#include "omgeant/omcntup.inc"
#include "omgeant/bg_evec.inc"
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id: $'/
#endif
C
      INTEGER   id,io,nent,ierr
      CHARACTER lenv*3,lent*3
C
C     ------------------------------------------------------------------
C
      LOUT=6

      io=2
C
      CALL HCDIR('//'//NTUPNAM(io),' ')
      id=9
      CALL HNOENT(id,nent)
C      write(6,*) ' nent=',nent
C      IF(nent.EQ.0) THEN
C         id=1
C         CALL HNOENT(id,nent)
C      ENDIF
C
      IF(IEV.GT.nent) THEN
         WRITE(LOUT,1000) IEV,nent
         IRET=1
         GO TO 999
      ENDIF
C
      CALL HGNT(id,IEV,ierr)
C      write(6,*) ' ierr,IEVPROC...',ierr,IEVPROC,NTRABG
      IF(ierr.NE.0) THEN
         IRET=2
         WRITE(LOUT,1100) IEV
         GO TO 999
      ENDIF
C
      IRET=0
C
 999  RETURN
 1000 FORMAT(' *** OMNTUPRD error: tries to read event ',I8
     +      ,' while the ntuple contains only ',I8)
 1100 FORMAT(' *** OMNTUPRD error: reading event ',I8)
      END
