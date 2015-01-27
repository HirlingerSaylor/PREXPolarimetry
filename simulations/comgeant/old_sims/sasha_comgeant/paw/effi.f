      SUBROUTINE EFFI(ID1,ID2,ID3)
      INTEGER ID1,ID2,ID3
C
      DIMENSION EF(5000),EREF(5000)
      CHARACTER*80 TIT,tite
      LOGICAL HEXIST
C      WRITE(6,*) 'start EFFI: ID1,ID2,ID3 ', ID1,ID2,ID3
C
      WRITE(6,*) id1,id2,id3
      IF(.NOT.HEXIST(ID1)) THEN
        WRITE(6,*) ' *** No histogram ID=',ID1
        GO TO 999
      ENDIF
      IF(.NOT.HEXIST(ID2)) THEN
        WRITE(6,*) ' *** No histogram ID=',ID2
        GO TO 999
      ENDIF
      IF(.NOT.HEXIST(ID3)) THEN
        WRITE(6,*) ' *** No histogram ID=',ID3
        GO TO 999
      ENDIF
C
      CALL HGIVE(ID1,TIT,NX1,X11,X12,NY1,Y1,Y2,NWT,LLOC)
      CALL HGIVE(ID2,TIT,NX2,X21,X22,NY2,Y1,Y2,NWT,LLOC)
      IF(NY1.GT.0.OR.NY2.GT.0) THEN
        WRITE(6,*) ' *** Error: only 1-dim. histograms are allowed'
        GO TO 999
      ENDIF
      IF(NX1.NE.NX2.OR.X11.NE.X21.OR.X12.NE.X22) THEN
        WRITE(6,*) ' *** Error: not identical input histograms'
        GO TO 999
      ENDIF
C
      CALL HBARX(ID3)
C
      WRITE(6,*) id1,id2,id3,nx1
      DO i=1,NX1
        a1=HI(ID1,i)
        a2=HI(ID2,i)
        EF(i)=0.
        EREF(i)=0.
          WRITE(6,*) i,a1,a2,EF(i)
        IF(a1.GT.0..AND.a2.GT.0.) THEN
          EF(i)=a2/a1
          EREF(i)=SQRT(EF(i)*(1.-EF(i))/a1)
        ENDIF
      END DO
C
C      WRITE(6,*) 'EFFI: ID1,ID2,ID3 ', ID1,ID2,ID3
C
      CALL HPAK(ID3,EF)
      CALL HPAKE(ID3,EREF)

C      WRITE(6,*) 'EFFI: END'
C
  999 RETURN
      END
