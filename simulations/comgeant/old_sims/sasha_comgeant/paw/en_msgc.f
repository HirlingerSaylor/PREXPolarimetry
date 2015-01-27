      REAL FUNCTION EN_MSGC(IPAR)
C
C===    Gives the particle momentum measured by the MSGC for particle IPAR
C
      IMPLICIT NONE
      INTEGER  IPAR
      INCLUDE ?
C
      INTEGER ihi,nfind,k,itrh,idi,idg1,i1,i2,jhi,iddd,idms(2,2),jm
      REAL pp,pp1,y(2),sly(2),zl1,zl2,zl3,zz(4),del,s1,s2,c2
      INTEGER ifirst
      DATA ifirst/1/
      DATA zz/0.,422.8,556.,576./
C
C     -----------------------------------------------------------------
C
      ifirst=0
      idms(1,1)=1
      idms(2,1)=2
      idms(1,2)=3
      idms(2,2)=4
      zl1=zz(2)-zz(1)
      zl2=zz(3)-zz(2)
      zl3=zz(4)-zz(3)
C
      EN_MSGC=0.
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      nfind=0
C
      IF(npdig.GT.0) THEN
C
         DO idi=1,ndig
            k=ip1dig(idi)
            iddd=IBITS(k,0,16)
            k=ip2dig(idi)
            idg1=IBITS(k,0,16)
            jm=0
            IF(iddd.EQ.idms(1,1).OR.iddd.EQ.idms(2,1)) jm=1
            IF(iddd.EQ.idms(1,2).OR.iddd.EQ.idms(2,2)) jm=2
            IF(jm.GT.0) THEN
              
C               IF(ICHA.EQ.0.OR.idg1.EQ.ICHA) THEN
               i1=1
               IF(idi.GT.1) THEN
                  k=ip1dig(idi-1)
                  i1=IBITS(k,16,16)
               ENDIF
               k=ip1dig(idi)
               i2=IBITS(k,16,16)
C              WRITE(6,*) 'idi,idg1=',idi,idg1,jpdig(i1)
C
               DO jhi=i1,i2-1
                  ihi=jpdig(jhi)
                  IF(ihi.GT.0.AND.ihi.LE.nhit) THEN
                     k=ip1hit(ihi)
                     itrh=IBITS(k,0,16)
C               WRITE(6,*) 'iddd,idg1=',iddd,idg1,itrh,itra(itrh)
                     IF(itrh.GT.0.AND.itrh.LE.ntra) THEN
                        IF(itrh.EQ.IPAR) THEN
                           nfind=nfind+1
                           y(jm)=hit(2,ihi)
                           sly(jm)=hitsld(2,ihi)
                           pp=phit(ihi)
C               WRITE(6,*) 'nfind,jm,y(jm),pp',nfind,jm,y(jm),pp
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
C
            ENDIF

         ENDDO
C
         IF(nfind.EQ.2) THEN
            del=y(1)*(1.+zl2/zl3)-y(2)*zl2/zl3
            s1=del/SQRT(del**2+zl1**2)
            s2=(y(2)-y(1))/SQRT(zl3**2+(y(2)-y(1))**2)
            c2=SQRT(1.-s2**2)
C            write(6,*) 'del,s1,s2=',del,s1,s2
            IF(s1-s2.NE.0.) THEN
C               EN_MSGC=1.5153*0.17/(s1-s2)
C               EN_MSGC=pp
C               EN_MSGC=y(1)-(y(2)-y(1))/zl3*zl2
C     +             -zl1*ptra(2,IPAR)/ptra(3,IPAR)
               EN_MSGC=
     +             y(1)-sly(1)*(zl2-164*(1./(1+c2)-.5))
     +       -ptra(2,IPAR)/ptra(3,IPAR)*zl1
            ENDIF
         ENDIF
C
      ENDIF
C
 999  CONTINUE
      END



