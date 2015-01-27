      REAL FUNCTION SIGJPSI(ID)
C
C===        Analyze J/Psi + proton events and fills histograms
C===          ID - start histogram ID
C
      IMPLICIT NONE
      INTEGER  ID
C
      INCLUDE ?
C
      LOGICAL HEXIST
C
      INTEGER ihi,nfind,k,itrh,iorh,idi,idg1,i1,i2,jhi,iddd,iphtyp,id1
C
      INTEGER nx,ny,nwt,lloc,ix,iy,imask
      REAL xx1,xx2,yy1,yy2
      CHARACTER ctit*80
C
C      REAL RDLIM(2,2)  ! radii limits (in the calorimeter) 
      VECTOR RDLIM(2,2)  ! radii limits (in the calorimeter) 
C                          (1-2,k) - min-max, k=1 - e+-,gamma , k=2 -hadrons
C
      INTEGER    mxtrec,mxdet
      PARAMETER (mxtrec=10   ! max number of tracks to reconstruct
     +          ,mxdet=10)   ! max number of detectors to consider
      INTEGER itv,i,j,idet,itr,ipsirec,iprec,it,ipty,icpart
     +       ,ntrec,itrec(mxtrec)  ! nuber of tracks to reconstruct
     +       ,ndetr,idetr(mxdet)   ! number of detectors to use
     +       ,itrhit(mxdet,mxtrec) ! =ihit (or 0) for a give detector/track
     +       ,irec(mxtrec)         ! >0 - the track was reconstructed
C
      REAL    prec(5,mxtrec)       ! momenum/energy of the reconstructed tracks
     +    ,  eprec(mxtrec)         ! error on the full momentum
     +    ,hitcrec(2,mxtrec)       ! hits in the calorimeter for the reconstructed tracks
     +       ,rrec(mxtrec)         ! radius in the calorimeter for the reconstructed tracks
     +    , slorec(2,mxtrec)       ! slopes of the reconstructed particles
     +    ,eslorec(2,mxtrec)       ! errors of the slopes
     +    , xyzver(3)              ! reconstructed interaction point
     +    ,exyzver(3)              ! errors of the interaction point
     +    ,amass(100)              ! particle masses
C
      INTEGER ndpro,idpro(5)  ! decay products
      INTEGER idetcal,itrpsi,itrrp
C
      REAL pp,pph,rr,aa(2),xy(2),zdist,de,am,amjp,pjp(5),tt,ttrec,pps
      REAL resbxy                  ! beam x-y rexolution
     +    ,resencal(2,2)           ! calor. energy  resolution for the central and peripheral parts
     +    ,rspaccal(2)             ! calor. spacial resolution for the central and peripheral parts
     +    ,rcalor                  ! calorimeter - the radius of the central (high resol.) part
C
      REAl rmetric(4),plab(4,10),pcm(4),bet(4),costh
      INTEGER ifirst
      DATA ifirst/1/
      DATA rmetric/-1.,-1.,-1.,1./
      DATA idetcal/1/   ! calorimeter detector #
      DATA itrpsi/2/    ! J/Psi
      DATA itrrp/3/     ! recoil proton
      DATA resbxy/0.02/
      DATA resencal/0.030,0.010,0.060,0.015/
      DATA rspaccal/0.2,0.4/
      DATA rcalor/50./
C
C     -----------------------------------------------------------------
C
      SIGJPSI=0.
C
      DO i=1,100
         amass(i)=0.
      ENDDO
      amass(1)=0.
      amass(2)=0.000511
      amass(3)=0.000511
      amass(14)=0.938
      amass(83)=3.097
C
      IF(ifirst.EQ.1) THEN
         DO i=0,20
            id1=ID+i
            IF(HEXIST(id1)) THEN
               CALL HRESET(id1,'    ')
            ENDIF
         ENDDO
      ENDIF
      ifirst=0
C
      pps=SQRT(ptra(1,2)**2+ptra(2,2)**2+ptra(3,2)**2)
      CALL HF1(ID+16,pps,1.)
C
      DO it=1,MIN(ntra,10)
         iphtyp=itra(it)
         plab(4,it)=amass(iphtyp)**2
         DO i=1,3
            plab(i,it)=ptra(i,it)
            IF(it.EQ.1) plab(i,it)=-plab(i,it)
            plab(4,it)=plab(4,it)+plab(i,it)**2
         ENDDO
         plab(4,it)=SQRT(plab(4,it))
      ENDDO
C
C      write(6,FMT='(4F10.4)') plab
      tt=0.
      DO i=1,4
         tt=tt+rmetric(i)*(plab(i,2)-plab(i,1))**2
      ENDDO
      id1=ID+12
      CALL HF1(id1,tt,1.)
C
C      write(6,*) 'tt=',tt
C
C---     Decay angle in CM
C
      DO i=1,3
         bet(i)=plab(i,2)/plab(4,2)
      ENDDO
      bet(4)=1./SQRT(1.-(bet(1)**2+bet(2)**2+bet(3)**2))
      CALL GLOREN(bet(1),plab(1,4),pcm(1))
      costh=pcm(3)/SQRT(pcm(1)**2+pcm(2)**2+pcm(3)**2)
      id1=ID+14
      CALL HF1(id1,costh,1.)
C      write(6,*) 'bet=',bet
C
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      nfind=0
C
      IF(npdig.LT.1) GO TO 999
      IF(ntra.LT.3) GO TO 999
      IF(ntra.LT.itrpsi.OR.ntra.LT.itrrp) GO TO 999
C
      itv=itve(itrpsi)
      IF(nver.LT.itv) GO TO 999
C
      ndpro=ntdv(itv)
      IF(ndpro.LT.2) GO TO 999
C
      ntrec=0
      IF(ndpro+1.GT.mxtrec) THEN
         WRITE(6,*) ' *** Error - too many tracks to reconstruct '
     +             ,ndpro+1
         GO TO 999
      ENDIF
      DO i=1,ndpro
         idpro(i)=itdv(itv)-1+i
         ntrec=ntrec+1
         itrec(ntrec)=idpro(i)
      ENDDO
      ntrec=ntrec+1
      itrec(ntrec)=itrrp
C
      ndetr=1
      idetr(ndetr)=idetcal
C
      DO itr=1,ntrec
         DO idet=1,ndetr
            itrhit(idet,itr)=0
         ENDDO
      ENDDO
C
C---     Connect hits <==> tracks
C
      DO idi=1,ndig
         k=ip1dig(idi)
         iddd=IBITS(k,0,16)
         k=ip2dig(idi)
         idg1=IBITS(k,0,16)
         idet=0
         DO i=1,ndetr
            IF(idetr(i).EQ.iddd) idet=i
         ENDDO
         IF(idet.GT.0) THEN 
            i1=1
            IF(idi.GT.1) THEN
               k=ip1dig(idi-1)
               i1=IBITS(k,16,16)
            ENDIF
            k=ip1dig(idi)
            i2=IBITS(k,16,16)
C            WRITE(6,*) 'idg1=',idg1,jpdig(i1)
C
            DO jhi=i1,i2-1
               ihi=jpdig(jhi)
               IF(ihi.GT.0.AND.ihi.LE.nhit) THEN
                  k=ip1hit(ihi)
                  itrh=IBITS(k,0,16)
                  iorh=IBITS(k,16,16)
C                  WRITE(6,*) 'idg1=',idg1,itrh,itra(itrh)
                  IF(itrh.GT.0.AND.itrh.LE.ntra.AND.iorh.EQ.0) THEN
                     itr=0
                     DO i=1,ntrec
                        IF(itrec(i).EQ.itrh) itr=i
                     ENDDO
                     IF(itr.GT.0) THEN
                        itrhit(idet,itr)=ihi
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
C
         ENDIF
      ENDDO
C
C---     Did J/Psi products hit the calorimeter?
C
C      WRITE(6,*) 0,(itrhit(1,i),i=1,3)
      ipsirec=1
      DO itr=1,ndpro
         IF(itrhit(idetcal,itr).EQ.0) ipsirec=0
      ENDDO
      IF(ipsirec.EQ.0) GO TO 999
C
      CALL RANNOR(aa(1),aa(2))
      DO i=1,2
         xyzver(i)=vert(i,1)+aa(i)*resbxy
      ENDDO
      xyzver(3)=vert(3,1)
C      WRITE(6,*) 1,xyzver
C
C---     Reconstruct the tracks
C
      DO itr=1,ntrec
         it=itrec(itr)
         iphtyp=itra(it)
         ipty=2
         IF(iphtyp.LE.3) ipty=1    ! =1 - electron
         irec(itr)=0
         ihi=itrhit(idetcal,itr)
         pp=SQRT(ptra(1,it)**2+ptra(2,it)**2+ptra(3,it)**2)
C
         IF(ihi.NE.0) THEN
            rr=SQRT(hit(1,ihi)**2+hit(2,ihi)**2)
            pph=phit(ihi)
C
            IF(rr.GT.RDLIM(1,ipty).AND.rr.LT.RDLIM(2,ipty)) THEN
               icpart=1
               IF(rr.GT.rcalor) icpart=2    !   central part=1, peripheral=2
C
               zdist=zhit(ihi)-xyzver(3)
               CALL RANNOR(aa(1),aa(2))
               DO i=1,2
                  IF(ipty.EQ.2) aa(i)=aa(i)*5.
                  xy(i)=hit(i,ihi)+aa(i)*rspaccal(icpart)
                  slorec(i,itr)=(xy(i)-vert(i,1))/zdist
                  eslorec(i,itr)=(ptra(i,it)/pp-slorec(i,itr))
                  hitcrec(i,itr)=hit(i,ihi)
               ENDDO
               rrec(itr)=SQRT(hitcrec(1,itr)**2+hitcrec(2,itr)**2)
C
               CALL RANNOR(aa(1),aa(2))
               IF(ipty.EQ.1) THEN
                  de=(resencal(1,icpart)/SQRT(pph)
     +               +resencal(2,icpart))*aa(1)
                  prec(5,itr)=pph*(1.+de)
               ELSE
                  prec(5,itr)=1.
               ENDIF
               eprec(itr)=prec(5,itr)-pph
C
C               WRITE(6,*) itr,it,iphtyp,am,zdist
               prec(3,itr)=prec(5,itr)
     +                    /SQRT(1.+slorec(1,itr)**2+slorec(2,itr)**2)
               DO i=1,2
                  prec(i,itr)=slorec(i,itr)*prec(3,itr)
               ENDDO
C
               am=0.
               IF(iphtyp.GT.0.AND.iphtyp.LE.100) am=amass(iphtyp)
               prec(4,itr)=SQRT(prec(5,itr)**2+am**2)
C
               irec(itr)=1
               IF(ipty.EQ.1) THEN
                  IF(phit(ihi).GT.0.) THEN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
C      WRITE(6,*) 3,(irec(i),i=1,3)
C
C---     Reconstruct J/Psi mass
C
      DO i=1,4
         pjp(i)=0.
      ENDDO
C
      ipsirec=1
      DO itr=1,ndpro
         IF(irec(itr).NE.0) THEN
            DO i=1,4
               pjp(i)=pjp(i)+prec(i,itr)
            ENDDO
         ELSE
            ipsirec=0
         ENDIF
      ENDDO
      IF(ipsirec.EQ.0) GO TO 999
C
      amjp=pjp(4)**2-pjp(1)**2-pjp(2)**2-pjp(3)**2
      IF(amjp.GT.0.) amjp=SQRT(amjp)
C
      SIGJPSI=amjp
C
      ttrec=0.
      DO i=1,4
         ttrec=ttrec+rmetric(i)*(pjp(i)-plab(i,1))**2
      ENDDO
C
      id1=ID
      iphtyp=itra(itrpsi)
      am=amass(iphtyp)
      IF(am.GT.0.) THEN
         CALL HF1(id1,(amjp-am)/am,1.)
      ENDIF
      id1=ID+1
      CALL HF1(id1,rrec(3),1.)
      id1=ID+2
      CALL HF1(id1,tt,1.)
      id1=ID+3
      CALL HF2(id1,rrec(1),rrec(2),1.)
      id1=ID+4
      CALL HF1(id1,costh,1.)
      id1=ID+5
      CALL HF1(id1,ttrec-tt,1.)
      id1=ID+6
      CALL HF1(id1,pps,1.)
C
 999  CONTINUE
      END
C
      INCLUDE 'gloren.f'
