      REAL FUNCTION ENHITDIG(IDET,ICHA,IPAR,ID,IPROJ,IEN,ANORM,ECUT)
C
C===    Fills the histogram ID with the energies of the particles type IPAR
C===        ex:      IPAR=3  - e- (original particle)
C===                      -3 - e- (both original and secondary particles) 
C===                       0 - all particles
C===              hitting the detector IDET (IDET=0 - all detectors), channel ICHA (=0 - all channels)
C===              IPROJ :=ix+iy*10 <10 iy>0 - 2-dim, ix and iy stand for:
C===                     ix/iy   variable
C===                      1       momentum
C===                      2       radius of the hit
C===                      3       x of the hit
C===                      4       y of the hit
C===                      5       z of the hit (in DRS)
C===                      6       slope x of the hit (DRS)
C===                      7       slope y of the hit (DRS)
C===              IEN = 0    - weight 1.
C===                    1    - weight=energy
C===                    2    - weight=energy/(2*pi*rr/arc cm), arc=2cm
C===                    3    - weight=1./(2*pi*rr/ arc cm)
C===              ANORM      - multiplicative weight (for example normalizing the result to 1uA of beam)
C===              ECUT       - energy cut
C
      IMPLICIT NONE
      INTEGER  IDET,ICHA,IPAR,ID,IPROJ,IEN
      REAL     ECUT,ANORM
      INCLUDE ?
C
      LOGICAL HEXIST
C
      INTEGER ihi,nfind,k,itrh,iorh,idi,idg1,i1,i2,jhi,iddd,iphtyp
      INTEGER    mxp
      PARAMETER (mxp=7)
      REAL pp,rr,wg,parh(mxp),arc
C
      INTEGER nx,ny,nwt,lloc,ix,iy,imask
      REAL xx1,xx2,yy1,yy2,dx
      CHARACTER ctit*80
C
      INTEGER ifirst
      DATA ifirst/1/
      DATA arc/0./
C
C     -----------------------------------------------------------------
C
      ENHITDIG=0.
C
      IF(HEXIST(id)) THEN
         IF(ifirst.EQ.1) THEN
            CALL HRESET(ID,'    ')
            CALL HGIVE(id,ctit,nx,xx1,xx2,ny,yy1,yy2,nwt,lloc)
            dx=100.
            IF(nx.GT.0) dx=(xx2-xx1)/nx
            arc=dx
            IF(ny.GT.0) ny=1
            k=0
            IF(IPROJ.GT.10) k=1
            IF(ny.NE.k) THEN
               WRITE(6,*) ' *** Error: wrong dimension of ID=',id
               GO TO 999
            ENDIF
         ENDIF
      ENDIF
      ifirst=0
C
      iy=IPROJ/10
      ix=IPROJ-iy*10
C      write(6,*) 'ix,iy',ix,iy
      IF(ix.LE.0.OR.ix.GT.mxp.OR.iy.LT.0.OR.iy.GT.mxp) THEN
         IF(ifirst.NE.0) 
     +      WRITE(6,*) ' *** Error: wrong parameter IPROJ=',IPROJ
         GO TO 999
      ENDIF
C
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
            IF(IDET.EQ.0.OR.iddd.EQ.IDET) THEN
               IF(ICHA.EQ.0.OR.idg1.EQ.ICHA) THEN
                  i1=1
                  IF(idi.GT.1) THEN
                     k=ip1dig(idi-1)
                     i1=IBITS(k,16,16)
                  ENDIF
                  k=ip1dig(idi)
                  i2=IBITS(k,16,16)
C               WRITE(6,*) 'idg1=',idg1,jpdig(i1)
C
                  DO jhi=i1,i2-1
                     ihi=jpdig(jhi)
                     IF(ihi.GT.0.AND.ihi.LE.nhit) THEN
                        k=ip1hit(ihi)
                        itrh=IBITS(k,0,16)
                        iorh=IBITS(k,16,16)
C               WRITE(6,*) 'idg1=',idg1,itrh,itra(itrh)
                        IF(itrh.GT.0.AND.itrh.LE.ntra) THEN
                           iphtyp=itra(itrh)
                           IF(iorh.NE.0) iphtyp=iorh
                           IF(IPAR.EQ.0.OR.
     +                       (IPAR.GT.0.AND.iphtyp.EQ.IPAR
     +                                 .AND.iorh.EQ.0).OR.
     +                       (IPAR.LT.0.AND.iphtyp.EQ.IABS(IPAR))) THEN
                              pp=phit(ihi)
                              imask=1
                              IF(hit(2,ihi).LT.0.) THEN
                                 IF(ABS(hit(1,ihi)).LT.30.) imask=0
                              ELSE
                                 IF(ABS(hit(1,ihi)).LT.20.) imask=0
                              ENDIF
                              IF(pp.GT.ECUT.AND.imask.GE.0) THEN      ! mask
                                 nfind=nfind+1
                                 ENHITDIG=ENHITDIG+1.
                                 rr=SQRT(hit(1,ihi)**2+hit(2,ihi)**2)
                                 wg=ANORM
                                 IF(IEN.EQ.1) THEN
                                    wg=pp*wg
                                 ELSE IF(IEN.EQ.2) THEN
                                    wg=pp/2./3.1415/rr*ARC*wg
                                 ELSE IF(IEN.EQ.3) THEN
                                    wg=1./2./3.1415/rr*ARC*wg
                                 ENDIF
C     
                                 parh(1)=pp
                                 parh(2)=rr
                                 parh(3)=hit(1,ihi)
                                 parh(4)=hit(2,ihi)
                                 parh(5)=hitd(3,ihi)
                                 parh(6)=hitsld(1,ihi)
                                 parh(7)=hitsld(2,ihi)
C
                                 IF(iy.EQ.0) THEN
                                    CALL HF1(ID,parh(ix),wg)
                                 ELSE
                                    CALL HF2(ID,parh(ix),parh(iy),wg)
                                 ENDIF
C
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
C
               ENDIF
            ENDIF

         ENDDO
      ENDIF
C
 999  CONTINUE
      END



