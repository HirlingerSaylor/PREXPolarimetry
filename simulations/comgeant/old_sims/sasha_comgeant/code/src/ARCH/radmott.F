      REAL FUNCTION RADMOTT(THETX,EEFRAC)
      IMPLICIT NONE
C
C===     Radiative Mott scattering on Fe
C        POLRAD code from Igor Akushevich
C        Terms:   1  - DIS
C                 2  - elastic on nucleus, radiative
C                 3  - quasielastic      , radiative
C                 4  - inelastic         , radiative
C
      REAL  THETX   ! - theta angle (rad) (or a function of theta)
     +     ,EEFRAC  ! - the fraction of the beam energy, left to the electron (radiative losses)
C
      COMMON/CCSMOTT/ CSDMOTT(4),ALMOTT
      REAL            CSDMOTT,ALMOTT
C---                  CSMOTT(1:3) = E0 electron (LAB), Z target, A target.
C---                  ALMOTT      - exponential slope of the variable
C
      REAL APOLRAD
C
#if defined OMGEANT_VERSION
      CHARACTER*80 VersionString
      DATA VersionString /
     & '$Id:$'/
#endif
C
      INTEGER npoint
      REAL thet,e0,zat,aat,sig,amh
     +    ,q2,xs,ys
C
      IF(ABS(ALMOTT).LT.1.E-6) THEN
         thet=THETX
      ELSE
         thet=-ALOG(THETX)/ALMOTT
      ENDIF
C
      e0= CSDMOTT(1)
      zat=CSDMOTT(2)
      aat=CSDMOTT(3)
C
      amh=0.938
      q2=4.*e0**2*EEFRAC*SIN(thet/2.)**2
      ys=1.-EEFRAC
      xs=q2/(2.*amh*ys*e0)
C
C---     GET THE CROSS SECTION IN 'barn/ster/GeV'
C
      sig=APOLRAD(e0,xs,ys)
C
      RADMOTT=sig*SIN(thet)  !  from d(sig)/d(cos_the) to d(sig)/d(thet)
      IF(ABS(ALMOTT).GT.1.E-6) THEN  ! in case of the LOG angle
         RADMOTT=RADMOTT/THETX/ALMOTT
      ENDIF
C
      WRITE(6,'(6E12.4)') E0,EEFRAC,thet,xs,ys,sig
      END
        
CDECK  ID>, POLRAD.
      REAL FUNCTION APOLRAD(ABEAM,AXS,AYS)
c
c      version 2.0  01.07.1996
c
      implicit real*8(a-h,o-z)
C
      REAL ABEAM,AXS,AYS
c
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich


      dimension tai(3),taip(3),tai2ll(3),si(2,3),si0(2,3),tls(2,3,4)
      parameter(nxfin=8)
      parameter(nyfin=4)
      dimension xmas(nxfin*nyfin),iiymas(nyfin)
      dimension ymas(nxfin*nyfin),ynet(20)
C      parameter(npoi=10)
C      dimension e00(npoi)
C      DATA E00 /.5, 1., 1.5, 2., 2.5, 3., 4., 5., 6. ,45./


	pl1=0.
	pn1=0.
	qn1=0.

C	 call titout('10.04.1997',bmom,tmom,pl1,pn1,qn1)


C      do 1 i=1,npoi
      bmom=ABEAM
      tmom=0.
      snuc=2.*(sqrt(tmom**2+amh*amh)*sqrt(bmom**2+aml2)+bmom*tmom)
C	xs= 0.5444E-03
      xs=AXS
C	ys=.5d0
      ys=AYS
	  y=snuc*xs*ys	 !  q2

C      write(9,998)xs,ys,snuc
C      print   998,xs,ys,snuc
      yma=1d0/(1d0+amp**2*xs/snuc)
      ymi=(amc2-amp**2)/(snuc*(1d0-xs))
      if(ys.gt.yma.or.ys.lt.ymi.or.xs.gt.1d0.or.xs.lt.0d0)then
C	write(9,66)
	goto 1
      endif
	 call conkin(snuc,amh)
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      call deltas(delta,delinf,tr)

      do 10 il=1,1
      do 10 in=1,1,2
       si(il,in)=0d0
       si0(il,in)=0d0
       tls(il,in,1)=0d0
       tls(il,in,2)=0d0
       tls(il,in,3)=0d0
       tls(il,in,4)=0d0

      isf1=1
      isf2=isf20
      isf3=1
	    if(in.eq.1)then
	       un=1.
	       pl=pl1
	       pn=0.
	       qn=0.
	       isf1=1
	       isf2=2
	    elseif (in.eq.3)then
	       un=0.
	       pl=pl1
	       pn=pn1
	       qn=0.
	       isf1=3
	       isf2=4
	    else
	       stop
	    endif

      do 30 ita=1,3
C      write(9,'(10(1h*),'' ita = '',i2,10(1h*))')ita

c
c     sib is born cross section with polarized initial
c     lepton and proton
c     sia is contribution of anomalous magnetic moment.
c
      if(ita.eq.1)then
	 call bornin(sib,sia)
      endif
c
c     tai(1),tai(2),tai(3) are contributions of radiative tails:
c		 1 - inelastic
c		 2 - elastic
c		 3 - quasielastic
c
      if(ita.eq.2) call conkin(snuc,amt)

	 call qqt(sib,tai(ita),taip(ita),tai2ll(ita))

      if(ita.eq.2) call conkin(snuc,amh)
  30  continue
      extai1=exp(alfa/pi*delinf)
      extai2=((sx-y/tara)**2/s/(s-y/tara))**tr
      extai3=((sx-y)**2/s/(s-y))**tr
      sig=sib*extai1*(1.+alfa/pi*(delta-delinf))+sia
     . +tai(1)+(tai(2)*extai2+tai(3)*extai3)/tara
      si(il,in)=si(il,in)+sig
      si0(il,in)=si0(il,in)+sib
      tls(il,in,1)=tls(il,in,1)+tai(1)
      tls(il,in,2)=tls(il,in,2)+tai(2)/tara
      tls(il,in,3)=tls(il,in,3)+tai(3)/tara
      tls(il,in,4)=tls(il,in,4)+(sib*alfa/pi*delta+sia)
      silev=tai(2)/tara/1d9*(x/ap)/pi/sx
      APOLRAD=silev*tara
C      if(in.eq.1)write(6,'(16g13.4)')
C     .	xs,ys,s/ap,silev,silev*tara
C     .	 ,tai(2)/tara!,tai(3)/tara,tai2ll(1)
c    .	 ,tai2ll(2)/tara,tai2ll(3)/tara,taip(1)
c    .	 ,taip(2)/tara,taip(3)/tara
C      if(in.eq.3)write(23,'(16g11.5)')
C     .	xs,ys,s,sib,sig,delta,tai(1),tai(2)/tara,tai(3)/tara
C     .	 ,tai2ll(1),tai2ll(2)/tara,tai2ll(3)/tara
C     .	 ,taip(1),taip(2)/tara,taip(3)/tara
 10   continue
c		 d factor

      ddf=ys*(2.-ys)/(ys*ys+2.*(1.-ys)*(1.+r1990(xs,y)))
c
c     si (si0) is dis cross section for the corresponding target
c     aspop is polarized asymmetry including radiative corrections
c     asbor is born polarized asymmetry
c     del is total radiative correction
c
       aspop=-si(1,3)/si(1,1)/ddf
       asbor=-si0(1,3)/si0(1,1)/ddf
       del= aspop-asbor
       tip=tls(1,3,1)/si(1,1)/ddf/asbor
       tep=tls(1,3,2)/si(1,1)/ddf/asbor
       tqp=tls(1,3,3)/si(1,1)/ddf/asbor
       tvp=tls(1,3,4)/si(1,1)/ddf/asbor
       tiu=tls(1,1,1)*si0(1,3)/si0(1,1)/si(1,1)/ddf/asbor
       teu=tls(1,1,2)*si0(1,3)/si0(1,1)/si(1,1)/ddf/asbor
       tqu=tls(1,1,3)*si0(1,3)/si0(1,1)/si(1,1)/ddf/asbor
       tvu=tls(1,1,4)*si0(1,3)/si0(1,1)/si(1,1)/ddf/asbor

c     if ( abs(asbor).gt.1d-10) del=(aspop-asbor)/asbor
      del= aspop-asbor
C      write(9,68)xs,si0(1,1),si0(1,2),si0(1,3),asbor
C      write(9,68)xs,si (1,1),si (1,2),si (1,3),aspop
C      write(7,69)xs,w2,y,asbor,aspop,del
C      write(16,'(2f6.3,0pf7.2,7f7.2)')xs,ys,tip,tep,tqp,
C     .tvp,tiu,teu,tqu,tvu
   1  continue



  66  format(1x,'kinematics')
  68  format(1x,f7.3,3e15.5,2pf15.5)
  69  format(1x,f7.3,2f8.1,2pf8.3,2f8.3)
 998  format(' ****** x = ',f7.5,' y = ',f7.5,' s = ',g13.4,' ******')
  72  format(//' ******   iteration procedure result table    ******')

1000  end


      block data
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/p/pi,pi2,alfa,i1(8),i2(8)
      data
     .amm/2.7928456d0/,amn/-1.913148d0/,chbar/.197328d0/,barn/.389379d6/
     .aml/.511000d-3/,aml2/.261112d-6/,al2/.522240d-6/,
     .amt/52.5432d0/,
     .tara/56d0/
     .tarz/26d0/
     .fermom/.26d0/
     .isf20/2/
     .pi/3.1415926d0/,pi2/9.869604d0/,alfa/.729735d-2/,amc2/1.151857d0/,
     .amp/.938272d0/,amh/.938272d0/,
     .i2/1,1,1,2,3,3,1,2/,i1/3,3,4,4,3,3,3,3/
      end


CDECK  ID>, CONKIN.
****************** conkin *************************************

      subroutine conkin(snuc,amtar)
c set of kinematical constants
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      amp=amtar
      ap=2.*amp
      amp2=amp**2
      ap2=2.*amp**2
      s=snuc*amp/amh
      x=s*(1.-ys)
      sx=s-x
      sxp=s+x
      ym=y+al2
      tpl=s**2+x**2
      tmi=s**2-x**2
      w2=amp2+s-y-x
      als=s*s-al2*ap2
      alx=x*x-al2*ap2
      alm=y*y+4.*aml2*y
      aly=sx**2+4.*amp2*y
      sqls=dsqrt(als)
      sqlx=dsqrt(alx)
      sqly=dsqrt(aly)
      sqlm=dsqrt(alm)
      allm=dlog((sqlm+y)/(sqlm-y))/sqlm
      axy=pi*(s-x)
      an=2.*alfa**2/sqls*axy*barn*amh/amp
c      tamin=(sx-sqly)/ap2
      tamax=(sx+sqly)/ap2
      tamin=-y/amp2/tamax
      as=s/2./aml/sqls
      bs=0.
      cs=-aml/sqls
	 ae=amp/sqls
	 be=0.
	 ce=-s/ap/sqls
      apq=-y*(ae-be)+ce*sx
      apn=(y+4.*aml2)*(ae+be)+ce*sxp
      dk2ks=as*ym+al2*bs+cs*x
      dksp1=as*s+bs*x+cs*ap2
      dapks=2.*(al2*(as*ae+bs*be)+ap2*cs*ce+ym*(as*be+bs*ae)
     .+s*(as*ce+cs*ae)+x*(bs*ce+cs*be))
      return
      end

CDECK  ID>, BORNIN.
****************** bornin *************************************

      subroutine bornin(sibor,siamm)
c
c     sibor is born cross section with polarized initial
c     lepton and polarized target
c     siamm is contribution of anomalous magnetic moment.
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
       common/print/ipri1
      dimension sfm0(8),tm(8)
       ipri1=1
      call strf(0d0,0d0,sfm0)
       ipri1=0
      tm(1)=-(2.*aml2-y)
      tm(2)=(-(amp2*y-s*x))/(2.*amp2)
      tm(3)=(2.*(apq*dk2ks-dapks*y)*aml)/amp
      tm(4)=apq/amp*(-(dk2ks*sx-2.*dksp1*y)*aml)/amp2
      tm(7)=(-(4.*aml2+3.*apn**2-3.*apq**2+y))/2.
      tm(8)=apq/amp*(-3.*(apn*sxp-apq*sx))/(2.*ap)
      ek=(3.*apq**2-y)/amp2
      tm(5)=-ek*tm(1)
      tm(6)=-ek*tm(2)
      ssum=0.
      do 1 isf=isf1,isf2,isf3
	ppol=1.
	if(isf.eq.3.or.isf.eq.4)ppol=-pn
	if(isf.ge.5)ppol=qn/6
      ssum=ssum+tm(isf)*sfm0(isf)*ppol
    1 continue
      sibor=ssum*an/y/y*2.
c
c  formula (4) of kukhto and shumeiko paper
c
cc    res1=amp*ww1*(y+4.*aml2)-ww2*(s+x)**2/4./amp
cc    siamm=alfa/pi*al2*allm*(sibor+an*res1/y**2)
      siamm=0.
      return
      end

CDECK  ID>, DELTAS.
****************** deltas *************************************

      subroutine deltas(delta,delinf,tr)
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      del1=-ym*(alm*allm**2/2.+2.*fspen(2d0*sqlm/(y+sqlm))-pi2/2.)/sqlm
      del2=(3.*y/2.+4.*aml2)*allm-2.

      sum=vacpol(y)

      aj0=2.*(ym*allm-1.)
      deltai=aj0*dlog((w2-amc2)/aml/dsqrt(w2))

      ss=x+y
      xx=s-y
      alss=ss**2-2.*w2*al2
      alxx=xx**2-2.*w2*al2
      sqlss=dsqrt(alss)
      sqlxx=dsqrt(alxx)
      allss=dlog((sqlss+ss)/(-sqlss+ss))/sqlss
      allxx=dlog((sqlxx+xx)/(-sqlxx+xx))/sqlxx
      dlm=dlog(y/aml2)
	 sfpr=dlm**2/2.-dlm*dlog(ss*xx/aml2/w2)
     . -(dlog(ss/xx))**2/2.+fspen((s*x-y*amp2)/ss/xx)-pi2/3.
      delta0=(ss*allss+xx*allxx)/2.+sfpr
      delta=deltai+delta0+del1+del2+sum
      delinf=(dlm-1.)*dlog((w2-amc2)**2/ss/xx)
      tr=alfa/pi*(dlm-1.)
      return
      end

CDECK  ID>, VACPOL.
      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/p/pi,pi2,alfa,i1(8),i2(8)
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alfa

      vacpol=suml+sumh

      end

CDECK  ID>, FSPENS.
****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end
CDECK  ID>, FSPEN.
****************** fspen **************************************

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end

CDECK  ID>, QQT.
****************** qqt ****************************************

      subroutine qqt(bo,tai,taip,tai2ll)
      implicit real*8(a-h,o-z)
      external rv2ln,rv2
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      ep=1d-8
      abb0=bo*pi/an/alfa
      tlnmin=dlog(tamin+y/sx)
      tlnmax=dlog(tamax+y/sx)
      tln0=dlog(y/sx)
      tlns=dlog(-y/s+y/sx)
      tlnp=dlog(y/x+y/sx)
      if(ita.ne.1)then
c      call qunc8(rv2,tamin,tamax,ep*abb0,ep,res2,er,nn2,fl2,3500)
c      call qunc8(rv2ln,tlnmin,tlnmax,ep*abb0,ep,res2,er,nn2,fl2,3500)
      call simpxx(tlnmin,tln0,200,1d-4,rv2ln,res1)
      call simpxx(tln0,tlnmax,100,1d-2,rv2ln,res2)
	  tai=an*alfa/pi*(res1+res2)
c	  write(*,9)tai,fl2,nn2
C	  write(9,9)tai   !,fl2,nn2
	  taip=0.

      else

	  call qqint( 500,26000,3,1,5d-4,tai)


      end if
c   9  format(' tai = ',e12.4,'   fl = ',f6.2,'     nn = ',i4)
   9  format(' tai = ',e12.4)

      end

CDECK  ID>, TAILS.
****************** tails **************************************

       subroutine tails(ta,tm)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
       dimension tm(8,6),ajm2(2),ajm3(3),ii(8)
      data ii/1,2,3,4,1,2,5,6/
      b2=(-aly*ta+sxp*sx*ta+2.*sxp*y)/2.
      b1=(-aly*ta-sxp*sx*ta-2.*sxp*y)/2.
      c1=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(s*ta+y)**2)
      c2=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(ta*x-y)**2)
       bb=1./sqly
       sc1=dsqrt(c1)
       sc2=dsqrt(c2)
       bi12=(sxp*(sx*ta+2.*y))/(sc1*sc2*(sc1+sc2))
       bi1pi2=1./sc2+1./sc1
       bis=-b1/sc1/c1+b2/sc2/c2
       bir=b2/sc2/c2+b1/sc1/c1
       b1i=-b1/aly/sqly
       b11i=(3.*b1**2-aly*c1)/2./aly**2/sqly
	 sps=as+bs
	 spe=ae+be
	 ccpe=(ae-be)*ta+2.*ce
	 ccps=(as-bs)*ta+2.*cs
      sis=(2.*bi1pi2*sps+bir*sps*ta+bis*ccps)/2.
      sir=( (2.*bi12*sps*ta+bir*ccps+bis*sps*ta))/2.
      si12=(bi12*ccps+bi1pi2*sps)/2.
      eis=(2.*bi1pi2*spe+bir*spe*ta+bis*ccpe)/2.
      eir=( (2.*bi12*spe*ta+bir*ccpe+bis*spe*ta))/2.
      ei12=(bi12*ccpe+bi1pi2*spe)/2.
      ois=((2.*bi1pi2+bir*ta)*(ccpe*sps+ccps*spe)+(ccpe*ccps+
     . spe*sps*ta**2)*bis+8.*bb*spe*sps+4.*bi12*spe*sps*ta**2)/
     . 4.
      oir=( ((2.*bi12+bis)*(ccpe*sps+ccps*spe)*ta+(ccpe*ccps+
     . spe*sps*ta**2)*bir+4.*bi1pi2*spe*sps*ta))/4.
      oi12=((ccpe*ccps+spe*sps*ta**2)*bi12+(ccpe*sps+ccps*spe)*
     . bi1pi2+4.*bb*spe*sps)/4.
      eeis=((ccpe**2+spe**2*ta**2)*bis+8.*bb*spe**2+4.*bi12*spe
     . **2*ta**2+4.*bi1pi2*ccpe*spe+2.*bir*ccpe*spe*ta)/4.
      eeir=( ((ccpe**2+spe**2*ta**2)*bir+4.*bi12*ccpe*spe*ta+4.
     . *bi1pi2*spe**2*ta+2.*bis*ccpe*spe*ta))/4.
      eei12=((ccpe**2+spe**2*ta**2)*bi12+4.*bb*spe**2+2.*bi1pi2
     . *ccpe*spe)/4.
      ei1pi2=(4.*bb*spe+bi12*spe*ta**2+bi1pi2*ccpe)/2.
      eei1i2=((ccpe**2+spe**2*ta**2)*bi1pi2+4.*(2.*ccpe-spe*ta)
     . *bb*spe+8.*b1i*spe**2+2.*bi12*ccpe*spe*ta**2)/4.
      eb=((ccpe-spe*ta)*bb+2.*b1i*spe)/2.
      eeb=((ccpe-spe*ta)**2*bb+4.*(ccpe-spe*ta)*b1i*spe+4.*b11i
     . *spe**2)/4.
       call ffu(1,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .,eis,eir,ei12,ei1pi2,ta)
       call ffu(2,eb,eis,eir,ei12,ei1pi2,oir,ois,oi12
     .,eeis,eeir,eei12,eei1i2,ta)
       call ffu(3,eeb,eeis,eeir,eei12,eei1i2,0d0,0d0,0d0
     .,0d0,0d0,0d0,0d0,ta)
       ajm2(1)=apq/amp
       ajm2(2)=-1./amp
       ajm3(1)=(y-3.*apq**2)/amp2
       ajm3(2)=6.*apq/amp2
       ajm3(3)=-3./amp2
       do 15 i=1,8
       do 13 l=1,6
   13  tm(i,l)=0
       do 10 k=1,i2(i)
       ajk=1.
       if(i.eq.4.or.i.eq.8)ajk=ajm2(k)
       if(i.eq.5.or.i.eq.6)ajk=ajm3(k)
       do 10 j=k,i1(i)+k-1
       tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,k)*ajk
       if((i.eq.5.or.i.eq.6).and.k.eq.2)
     . tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,1)*ta/amp2
  10   continue
  15   continue
       return
       end

CDECK  ID>, FFU.
****************** ffu ****************************************

       subroutine ffu(n,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .	      ,eis,eir,ei12,ei1pi2,ta)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
      hi2=aml2*bis-ym*bi12
      shi2=aml2*sis-ym*si12
      ehi2=aml2*eis-ym*ei12
      ohi2=aml2*ois-ym*oi12
       goto(10,20,30)n
  10   continue
      tm3(3,1,n)=(8.*(apq*dk2ks-dapks*y)*aml*hi2)/amp
      tm3(3,2,n)=(-2.*((2.*(bi12*dk2ks*ta-2.*shi2)*apq+(2.*shi2-
     . sir*y+sis*ym)*apn+4.*dapks*hi2*ta)-4.*((2.*ei12-eis)*
     . dk2ks-(si12-sis)*apn)*aml2)*aml)/amp
      tm3(3,3,n)=(2.*(((2.*si12+sir-sis)*apn*ta-2.*dk2ks*ei12*ta
     . -6.*ohi2-oir*y+ois*ym)-4.*aml2*oi12)*aml)/amp
      tm3(3,4,n)=(2.*(2.*oi12-oir+ois)*aml*ta)/amp
      tm3(5,1,n)=-2.*(4.*aml2+3.*apn**2-3.*apq**2+y)*hi2
      tm3(5,2,n)=-2.*(6.*aml2*apn*eir-3.*apn**2*bi12*ta+3.*apn*
     . apq*bi1pi2+6.*apq*ehi2+hi2*ta)
      tm3(5,3,n)=-(24.*aml2*eei12-6.*apn*ei1pi2-6.*apq*ei12*ta-
     . 2.*bb-bi12*ta**2)
  20   continue
      tm3(4,1,n)=(-4.*(dk2ks*sx-2.*dksp1*y)*aml*hi2)/amp2
      tm3(4,2,n)=(((2.*(sxp-2.*sx)*shi2+2.*bi12*dk2ks*sx*ta+8.*
     . dksp1*hi2*ta-sir*sxp*y+sis*sxp*ym)-4.*(2.*bi12*dk2ks-bis*
     . dk2ks-si12*sxp+sis*sxp)*aml2)*aml)/amp2
      tm3(4,3,n)=((((sxp*ta-ym)*sis-(sxp*ta-y)*sir+2.*bi12*dk2ks
     . *ta+6.*shi2-2.*si12*sxp*ta)+4.*aml2*si12)*aml)/amp2
      tm3(4,4,n)=(-(2.*si12-sir+sis)*aml*ta)/amp2
      tm3(6,1,n)=(-3.*(apn*sxp-apq*sx)*hi2)/amp
      tm3(6,2,n)=(-3.*(2.*(apn*bir+eir*sxp)*aml2-(2.*bi12*sxp*ta
     . -bi1pi2*sx)*apn+(bi1pi2*sxp+2.*hi2)*apq+2.*ehi2*sx))/(2.*
     . amp)
      tm3(6,3,n)=(-3.*(8.*aml2*ei12-apn*bi1pi2-apq*bi12*ta-ei12*
     . sx*ta-ei1pi2*sxp))/(2.*amp)
  30   continue
      tm3(1,1,n)=-4.*(2.*aml2-y)*hi2
      tm3(1,2,n)=4.*hi2*ta
      tm3(1,3,n)=-2.*(2.*bb+bi12*ta**2)
      tm3(2,1,n)=(((sxp**2-sx**2)-4.*amp2*y)*hi2)/(2.*amp2)
      tm3(2,2,n)=(2.*aml2*bir*sxp-4.*amp2*hi2*ta-bi12*sxp**2*ta+
     . bi1pi2*sxp*sx+2.*hi2*sx)/(2.*amp2)
      tm3(2,3,n)=(2.*(2.*bb+bi12*ta**2)*amp2+4.*aml2*bi12-bi12*
     . sx*ta-bi1pi2*sxp)/(2.*amp2)
       return
       end


CDECK  ID>, QQINT.
****************** qqint **************************************

      subroutine qqint(mi,ma,it,ir,ot,res)
      implicit real*8(a-h,o-z)
      external rv2di,rv2dia
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/tttrr/ttr
      dimension am(2),bm(2),tlm(4),rrm(4),wrk(500)
      tlm(1)=log(xs+tamin)
      tlm(4)=log(xs+tamax)
      tlm(2)=log(xs-y/s)
      tlm(3)=log(xs+y/x)
      rrm(1)=1d-10
      rrm(4)=w2-amc2-1d-10
      rrm(2)=w2-1.5**2
      rrm(3)=w2-1.215**2
      res=0.
      do 10 i=1,it
	 ii=i+1
	 if(i.eq.it)ii=4
	 am(1)=tlm(i)
	 bm(1)=tlm(ii)
	   do 10 j=1,ir
	     jj=j+1
	     if(j.eq.ir)jj=4
	     if(rrm(jj).le.rrm(j))then
		rrm(jj)=rrm(j)
		tai=0.
C		goto 11
	     endif
	     am(2)=rrm(j)
	     bm(2)=rrm(jj)
	     id=1
	     mir=mi
      call d01fce(2,am,bm,mir,ma,rv2di,ot,otr,500,wrk,re,id)
c     call simpdo(am(1),bm(1),1d-2,100,am(2),bm(2),1d-3,200,rv2dia,re)
	     tai=an*alfa/pi*re
	     res=res+tai
C  11  write(9,'(1x,''tai:'',2i3,2g13.4,2i5)')i,j,tai,otr,mir,id
c     write(*,'(1x,''tai:'',2i3,2g13.4,2i5)')i,j,tai,otr,mir,id
  10  continue
C      write(9,'(1x,''result:'',g13.4)')res
c     write(*,'(1x,''result:'',g13.4)')res
      return
      end

CDECK  ID>, RV2DI.
****************** rv2di **************************************

      double precision function rv2dia(xx)
      implicit real*8(a-h,o-z)
      common/simpc/x
      dimension z(2)
      z(1)=x
      z(2)=xx
      rv2dia=rv2di(2,z)
      end

      double precision function rv2di(ndim,z)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tttrr/ttr
      dimension sfm(8),sfm0(8),tm(8,6)
      dimension z(ndim)

      ta=ddexp(z(1))-xs
      r=z(2)/(1.+ta)
      call strf(0d0,0d0,sfm0)
      call tails(ta,tm)
      call strf(ta ,r,sfm)

      podinl=0.
      do 11 isf=isf1,isf2,isf3
	ppol=1.
	if(isf.eq.3.or.isf.eq.4)ppol=-pn
	if(isf.ge.5)ppol=qn/6
      do 1 irr=1,i1(isf)+i2(isf)-1
	pp=sfm(isf)
      if(irr.eq.1.and.ita.eq.1)pp=pp-sfm0(isf)*(1.+r*ta/y)**2
      pres=pp*r**(irr-2)/(y+r*ta)**2/2.
      podinl=podinl-tm(isf,irr)*pres*ppol
    1 continue
   11 continue
      rv2di=podinl*(xs+ta)/(1.+ta)
c     rv2di=0.
   3  format(3e12.6,i6)
      return
      end

CDECK  ID>, RV2LN.
****************** rv2ln **************************************

      double precision function rv2ln(taln)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      external podinl
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/amf2/taa,tm(8,6),sfm0(8)


      common/p/pi,pi2,alfa,i1(8),i2(8)


      ta=ddexp(taln)-y/sx
      taa=ta
cccc  call strf(0d0,0d0,sfm0)
       call tails(ta,tm)
      rmin=1d-8
      rmax=(w2-amc2)/(1.+ta)
      if(ita.eq.1)then
c	 call qvnc8(podinl,rmin,rmax,1d-4,1d-9,res,er,nn,fl,3500)
	 call dqn32(rmin,rmax,podinl,res)
      else
	 aa=amt/amp
	 if(ita.eq.3)aa=1.
	 res=podinl((sx-y/aa)/(1d0+ta/aa))/(1.+ta/aa) /aa**2
      endif
      rv2ln=res*(y/sx+ta)

      return
      end

CDECK  ID>, PODINL.
****************** podinl *************************************

      double precision function podinl(r)
c
c     integrand (over r )
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/amf2/ta,tm(8,6),sfm0(8)
      dimension sfm(8),iel(8)
      data iel/0,2,1,2,2,4,2,3/
	call strf(ta ,r,sfm)
      podinl=0.
      do 11 isf=isf1,isf2,isf3
	ppol=1.
	if(isf.eq.3.or.isf.eq.4)ppol=-pn
	if(isf.ge.5)ppol=qn/6
	if(ita.eq.2)ppol=ppol*(amt/amp)**iel(isf)
      do 1 irr=1,i1(isf)+i2(isf)-1
	pp=sfm(isf)
      if(irr.eq.1.and.ita.eq.1)pp=pp-sfm0(isf)*(1.+r*ta/y)
      pres=pp*r**(irr-2)/(y+r*ta)**2/2.
      podinl=podinl-tm(isf,irr)*pres*ppol
    1 continue
   11 continue

      aa=1.
      t=y+ta*(sx-y/aa)/(1d0+ta/aa)
      sxsx=(s**2+x**2)/s/x
      si=sxsx*(sfm(1)- .5 *sfm(2))
c     write (10,'(1x,2f15.8,5f13.4)')ta,t,podinl,si
      t1=tm(1,1)/r+tm(1,2)+r*tm(1,3)
      t2=tm(2,1)/r+tm(2,2)+r*tm(2,3)
c     write (10,'(1x,2f15.8,5f13.4)')ta,t,sxsx,t1,t2

       return
      end

CDECK  ID>, RV2.
****************** rv2 ****************************************

      double precision function rv2(ta)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      external podinl
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/amf2/taa,tm(8,6),sfm0(8)


      common/p/pi,pi2,alfa,i1(8),i2(8)


      taa=ta
cccc  call strf(0d0,0d0,sfm0)
       call tails(ta,tm)
      rmin=1d-8
      rmax=(w2-amc2)/(1.+ta)
      if(ita.eq.1)then
c	 call qvnc8(podinl,rmin,rmax,1d-4,1d-9,res,er,nn,fl,3500)
	 call dqn32(rmin,rmax,podinl,res)
      else
	 aa=amt/amp
	 if(ita.eq.3)aa=1.
	 res=podinl((sx-y/aa)/(1d0+ta/aa))/(1.+ta/aa) /aa**2
      endif
      rv2=res
      return
      end


CDECK  ID>, STRF.
      subroutine strf(ta,rr,sfm)
c
c     the programm calculates deep inelastic (ita=1),
c     elastic (ita=2), quasielastic (ita=3) structure functions
c     in kinematical point (ta,rr).
c	   rr=sx-tt,
c	   ta=(t-y)/rr,
c    where tt=t+amf2-amp2, amf2 is invarint mass of final hadrons
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
       common/print/ipri1
      dimension sfm(8)
      t=y+rr*ta
      if(ita.eq.1)then
      tt=sx-rr
      amf2=tt-t+amp2
      aks=t/tt
      anu=tt/ap
      epsi=ap2/tt
      g1=0.d0
      g2=0.d0
      if(un.gt.1d-10)then
       f1=f1sfun(aks,t)
       f2=f2sfun(aks,t)
      else
       f1=0.d0
       f2=0.d0
      endif
      if(abs(pn).gt.1d-10)then
       g1=g1sfun(aks,t)
       g2=g2sfun(aks,t)
      else
       g1=0.d0
       g2=0.d0
      endif
       b1=0.d0
       b2=0.d0
       b3=0.d0
       b4=0.d0

      goto 10
      endif
c
c   tarn,tarz,tara are n,z,a of the nucleus
c
      epsi=ap2/t
      tarn=tara-tarz
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   gep,gmp are electric and magnetic form factors of proton
c   s.i.bilenkaya et al pisma zhetf 19(1974)613
c
      call ffpro(t,gep,gmp)
      if(ita.eq.2)then
      tau=t/4./amt**2
      tau1=1.+tau
	  call ffco(t,ff)
	   f1=0d0
	   f2=4.*amp2*tau*(tarz*ff)**2
	   g1=0d0
	   g2=0d0
       else if(ita.eq.3)then
	   tau=t/4./amp**2
	   tau1=1.+tau
	   call ffquas(t,geun,gmun,gepo,gmpo)
	   f1=2.*amp2*tau*gmun**2
	   f2=4.*amp2*tau*(geun**2+tau*gmun**2)/tau1
	   g1=amp2*tau*2.*gmpo*(gepo+tau*gmpo)/tau1
	   g2=2.*amp2*tau**2*gmpo*(gepo-gmpo)/tau1
      endif
  10  continue
      sfm(1)=un*f1+qn/6.*b1
      sfm(2)=epsi*(un*f2+qn/6.*b2)
      sfm(3)=epsi*(g1+g2)
      sfm(4)=epsi**2*g2
      sfm(5)=epsi**2*b1
      sfm(6)=epsi**3*(b2/3.+b3+b4)
      sfm(7)=epsi*(b2/3.-b3)
      sfm(8)=epsi**2*(b2/3.-b4)
      return
      end


CDECK  ID>, F1SFUN.
********************** f1sfun ***********************************
      double precision function f1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      f2=f2sfun(aks,t)
      anu=t/ap/aks
      f1sfun=amp*(1.+anu**2/t)*f2/anu/(1.+r1990(aks,t))
      end

CDECK  ID>, F2SFUN.
********************** f2sfun ***********************************
      double precision function f2sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      f2p=df2h8(t,aks)
      f2d=df2d8(t,aks)
c      f2=f2d*(-.45*aks+1.07-.45*ddexp(-44.00d0*aks))
       f2sfun=f2d*ranucl(aks,tara)
      end


CDECK  ID>, R1990.
********************** r1990 ************************************

      double precision function r1990(aks,tc)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      t=tc
      if(tc.lt..35d0)t=0.35
      teta=1.+12.*t/(1.+t)*(.125**2/(aks**2+.125**2))
      zn=teta/log(t/.04)
      ra=.0672*zn+.4671/(t**4+1.8979**4)**(.25d0)
      rb=.0635*zn+.5747/t-.3534/(t**2+.09)
      rc=.0599*zn+.5088/sqrt((t-5.*(1.-aks)**5)**2+2.1081**2)
      rrr=(ra+rb+rc)/3.
      r1990=rrr
      return
      end


CDECK  ID>, G1SFUN.
********************** g1sfun ***********************************
      double precision function g1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20

      end


CDECK  ID>, G2SFUN.
********************** g2sfun ***********************************
      double precision function g2sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      g2sfun=0.
      end
CDECK  ID>, RANUCL.
       double precision function ranucl(x,an)
      implicit real*8(a-h,o-z)
c
c  'a' dependence of the structure functions ratio
c   r=f2(nucleus)/f2(deuteron) is represented  by formulas
c   of barshay and rein, z.phys.c 46 (1990) 215,
c   with the normalization parameter bnorm=1
c
c   'x' dependence of the ratio r and normalisation parameters
c   pam(1-3) for a  dependence	were found from the fit of
c   he, li, c, ca, cu, xe and pb data in the range 0.001< x< 0.7
c   done in june 1995 by g.smirnov: phys. lett. b364 (1995) 87.
c
       dimension pam(3), pamer(3)
       data   pam/ 0.130, 0.456, 0.773/
       data pamer/ 0.004, 0.017, 0.020/
       data b1,b2,b3,b4,b5/1., 1.145, 0.93, 0.88, 0.59/
c
       if(an.le.4.5)then
	 pam(1)=0.02934
	 pam(2)=0.1046
	 pam(3)=0.2360
       endif

       bnorm = 1.0
	d1 = an**(1./3.)
	d2 = an**(2./3.)
	d3 = an
	d4 = an**(4./3.)
	d5 = an**(5./3.)

	ser = 1. - b1/d1 -b2/d2 + b3/d3 + b4/d4 - b5/d5
c
	del1 = bnorm*ser
c
       aa = x**(pam(1)*del1)
       bb = 1.0 + pam(2)*del1
       cc = x*pam(3)*del1
c
       ranucl = aa*bb*(1.0  - cc)

      x0=0.8
      cb=15.
      db=1.
c      ranucl=ranucl+(-x0**pam(1)*(1.+pam(2))*(1.-pam(3)*x0)+1.)
c     .    *exp(-cb*(x0-x)**db)
c      if(x.gt.0.79 .and. ranucl.gt.1.)ranucl=1.

       end


CDECK  ID>, FFPRO.
****************** ffpro **************************************

      subroutine ffpro(t,gep,gmp)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      gep=1.2742/(1.+t/0.6394**2)-.2742/(1.+t/1.582**2)
      gmp=(1.3262/(1.+t/0.6397**2)-.3262/(1.+t/1.3137**2))*amm
c     gep=1./((1.+.61*t)*(1.+2.31*t)*(1.+.04*t))
c     gmp=amm*gep
      end

CDECK  ID>, FFHE3.
****************** ffhe3 **************************************

      subroutine ffhe3(t,ge,gm)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      tf=t/chbar**2
       qf=sqrt(tf)
      a=.675
      b=.366
      c=.836
      am=.654
      bm=.456
      cm=.821
      d=-6.78d-3
      p=.9
      q0=3.98
      f0=ddexp(-a**2*qf**2) - b**2*qf**2*ddexp(-c**2*qf**2)
      fm=ddexp(-am**2*qf**2) - bm**2*qf**2*ddexp(-cm**2*qf**2)
      df=d*ddexp(-((qf-q0)/p)**2)
      ge=(f0+df)*tarz
      gm=fm*tara      * (-2.13)
      end

CDECK  ID>, FFCO.
****************** ffco **************************************

      subroutine ffco(t,ff)
c
c   j.e. bailey et al nucl. phys. b151(1979)367
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
	   tf=t/chbar**2
	   am=(tarz-2.)/3.
	   xm=1.07*dsqrt(tf)*tara**(1./3.)
	   um=3.*(2.+5.*am)/(2.*(2.+3.*am))
	   ff=(1.-am*xm**2/(2.*um*(2.+3.*am)))*ddexp(-xm**2/(4d0*um))
	   ff=max(ff,0d0)
      end

CDECK  ID>, FFQUAS.
****************** ffquas **************************************

      subroutine ffquas(t,geun,gmun,gepo,gmpo)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      call ffpro(t,gep,gmp)
      tf=t/chbar**2
      tau=t/4./amp**2
      tau1=1.+tau
c
c   t. de forest and d.j. valecka adv. in phys. 15(1966) no.57
c
	supele=1.
	supmag=1.
	     if(t.le.(2.d0*fermom)**2)then
	       sqrat=dsqrt(t)/fermom
	       supele=0.75*sqrat-sqrat**3/16.
	       supmag=supele
	     endif
	geun=gep*dsqrt(supele*tarz)
	tarn=tara-tarz
	gmun=gep*dsqrt(supmag*(tarz*amm**2+tarn*amn**2))
	gepo=0.
	gmpo=0.
	end

CDECK  ID>, SCHAF.
****************** schaf **************************************

      subroutine schaf(aks,g1p,g1n,f2p,f2n,a10,au0)
c
c     g1p,g1n are polarized structure function from
c     a.schafer phys.lett. b208(1988)175
c
      implicit real*8(a-h,o-z)
      data a7/7.00d0/

      xuv =2.75033*aks**(.588d0)*(1d0-aks)**(2.69d0)
      xdv =8.52617*aks**(1.03d0)*(1d0-aks)**(6.87d0)
      a0=xuv-xdv/2d0
      a1=1.5d0*xdv
      f2p=4./9.*a0+2./9.*a1
      f2n=1./9.*a0+1./3.*a1


      fu0=1./(1.+au0*aks**(-.588)*(1.-aks)**2)
      fu1=1./(1.+a10*au0*aks**(-.588)*(1.-aks)**2)
      fd0=1./(1.+a7*au0*aks**(-1.03)*(1.-aks)**2)
      fd1=1./(1.+a7*a10*au0*aks**(-1.03)*(1.-aks)**2)
      g1p=(4./9.*a0*fu0+2./27.*(-2./3.*fd1+1./3.*fd0)*a1+
     . 4./27.*(-2./3.*fu1+1./3.*fu0)*a1)/2./aks
      g1n=(1./9.*a0*fd0+8./27.*(-2./3.*fu1+1./3.*fu0)*a1+
     . 1./27.*(-2./3.*fd1+1./3.*fd0)*a1)/2./aks

      return
      end

CDECK  ID>, SUPST.
********************** supst ************************************

      double precision function supst(t)
      implicit real*8(a-h,o-z)
      data chbar/.197328d0/
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   s.stein et al phys. rev. 12(1975)1884 (appendix 1)
c
	       sqtf=dsqrt(tf)
	       delff=(datan(sqtf/.93d0)-2.*datan(sqtf/3.19d0)+
     .	  datan(sqtf/5.45d0))*1.580d0/sqtf
	       delff=dmax1(0.d0,delff)
	       supele=1.-delff**2
	       supst=dmax1(0.d0,supele)
       return
      end

CDECK  ID>, PORTN.
********************** portn ************************************

      double precision function portn(x,x10,x20)
      implicit real*8(a-h,o-z)
      x1=min(x10,x20)
      x2=max(x10,x20)
      if(x.le.x1) portn=0.
      if(x.ge.x2) portn=1.
      if(x.gt.x1.and.x.lt.x2)portn=(x-x1)**2*(3.*x2-2.*x-x1)/(x2-x1)**3
      return
      end

CDECK  ID>, POLINO.
********************** polino ************************************

      double precision function polino(x,c,l)
      implicit real*8(a-h,o-z)
      dimension c(l)
      sa=c(l)
      do 45 ia=2,l
  45  sa=sa*x+c(l-ia+1)
      polino=sa
      return
      end

CDECK  ID>, BETA.
      function beta(z,w)
      implicit real*8(a-h,o-z)
      beta=exp(gammln(z)+gammln(w)-gammln(z+w))
      return
      end
      double precision function gammln(xx)
      implicit real*8(a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *	 -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
       x=x+one
       ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
      double precision function dgamma(xx)
      implicit real*8(a-h,o-z)
	dgamma=exp(gammln(xx))
      end

CDECK  ID>, DDEXP.
****************** ddexp **************************************

      double precision function ddexp(x)
      implicit real*8(a-h,o-z)
	ddexp=0.
	if(x.gt.-50.)ddexp=exp(x)
      return
      end

CDECK  ID>, DF2D8.
      double precision function df2d8(dq2,dx)
*:=====================================================================:
*:								       :
*:	author:    m.dueren	   last update: 06.03.1991	       :
*:				   tested: yes			       :
*:								       :
*:	arguments: dq2,dx: double prec. input xbj,q2		       :
*:		   df2h8* double prec f2  output		       :
*:								       :
*:	called by: mkf2 					       :
*:								       :
*:	action:    calculate f2 structure function of the deuteron     :
*:		   nmc fit of dis-region with 8 parameters	       :
*:		   data of nmc (1992), slac,bcdms		       :
*:								       :
*:		   parametrized with a1,bi (i=1...4) as 	       :
*:								       :
*:		   f2_dis(x,q2) ~prop.				       :
*:		     [1/b(n1,n2)*x**n1*(1-x)**n2 + 1/3*n3*(1-x)**n4 ]  :
*:		     *s(x,q2)					       :
*:			  with	x = (q2+m_a)/(2m*nu + m_b**2)	       :
*:				ni= ai+bi*s			       :
*:				s = ln(ln(q2+m_a**2)/lambda)/ln(..q2_0):
*:		   reference:					       :
*:		   the new muon collaboration			       :
*:		   nuclear physics b 371 (1992) 3-31		       :
*:=====================================================================:
c
      implicit double precision (d)
c
c
c *** d1,...,d8 = 8 param of nmc, slac, bcdms (92)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c	     values: daw=1.512(gev2), dbw=0.351(gev2)
c	     ref:  bodek et al., p.r.d20(1979)1427.
c	     see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
c
      data d1,d2,d3,d4,d5,d6,d7,d8
     :	   ,d9,d10
     :	   ,daw,dbw
c
c     f2 from nmc, slac, bcdms - data (92)
     :	   /.764053,-.171774,3.48979,.611064,.946086
     :	   ,1.06818,13.8472,-2.40967
c     resonance-region
     :	   ,.89456,.16452
     :	   ,1.512,  .351 /
c
c
      df2d8=1.d-30
      dw2 = .8803686078d0+dq2*(1.d0/dx-1.d0)
      dwmas = dsqrt(dw2)
      ddw = (dwmas-1.03d0)
c
c *** ddw = w - (resonance threshold - smearing 0f 0.05 gev)
c *** lamda(qcd) = 0.2 gev
c
      if(ddw.le.0.d0) return
c
c *** use weizmann variable for low q2
c *** values: daweiz=1.512(gev2), dbweiz=0.351(gev2)
c *** ref:    bodek et al., p.r.d20(1979)1427.
c ***	      see p.1495, eq(5.1) and table viii
c
      dq2w=dq2+dbw
      dxw=dq2w/(dq2/dx+daw)
c
      dsbar = dlog(dlog(dq2w/.04d0)) - 1.404555751d0
c
      detav1 = d1+d2*dsbar
      detav2 = d3+d4*dsbar
      detas1 = d5+d6*dsbar
      detas2 = d7+d8*dsbar
c
      dxw1=1.d0-dxw
      de1=detav1
      de2=detav2+1.d0
c
c *** supression bue to "quark sum rule"
c ***	    sup.fact.= 1 - gd**2, with gd=nucleon dipole f.f.
c *** further supression due to low w phase space volume
c
      den1=(1.d0+dq2/.71d0)
      dgd2=1.d0/den1**4
      dssum = (1.d0-dgd2)
      ds = dssum * (1.d0-dexp(-4.177d0*ddw))
c
      df2d8 =
     :	   ( .83333333333333*dgamma(de1+de2)/dgamma(de1)/dgamma(de2)
     :	   *dxw**detav1*dxw1**detav2
     :	   + .33333333333333*detas1*dxw1**detas2  ) * ds
c
c *** resonance contribution
c
      dres = 0.d0
c      dres2=0.d0
c      dres3=0.d0
c
c *** lorentzian resonance ( small fermi-smearing effect)
c ***			     gamma(fermi)=0.0883 gev
c *** gamma(d2) = sqrt ( gamma(h2)**2 + gamma(fermi)**2 )
c 1.232**2 = 1.518
c 1.520**2 = 2.310
c 1.681**2 = 2.826
c 1.232**2 * 0.15**2 = 0.0342
c 1.520**2 * 0.14**2 = 0.0453
c 1.681**2 * 0.14**2 = 0.0554
c
      if (dwmas .le. 2.3d0) then
	 dres1 = d9**2*dexp(-(dwmas-1.232d0)**2/.0053d0)
     :	      /den1**3
c      dres2 = d10**2/( (dw2-2.310d0)**2 + 0.0453d0 )
c    :	      * dgd2
c      dres3 = d11**2/( (dw2-2.826d0)**2 + 0.0554d0 )
c    :	      * dgd2
c      endif
c
c *** background under resonances
c
c    mp**2 = 0.8803686078
c    mp**2-m(pi)**2=0.8608892416
c *** dqs = momentum of one pion, decaying from the resonance, in cm
c frame
c
	 dw2m = (dwmas+.05d0)**2
	 dqs=dsqrt((dw2m+0.8608892416d0)**2/4.d0/dw2m-0.8803686078d0)
	 dbg = (d10**2*dqs )
     :	      * dexp(-0.5d0*ddw**2) / den1
	 dres=(dres1 + dbg) * dssum
      endif
c
c *** total f2 of d2
c
      df2d8 = df2d8 + dres
c
      if (df2d8.gt.0d0) return
      df2d8 = 1.d-30
      return
      end

CDECK  ID>, DF2H8.
      double precision function df2h8(dq2,dx)
*:=====================================================================:
*:								       :
*:	author:    m.dueren	   last update: 06.03.1991	       :
*:				   tested: yes			       :
*:								       :
*:	arguments: dq2,dx: double prec. input xbj,q2		       :
*:		   df2h8* double prec f2  output		       :
*:								       :
*:	called by: mkf2 					       :
*:								       :
*:	action:    calculate f2 structure function of the proton       :
*:		   nmc fit of dis-region with 8 parameters	       :
*:		   data of nmc (1992)slac,bcdms 		       :
*:								       :
*:		   parametrized with a1,bi (i=1...4) as 	       :
*:								       :
*:		   f2_dis(x,q2) ~prop.				       :
*:		     [1/b(n1,n2)*x**n1*(1-x)**n2 + 1/3*n3*(1-x)**n4 ]  :
*:		     *s(x,q2)					       :
*:			  with	x = (q2+m_a)/(2m*nu + m_b**2)	       :
*:				ni= ai+bi*s			       :
*:				s = ln(ln(q2+m_a**2)/lambda)/ln(..q2_0):
*:		   reference:					       :
*:		   the new muon collaboration			       :
*:		   nuclear physics b 371 (1992) 3-31		       :
*:=====================================================================:
c
      implicit double precision (d)
*
c *** d1,...,d8 = 8 param of nmc, slac, bcdms (92)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c	     values: daw=1.512(gev2), dbw=0.351(gev2)
c	     ref:  bodek et al., p.r.d20(1979)1427.
c	     see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
*
      data d1,d2,d3,d4,d5,d6,d7,d8
     :	   ,d9,d10,d11,d12,d13,d14
     :	   ,d15,d16
     :	   ,daw,dbw
c     f2 from nmc, slac, bcdms	data '92 (final)
     :	   /.886627,-.11191,3.3951,1.04064,1.02702,1.40335,12.4577,-
     :	   .100622
c     resonance-region:
     :	   ,.1179, .044735, .038445, .27921, 8.8228d-5, 6.2099d-5
     :	   ,1.421,1.2582
     :	   ,1.642, .376/
c
      df2h8 =1.d-30
      dw2 = .8803686078d0+dq2*(1.d0/dx-1.d0)
      dwmas = dsqrt(dw2)
      ddw = (dwmas-1.08d0)
c
c *** ddw = w - (resonance threshold)
c *** lamda(qcd) = d2: 0.20 gev
c *** lamda(qcd) = h2: 0.15 gev
c
      if(ddw.le.0.d0) return
c
c *** use weizmann variable for low q2
c *** values = d2 : daweiz=1.512(gev2), dbweiz=0.351(gev2)
c *** values = h2 : daweiz=1.642(gev2), dbweiz=0.376(gev2)
c *** ref:    bodek et al., p.r.d20(1979)1427.
c ***	      see p.1495, eq(5.1) and table viii
c
      dq2w=dq2+dbw
      dxw=dq2w/(dq2/dx+daw)
      dsbar=dlog(dlog(dq2w/.0225d0)) - 1.538942135d0
c      dsbar=dlog(dlog(dq2w/.0225d0)/dlog((2.d0+dbw)/.0225d0))
c      dsbar = dlog(dlog(dq2w/.04d0)) - 1.404555751d0
      detav1 = d1+d2*dsbar
      detav2 = d3+d4*dsbar
      detas1 = d5+d6*dsbar
      detas2 = d7+d8*dsbar
c
      dxw1=1.d0-dxw
      de1=detav1
      de2=detav2+1.d0
c
c *** supression due to "quark sum rule"
c ***	    sup.fact.= 1 - gd**2, with gd=nucleon dipole f.f.
c *** further supression due to low w phase space volume
c
      den1=(1.d0+dq2/.71d0)
      dgd2=1.d0/den1**4
      dssum = (1.d0-dgd2)
      dsthr = 1.d0
      if( ddw .le. 5.0d0) then
	 dtemp = dexp(ddw*3.98213222d0)
	 dsthr = (dtemp-1.d0)/(dtemp+1.d0)
      endif
      ds = dssum * dsthr
c
      df2h8 =
     :	   ( .83333333333333*dgamma(de1+de2)/dgamma(de1)/dgamma(de2)
     :	   *dxw**detav1*dxw1**detav2
     :	   + .33333333333333*detas1*dxw1**detas2  ) * ds
c
c *** resonance region
c
      dres = 0.d0

      if(ddw .le. 5.0d0) then
c
c *** >>> + background under the resonance
c
c *** appropriate kinematic variables
c     ... dqs = momentum of one pion in pi-nucleon c.m.-system
c		 in the case of single pion production
c *** n.b.  dqs = 0  at  w = 1.08gev (inelastic threshold)
c	 mp**2 = 0.8803686078
c	 mp**2-m(pi)**2=0.8608892416
c
	 dqs2 = (dw2+0.8608892416d0)**2/4.d0/dw2 - 0.8803686078d0
	 dqs = dsqrt(dqs2)
c
c *** >>> + resonance shape
c
c *** lorentzian resonance
c *** this includes the w**2-dependence of the res.width.
c
c *** appropriate kinematic variables
c  1) correction to res.width due to resonance threshold
c     ... dqs = momentum of one pion in pi-nucleon c.m.-system
c		 in the case of single pion production
c  2) correction to res.width due to the q2-dependence
c     ... dks = momentum of virtual photon in pi-n c.m.-system
c
	 if(ddw .le. 1.d0) then
c
	    dks2 =
     :		 (dw2+0.8803686078d0+dq2)**2/4.d0/dw2 - 0.8803686078d0
	    dks = dsqrt(dks2)
c
c *** resonance form factor (effective guess!)
c
	    dresff = 1. / den1**(d15**2)
c
c *** 1236
c     wres**2	      = 1.232**2	    = 1.518
c     (wres*gamma)**2 = 1.232**2 * 0.119**2 = 0.02149
c
	    dw2r = 1.518d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    dqsr = dsqrt(dqsr2)
	    dksr2 =
     :		 (dw2r+0.8803686078+dq2)**2/4./dw2r - 0.8803686078
	    ddcorr = (dqs/dqsr) * (1.+.16/dqsr2)/(1.+.16/dqs2)
	    dncorr = ddcorr * (dksr2+.16)/(dks2+.16)
	    ddcorr = ddcorr**2
	    dres1 = d9**2 * dncorr
     :		 /( (dw2-1.518d0)**2 + 0.02149d0*ddcorr )
c
c *** 1520
c     wres**2	      = 1.520**2	    = 2.310
c     (wres*gamma)**2 = 1.520**2 * 0.097**2 = 0.02127
c     n.b. q2-dependence of the width is neglected
c
	    dw2r = 2.310d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    ddcorr = dqs2/dqsr2
	    dres2 = d10**2 * ddcorr
     :		 / ( (dw2-2.310d0)**2 + 0.02127d0*ddcorr )
c
c *** 1681
c     wres**2	      = 1.681**2	    = 2.826
c     (wres*gamma)**2 = 1.681**2 * 0.105**2 = 0.03115
c     n.b. q2-dependence of the width is neglected
c
	    dw2r = 2.826d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    dqsr = dsqrt(dqsr2)
	    ddcorr = ( dqs/dqsr )**3
	    dres3 = d11**2 * ddcorr
     :		 / ( (dw2-2.826d0)**2 + 0.03115d0*ddcorr )
c
c
c *** sum of all resonances
c	  * resonance form factor(q2-dependence)
c
	    dres = (dres1 +dres2 +dres3)*dresff
c
c *** end of resonance calculation (only if   ddw < 1.0 gev)
c
	 endif
c
c *** background under the resonances
c     n.b. exp(-0.92**2*3.5393817) = 0.05
c	 dbgff = dq2val/dxval /den1**(dp(16)**2)
c
	 dbgff = 1. / den1**(d16**2)
	 dbg = (d12**2*dsqrt(dqs) +d13**2*dqs +d14**2*dqs2 )
     :	      * dbgff * dexp(-0.5d0*ddw**2)
c
c *** (resonance region) = ( (resonance) + (background) )
c			   * dssum(=suppression)
c
	 dres = (dres + dbg) * dssum
c
c *** end of resonance calculation
c
      endif
c
c *** (total) = (qcd part) + (resonance region)
c
      df2h8 = df2h8 + dres
c
      if(df2h8 .gt. 0.d0) return
c
      df2h8=1.d-30
      return
      end

CDECK  ID>, TITOUT.
************ titout **************************
      subroutine titout(date,bmom,tmom,pl1,pn1,qn)
      implicit real*8(a-h,o-z)
      character*10 date
      character*12 switch(100)
      dimension lsw(100)

      open(unit=16,file='tails.dat')
      open(unit=9,file='all.dat')
      open(unit=7,file='asm.dat')
      open(unit=21,file='allu.dat')
      open(unit=23,file='allp.dat')
	 write(7,1)date
	 write(9,1)date
	 write(16,1)date
	 write(16,2)
	 write(9,3)
2     format(/' the file contains information about quantities'/
     .' delta (see (42) of smc/7/93)')
3     format(/' the file contains work information about'/
     .' contribution of tails')
	 write( 9,43)
	 write(16,43)
1     format(1x,'program polrad20 version from ',a10)
	 write( 7,4)
4     format(/' the file gives born asymmetry, observed asymmetry'/
     .' and radiative correstion')

	 write( 7,43)

43    format(/' the following switches are active')
	 nsw=0
	 nsw=nsw+1
	 switch(nsw)='polrad'
	 lsw(nsw)=6
	 nsw=nsw+1
	 switch(nsw)='strffun'
	 lsw(nsw)=7
	 nsw=nsw+1
	 switch(nsw)='integrat'
	 lsw(nsw)=8
	 nsw=nsw+1
	 switch(nsw)='polrad_add'
	 lsw(nsw)=10
	 nsw=nsw+1
	 switch(nsw)='exact'
	 lsw(nsw)=5
	 nsw=nsw+1
	 switch(nsw)='elect'
	 lsw(nsw)=5
	 nsw=nsw+1
	 switch(nsw)='long'
	 lsw(nsw)=4
	 nsw=nsw+1
	 switch(nsw)='kin_hermes'
	 lsw(nsw)=10
	 nsw=nsw+1
	 switch(nsw)='targ_o'
	 lsw(nsw)=6
	 nsw=nsw+1
	 switch(nsw)='f2nmc_d8'
	 lsw(nsw)=8
	 nsw=nsw+1
	 switch(nsw)='g2_eq_0'
	 lsw(nsw)=7
	 nsw=nsw+1
	 switch(nsw)='pol_asym'
	 lsw(nsw)=8

	 isw=0
	 do while(isw.lt.nsw)
	    isw0=isw+1
	    icc=3
	    do while(icc.le.71)
	       if(isw.eq.nsw)goto 449
	       isw=isw+1
	       icc=icc+lsw(isw)+1
	    enddo
	    isw=isw-1
449	    write(7,420)(switch(iisw)(1:lsw(iisw)+1),iisw=isw0,isw)
	    write(9,420)(switch(iisw)(1:lsw(iisw)+1),iisw=isw0,isw)
	    write(16,420)(switch(iisw)(1:lsw(iisw)+1),iisw=isw0,isw)
	 enddo
420	format(3x,15a)


	 write( 7,6)
	 write( 7,7)
	 write( 7,8)
	 write( 9,6)
	 write( 9,7)
	 write( 9,8)
	 write(16,6)
	 write(16,7)
	 write(16,8)
6     format(/' leptons are electrons')
7     format(' target is oxygen')
8     format(' target is longitudinally polarized')
      write(7,9)bmom,tmom,pl1,pn1,qn
9     format(' bmom = ',f5.1/' tmom = ',f5.1
     .	   /'    pl = ',f4.2,'    pn = ',f4.2,'   qn = ',f4.2)
      write( 7,11)
11    format(/' a is in %')
      write(9,9)bmom,tmom,pl1,pn1,qn
      write(16,9)bmom,tmom,pl1,pn1,qn
      write(16,10)
      write(16,12)
10    format(/' d is in %')
12    format(/'  x',6x,'y',3x,' d _i_p ',' d_e_p ',' d_q_p ',
     .' d_v_p ',' d_i_u ',' d_e_u ',' d_q_u ',' d_v_u ')
      write( 7,13)
13    format(/'      x',6x,'w2',7x,'q2',2x
     .,'a(born) ',' a(obs) ',' del(%) ')
      return
      end

CDECK  ID>, INTEGRAT.
c     there are modules necessary for integration over the
c     phase space of the emitted photon and over the variable
c     y in this patch. integrators dqunc8, dqvnc8, dqwnc8 are
c     identical and based on the newton-cotes method of the
c     eighth order.

CDECK  ID>, DQUNC8.

      subroutine dqunc8(fun,a,b,ab,rl,r,er,nn,fl,nx)


      implicit real*8 (a-h,o-z)
      dimension h(31),f(16),v(8,30),z(8,30),x(16)

      parameter (l=6)
      parameter (wd = 5.0d-1)
      parameter (wn = 2.79082892416225747d-1)
      parameter (w1 = 1.66151675485008798d+0)
      parameter (w2 =-2.61869488536155201d-1)
      parameter (w3 = 2.96183421516754830d+0)
      parameter (w4 =-1.28112874779541430d+0)

c+seq, automat.

      lm=30
      n=nx-1216
      if(n.lt.200)n=200
      r=0d0
      fl=r
      er=fl
      nn=0
      if(a.eq.b)return
      c=fl
      ar=fl


      k=nn
      m=1
      xn=a
      x(16)=b
      p=fl
      fn=fun(xn)
      s=(b-a)*625d-4
      x(8)=(a+b)*wd
      x(4)=(a+x(8))*wd
      x(12)=(x(8)+b)*wd
      x(2)=(a+x(4))*wd
      x(6)=(a+x(12))*wd
      x(10)=(x(4)+b)*wd
      x(14)=(x(12)+b)*wd
      do 25 j=2,16,2
   25 f(j)=fun(x(j))
      nn=9


   30 x(1)=(xn+x(2))*wd
      f(1)=fun(x(1))
      do 35 j=3,15,2
      x(j)=(x(j-1)+x(j+1))*wd
   35 f(j)=fun(x(j))
      nn=nn+8
      st=(x(16)-xn)*625d-4
      q=((fn+f(8))*wn+(f(1)+f(7))*w1+
     *(f(2)+f(6))*w2+(f(3)+f(5))*w3+
     *f(4)*w4)*st
      h(k+1)=((f(8)+f(16))*wn+(f(9)+f(15))*w1+
     *(f(10)+f(14))*w2+(f(11)+f(13))*w3+
     *f(12)*w4)*st
      w=q+h(k+1)
      d=w-p
      ar=ar+d


      e=dabs(d)/1023d0
      t=dmax1(ab,rl*dabs(ar))*(st/s)
      if(k.lt.1)go to 50
      if(k.ge.lm)go to 62
      if(nn.gt.n)go to 60
      if(e.le.t)go to 70


   50 m=2*m
      k=k+1
      do 52 i=1,8
      j=i+8
      v(i,k)=f(j)
   52 z(i,k)=x(j)
      p=q
      do 55 i=1,8
      j=9-i
      f(2*j)=f(j)
   55 x(2*j)=x(j)
      go to 30


   60 n=2*n
      lm=l
      fl=fl+(b-xn)/(b-a)
      go to 70


   62 fl=fl+1d0


   70 r=r+w
      er=er+e
      c=c+d/1023d0
   72 if(m.eq.2*(m/2))go to 75
      m=m/2
      k=k-1
      go to 72
   75 m=m+1
      if(k.le.0)go to 80
      p=h(k)
      xn=x(16)
      fn=f(16)
      do 78 i=1,8
      f(2*i)=v(i,k)
   78 x(2*i)=z(i,k)
      go to 30


   80 r=r+c
      if(er.eq.0d0)return
      ar=dabs(r)


   82 q=ar+er
      if(q.ne.ar)return
      er=2d0*er
      go to 82
      end

ccc-----------------------------------------------------------------------------




CDECK  ID>, DQVNC8.

      subroutine dqvnc8(fun,a,b,ab,rl,r,er,nn,fl,nx)


      implicit real*8 (a-h,o-z)
      dimension h(31),f(16),v(8,30),z(8,30),x(16)

      parameter (l=6)
      parameter (wd = 5.0d-1)
      parameter (wn = 2.79082892416225747d-1)
      parameter (w1 = 1.66151675485008798d+0)
      parameter (w2 =-2.61869488536155201d-1)
      parameter (w3 = 2.96183421516754830d+0)
      parameter (w4 =-1.28112874779541430d+0)

c+seq, automat.

      lm=30
      n=nx-1216
      if(n.lt.200)n=200
      r=0d0
      fl=r
      er=fl
      nn=0
      if(a.eq.b)return
      c=fl
      ar=fl


      k=nn
      m=1
      xn=a
      x(16)=b
      p=fl
      fn=fun(xn)
      s=(b-a)*625d-4
      x(8)=(a+b)*wd
      x(4)=(a+x(8))*wd
      x(12)=(x(8)+b)*wd
      x(2)=(a+x(4))*wd
      x(6)=(a+x(12))*wd
      x(10)=(x(4)+b)*wd
      x(14)=(x(12)+b)*wd
      do 25 j=2,16,2
   25 f(j)=fun(x(j))
      nn=9


   30 x(1)=(xn+x(2))*wd
      f(1)=fun(x(1))
      do 35 j=3,15,2
      x(j)=(x(j-1)+x(j+1))*wd
   35 f(j)=fun(x(j))
      nn=nn+8
      st=(x(16)-xn)*625d-4
      q=((fn+f(8))*wn+(f(1)+f(7))*w1+
     *(f(2)+f(6))*w2+(f(3)+f(5))*w3+
     *f(4)*w4)*st
      h(k+1)=((f(8)+f(16))*wn+(f(9)+f(15))*w1+
     *(f(10)+f(14))*w2+(f(11)+f(13))*w3+
     *f(12)*w4)*st
      w=q+h(k+1)
      d=w-p
      ar=ar+d


      e=dabs(d)/1023d0
      t=dmax1(ab,rl*dabs(ar))*(st/s)
      if(k.lt.1)go to 50
      if(k.ge.lm)go to 62
      if(nn.gt.n)go to 60
      if(e.le.t)go to 70


   50 m=2*m
      k=k+1
      do 52 i=1,8
      j=i+8
      v(i,k)=f(j)
   52 z(i,k)=x(j)
      p=q
      do 55 i=1,8
      j=9-i
      f(2*j)=f(j)
   55 x(2*j)=x(j)
      go to 30


   60 n=2*n
      lm=l
      fl=fl+(b-xn)/(b-a)
      go to 70


   62 fl=fl+1d0


   70 r=r+w
      er=er+e
      c=c+d/1023d0
   72 if(m.eq.2*(m/2))go to 75
      m=m/2
      k=k-1
      go to 72
   75 m=m+1
      if(k.le.0)go to 80
      p=h(k)
      xn=x(16)
      fn=f(16)
      do 78 i=1,8
      f(2*i)=v(i,k)
   78 x(2*i)=z(i,k)
      go to 30


   80 r=r+c
      if(er.eq.0d0)return
      ar=dabs(r)


   82 q=ar+er
      if(q.ne.ar)return
      er=2d0*er
      go to 82
      end

ccc-----------------------------------------------------------------------------




CDECK  ID>, D01FCE.
      subroutine d01fce(ndim, a, b, minpts, maxpts, functn, eps,
     * acc, lenwrk, wrkstr, finval, ifail)
      implicit real*8(a-h,o-z)
c     mark 8 release. nag copyright 1979.
c
c     adaptive multidimensional integration subroutine
c
c     *********  parameters for d01fce ****************************
c
c      input parameters
c
c     ndim    integer number of variables, must exceed 1 but
c	  not exceed 15.
c
c     a       real array of lower limits, with dimension ndim
c
c     b       real array of upper limits, with dimension ndim
c
c     minpts  integer minimum number of integrand values to be
c	  allowed, which must not exceed maxpts.
c
c     maxpts  integer maximum number of integrand values to be
c	  allowed, which must be at least
c	  2**ndim+2*ndim**2+2*ndim+1.
c
c     functn  externally declared user defined real function
c	  integrand. it must have parameters (ndim,z),
c	  where z is a real array of dimension ndim.
c
c     eps     real required relative accuracy, must be greater
c	  than zero
c
c     lenwrk  integer length of array wrkstr, must be at least
c	  2*ndim+4.
c
c     ifail   integer nag failure parameter
c	  ifail=0 for hard fail
c	  ifail=1 for soft fail
c
c      output parameters
c
c     minpts  integer number of integrand values used by the
c	  routine
c
c     wrkstr  real array of working storage of dimension (lenwrk).
c
c     acc     real estimated relative accuracy of finval
c
c     finval  real estimated value of integral
c
c     ifail   ifail=0 for normal exit, when estimated relative
c	    less integaccuracy rand values used.
c
c      ifail=1 if ndim.lt.2, ndim.gt.15, minpts.gt.maxpts,
c	    maxpts.lt.2**ndim+2*ndim*(ndim+1)+1, eps.le.0
c	    or lenwrk.lt.2*ndim+4.
c
c      ifail=2 if maxpts was too small for d01fce to obtain the
c	    required relative accuracy eps.  in this
c	    case d01fce returns a value of finval
c	    with estimated relative accuracy acc.
c
c      ifail=3 if lenwrk too small for maxpts integrand
c	    values.  in this case d01fce returns a
c	    value of finval with estimated accuracy
c	    acc using the working storage
c	    available, but acc will be greater
c	    than eps.
c
c     **************************************************************
c
c     .. scalar arguments ..
*     real eps, finval, acc
*     integer ifail, lenwrk, maxpts, minpts, ndim
c     .. array arguments ..
**    real a(ndim), b(ndim), wrkstr(lenwrk)
      dimension a(ndim), b(ndim), wrkstr(lenwrk)
c     .. function arguments ..
*     real functn
c     ..
c     .. local scalars ..
      character*8 srname
      double precision
     * abserr, df1, df2, difmax, f1, f2, f3, f4, half, lamda2,
     * lamda4, lamda5, one, ratio, rgncmp, rgnerr, rgnert, rgnval,
     * rgnvlt, rgnvol, rlndim, sum1, sum2, sum3, sum4, sum5, two,
     * twondm, weit1, weit2, weit3, weit4, weit5, weitp1, weitp2,
     * weitp3, weitp4, zero
      integer dvaxes, dvaxis, dvflag, funcls, ierror, j, k, maxaxs,
     * mxrgns, pointr, rgncls, rulcls, sbrgns, subrgn, subtmp,
     * tpontp, tpontr
c     .. local arrays ..
      dimension center(15), dif(15), oldcnt(15), width(15), z(15)
      integer dvcntl(15), dvcntr(15)
c     .. function references ..
*     real sqrt, x02aae
      integer p01aae, x02bbe
c     ..
      data srname /'  d01fce'/
      data zero, one, two, half /0.d0, 1.d0, 2.d0, 0.5d0/
c
c   subroutine initialisation and parameter checking
c
      if (ndim.lt.2 .or. ndim.gt.15) go to 560
      if (minpts.gt.maxpts) go to 560
      if (eps.le.zero) go to 560
      if (lenwrk.lt.2*ndim+4) go to 560
      funcls = 0
      finval = zero
      abserr = zero
      twondm = two**ndim
      rgnvol = twondm
      dvflag = 1
      fffff1 = float(x02bbe(one))
      fffff2 = 1.0/x02aae(0.0d0)
      maxaxs = int(dmin1(fffff1,fffff2))
c     maxaxs = int(amin1(float(x02bbe(one)),1.0/x02aae(0.0d0)))
      maxaxs = (maxaxs-ndim)/(ndim+1)
      mxrgns = lenwrk/(2*ndim+4)
      sbrgns = 0
      rgnvlt = zero
      rgnert = zero
      do 20 j=1,ndim
       center(j) = (a(j)+b(j))*half
       dif(j) = zero
       width(j) = (b(j)-a(j))*half
       dvcntl(j) = 1
       dvcntr(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   20 continue
c
c   end subroutine initialisation
c   basic rule initialisation
c
      rulcls = 2**ndim + 2*ndim*ndim + 2*ndim + 1
      funcls = rulcls
      if (maxpts.lt.rulcls) go to 560
      rlndim = ndim
      lamda2 = sqrt(9.0/70.0)
      lamda4 = sqrt(9.0/10.0)
      lamda5 = sqrt(9.0/19.0)
      weit1 = (12824.0-9120.0*rlndim+400.0*rlndim*rlndim)/19683.0
      weit2 = 980.0/6561.0
      weit3 = (1820.0-400.0*rlndim)/19683.0
      weit4 = 200.0/19683.0
      weit5 = 6859.0/19683.0/twondm
      weitp1 = (729.0-950.0*rlndim+50.0*rlndim**2)/729.0
      weitp2 = 245.0/486.0
      weitp3 = (265.0-100.0*rlndim)/1458.0
      weitp4 = 25.0/729.0
      ratio = (lamda2/lamda4)**2
c
c   end basic rule initialisation
      go to 100
c   divide subregion with largest error and prepare to use
c   basic rule on each portion
c
   40 subrgn = 1
      pointr = wrkstr(1)
      rgncls = rulcls
      rgnvol = twondm
      tpontr = pointr + 2
      do 60 j=1,ndim
       tpontr = tpontr + 2
       center(j) = wrkstr(tpontr-1)
       width(j) = wrkstr(tpontr)
       dvcntr(j) = 1
       dvcntl(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   60 continue
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.lt.0) go to 600
   80 dvaxis = dvaxes
      dvaxes = dvaxis/(ndim+1)
      dvaxis = dvaxis - (ndim+1)*dvaxes
      dvcntl(dvaxis) = 2*dvcntl(dvaxis)
      rgncls = rgncls*2
      if (dvaxes.gt.0) go to 80
      if (funcls+rgncls.gt.maxpts) go to 580
      if (rgncls/rulcls+sbrgns-1.gt.mxrgns) dvflag = 2
      funcls = funcls + rgncls
c      print *,funcls
      abserr = abserr - wrkstr(pointr)
      finval = finval - wrkstr(pointr+1)
c
c   begin basic rule
  100 do 120 j=1,ndim
       z(j) = center(j)
  120 continue
      sum1 = functn(ndim,z)
      sum2 = zero
      sum3 = zero
      do 140 j=1,ndim
       z(j) = center(j) - lamda2*width(j)
       f1 = functn(ndim,z)
       z(j) = center(j) + lamda2*width(j)
       f2 = functn(ndim,z)
       z(j) = center(j) - lamda4*width(j)
       f3 = functn(ndim,z)
       z(j) = center(j) + lamda4*width(j)
       f4 = functn(ndim,z)
       sum2 = sum2 + f1 + f2
       sum3 = sum3 + f3 + f4
       df1 = f1 + f2 - two*sum1
       df2 = f3 + f4 - two*sum1
       dif(j) = dif(j) + abs(df1-ratio*df2)
       z(j) = center(j)
  140 continue
      sum4 = zero
      do 200 j=2,ndim
       z(j-1) = center(j-1) - lamda4*width(j-1)
       do 160 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  160  continue
       z(j-1) = center(j-1) + lamda4*width(j-1)
       do 180 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  180  continue
       z(j-1) = center(j-1)
  200 continue
      sum5 = zero
      do 220 j=1,ndim
       z(j) = center(j) - lamda5*width(j)
  220 continue
  240 do 260 j=2,ndim
       if (z(j-1).lt.center(j-1)+width(j-1)) go to 280
       z(j-1) = center(j-1) - lamda5*width(j-1)
       z(j) = z(j) + two*lamda5*width(j)
  260 continue
      if (z(ndim).gt.center(ndim)+width(ndim)) go to 300
  280 sum5 = sum5 + functn(ndim,z)
      z(1) = z(1) + two*lamda5*width(1)
      go to 240
  300 rgnval = rgnvol*(weit1*sum1+weit2*sum2+weit3*sum3+weit4*
     * sum4+weit5*sum5)
      rgncmp = rgnvol*(weitp1*sum1+weitp2*sum2+weitp3*sum3+weitp4*
     * sum4)
      rgnerr = abs(rgnval-rgncmp)
c
c   end basic rule
c   store results of basic rule application
c
      rgnvlt = rgnvlt + rgnval
      rgnert = rgnert + rgnerr
      finval = finval + rgnval
      abserr = abserr + rgnerr
      if (dvflag.eq.0) go to 340
      if (dvflag.eq.2) go to 500
      pointr = mxrgns + sbrgns*(2*ndim+3) + 1
      sbrgns = sbrgns + 1
      wrkstr(sbrgns) = pointr
      subrgn = sbrgns
      tpontr = pointr + 2
      do 320 j=1,ndim
       tpontr = tpontr + 2
       wrkstr(tpontr-1) = center(j)
       wrkstr(tpontr) = width(j)
  320 continue
  340 wrkstr(pointr) = rgnert
      wrkstr(pointr+1) = rgnvlt
c   determine axis along which fourth difference is largest
      difmax = zero
      do 380 j=1,ndim
       if (difmax.gt.dif(j)) go to 360
       difmax = dif(j)
       dvaxis = j
  360	    dif(j) = zero
  380 continue
      tpontr = pointr + 2*(dvaxis+1)
      wrkstr(tpontr) = width(dvaxis)*half
      wrkstr(tpontr-1) = center(dvaxis) - wrkstr(tpontr)
      if (dvflag.ne.2) go to 400
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.gt.maxaxs) dvaxes = -1
      dvaxis = dvaxis + (ndim+1)*dvaxes
  400 wrkstr(pointr+2) = dvaxis
      if (dvflag.eq.1) go to 460
c   determine the position in the parially ordered list of
c   the subregion which replaces most recently divided subregion
  420 subtmp = 2*subrgn
      if (subtmp.gt.sbrgns) go to 480
      tpontr = wrkstr(subtmp)
      if (subtmp.eq.sbrgns) go to 440
      tpontp = wrkstr(subtmp+1)
      if (wrkstr(tpontr).ge.wrkstr(tpontp)) go to 440
      subtmp = subtmp + 1
      tpontr = tpontp
  440 if (rgnert.ge.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 420
c   when working storage is not used up, determine the
c   position in the partially ordered list for the description
c   of other portion(s) of most recently divided subregion
  460 subtmp = subrgn/2
      if (subtmp.lt.1) go to 480
      tpontr = wrkstr(subtmp)
      if (rgnert.le.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 460
  480 rgnvlt = zero
      rgnert = zero
      if (dvflag.eq.2) go to 540
      dvflag = 1 - dvflag
c   count to determine the next part of the recently divided
c   subregion for application of the basic rule
  500 center(1) = center(1) + two*width(1)
      dvcntr(1) = dvcntr(1) + 1
      do 520 j=2,ndim
       if (dvcntr(j-1).le.dvcntl(j-1)) go to 100
       dvcntr(j-1) = 1
       center(j-1) = oldcnt(j-1)
       dvcntr(j) = dvcntr(j) + 1
       center(j) = center(j) + two*width(j)
  520 continue
      if (dvcntr(ndim).le.dvcntl(ndim)) go to 100
      center(ndim) = oldcnt(ndim)
      if (dvflag.eq.2) go to 340
c
c   end ordering of basic rule results
c   make checks for possible termination of routine
c
  540 acc = abserr/abs(finval)
      if (acc.gt.eps .or. funcls.lt.minpts) go to 40
c
c   loop back to apply basic rule
c
c   termination point, set ifail and return
c
      ierror = 0
      go to 620
  560 ierror = 1
      go to 620
  580 ierror = 2
      go to 620
  600 ierror = 3
  620 minpts = funcls
      ifail = p01aae(ifail,ierror,srname)
      return
      end

      double precision function x02aae(x)
      implicit real*8(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
c+self,if=ibm.
cc     for ibm/360/370/3090
c      data z/z3380000000000000/
c      x02aae = z
c     for sun
      data z/1.1d-16/
      x02aae = z
c     * eps *
c     returns the value eps where eps is the smallest
c     positive
c     number such that 1.0 + eps > 1.0
c     the x parameter is not used
c     for icl 1900
c     x02aae = 2.0**(-37.0)
c+self,if=pc.
c     for pdp11
c      x02aae=2.d0**(-23.d0)
c+self.

      return
      end
c
      integer  function x02bbe(x)
      implicit real*8(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
*     real x
c     * maxint *
c     returns the largest integer representable on the computer
c     the x parameter is not used
c     for icl 1900
c      x02bbe = 8388607
c     for ibm,sun,vax,ibm pc/386/486
       x02bbe = 2147483647
c   for pdp11
c     x02bbe=32767
      return
      end

      integer function p01aae(ifail, error, srname)
c     mark 1 release.  nag copyright 1971
c     mark 3 revised
c     mark 4a revised, ier-45
c     mark 4.5 revised
c     mark 7 revised (dec 1978)
c     returns the value of error or terminates the program.
      integer error, ifail, nout
      character*8 srname
c     test if no error detected
      if (error.eq.0) go to 20
c     determine output unit for message
      call x04aae (0,nout)
c     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
c     hard failure
      write (nout,99999) srname, error
c     stopping mechanism may also differ
      stop
c     soft fail
c     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nout,99999) srname, error
   20 p01aae = error
      return
99999 format (1h0, 38herror detected by nag library routine , a8,
     * 11h - ifail = , i5//)
      end
      subroutine x04aae(i,nerr)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     *** note ***
c     this routine assumes that the value of nerr1 is saved
c     between calls.  in some implementations it may be
c     necessary to store nerr1 in a labelled common
c     block /ax04aa/ to achieve this.
c
c     .. scalar arguments ..
      integer i, nerr
c     ..
c     .. local scalars ..
      integer nerr1
c     ..
      data nerr1 /6/
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end

CDECK  ID>, DQG32.
      subroutine dqg32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
c
CDECK  ID>, QUNC8.
* * * * * * * * * * * * * * * * * * * * * * * * * * *
*						    *
      subroutine qunc8(fun,a,b,ab,rl,r,er,nn,fl,nx)
*						    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *

      implicit real*8 (a-h,o-z)
      dimension h(31),f(16),v(8,30),z(8,30),x(16)
      integer*4 m

      parameter (l=6)
      parameter (wd = 5.0d-1)
      parameter (wn = 2.79082892416225747d-1)
      parameter (w1 = 1.66151675485008798d+0)
      parameter (w2 =-2.61869488536155201d-1)
      parameter (w3 = 2.96183421516754830d+0)
      parameter (w4 =-1.28112874779541430d+0)


      lm=30
      n=nx-1216
      if(n.lt.200)n=200
      r=0d0
      fl=r
      er=fl
      nn=0
      if(a.eq.b)return
      c=fl
      ar=fl


      k=nn
      m=1
      xn=a
      x(16)=b
      p=fl
      fn=fun(xn)
      s=(b-a)*625d-4
      x(8)=(a+b)*wd
      x(4)=(a+x(8))*wd
      x(12)=(x(8)+b)*wd
      x(2)=(a+x(4))*wd
      x(6)=(a+x(12))*wd
      x(10)=(x(4)+b)*wd
      x(14)=(x(12)+b)*wd
      do 25 j=2,16,2
   25 f(j)=fun(x(j))
      nn=9


   30 x(1)=(xn+x(2))*wd
      f(1)=fun(x(1))
      do 35 j=3,15,2
      x(j)=(x(j-1)+x(j+1))*wd
   35 f(j)=fun(x(j))
      nn=nn+8
      st=(x(16)-xn)*625d-4
      q=((fn+f(8))*wn+(f(1)+f(7))*w1+
     *(f(2)+f(6))*w2+(f(3)+f(5))*w3+
     *f(4)*w4)*st
      h(k+1)=((f(8)+f(16))*wn+(f(9)+f(15))*w1+
     *(f(10)+f(14))*w2+(f(11)+f(13))*w3+
     *f(12)*w4)*st
      w=q+h(k+1)
      d=w-p
      ar=ar+d


      e=dabs(d)/1023d0
      t=dmax1(ab,rl*dabs(ar))*(st/s)
      if(k.lt.1)go to 50
      if(k.ge.lm)go to 62
      if(nn.gt.n)go to 60
      if(e.le.t)go to 70


   50 m=2*m
      k=k+1
      do 52 i=1,8
      j=i+8
      v(i,k)=f(j)
   52 z(i,k)=x(j)
      p=q
      do 55 i=1,8
      j=9-i
      f(2*j)=f(j)
   55 x(2*j)=x(j)
      go to 30


   60 n=2*n
      lm=l
      fl=fl+(b-xn)/(b-a)
      go to 70


   62 fl=fl+1d0


   70 r=r+w
      er=er+e
      c=c+d/1023d0
   72 if(m.eq.2*(m/2))go to 75
      m=m/2
      k=k-1
      go to 72
   75 m=m+1
      if(k.le.0)go to 80
      p=h(k)
      xn=x(16)
      fn=f(16)
      do 78 i=1,8
      f(2*i)=v(i,k)
   78 x(2*i)=z(i,k)
      go to 30


   80 r=r+c
      if(er.eq.0d0)return
      ar=dabs(r)


   82 q=ar+er
      if(q.ne.ar)return
      er=2d0*er
      go to 82
      end

CDECK  ID>, DQN32.
      subroutine dqn32(xl,xu,fct,y)
c
c
      double precision xl,xu,y,a,b,c,fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
c
CDECK  ID>, QVNC8.
* * * * * * * * * * * * * * * * * * * * * * * * * * *
*						    *
      subroutine qvnc8(fun,a,b,ab,rl,r,er,nn,fl,nx)
*						    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *

      implicit real*8 (a-h,o-z)
      dimension h(31),f(16),v(8,30),z(8,30),x(16)
      integer*4 m

      parameter (l=6)
      parameter (wd = 5.0d-1)
      parameter (wn = 2.79082892416225747d-1)
      parameter (w1 = 1.66151675485008798d+0)
      parameter (w2 =-2.61869488536155201d-1)
      parameter (w3 = 2.96183421516754830d+0)
      parameter (w4 =-1.28112874779541430d+0)


      lm=30
      n=nx-1216
      if(n.lt.200)n=200
      r=0d0
      fl=r
      er=fl
      nn=0
      if(a.eq.b)return
      c=fl
      ar=fl


      k=nn
      m=1
      xn=a
      x(16)=b
      p=fl
      fn=fun(xn)
      s=(b-a)*625d-4
      x(8)=(a+b)*wd
      x(4)=(a+x(8))*wd
      x(12)=(x(8)+b)*wd
      x(2)=(a+x(4))*wd
      x(6)=(a+x(12))*wd
      x(10)=(x(4)+b)*wd
      x(14)=(x(12)+b)*wd
      do 25 j=2,16,2
   25 f(j)=fun(x(j))
      nn=9


   30 x(1)=(xn+x(2))*wd
      f(1)=fun(x(1))
      do 35 j=3,15,2
      x(j)=(x(j-1)+x(j+1))*wd
   35 f(j)=fun(x(j))
      nn=nn+8
      st=(x(16)-xn)*625d-4
      q=((fn+f(8))*wn+(f(1)+f(7))*w1+
     *(f(2)+f(6))*w2+(f(3)+f(5))*w3+
     *f(4)*w4)*st
      h(k+1)=((f(8)+f(16))*wn+(f(9)+f(15))*w1+
     *(f(10)+f(14))*w2+(f(11)+f(13))*w3+
     *f(12)*w4)*st
      w=q+h(k+1)
      d=w-p
      ar=ar+d


      e=dabs(d)/1023d0
      t=dmax1(ab,rl*dabs(ar))*(st/s)
      if(k.lt.1)go to 50
      if(k.ge.lm)go to 62
      if(nn.gt.n)go to 60
      if(e.le.t)go to 70


   50 m=2*m
      k=k+1
      do 52 i=1,8
      j=i+8
      v(i,k)=f(j)
   52 z(i,k)=x(j)
      p=q
      do 55 i=1,8
      j=9-i
      f(2*j)=f(j)
   55 x(2*j)=x(j)
      go to 30


   60 n=2*n
      lm=l
      fl=fl+(b-xn)/(b-a)
      go to 70


   62 fl=fl+1d0


   70 r=r+w
      er=er+e
      c=c+d/1023d0
   72 if(m.eq.2*(m/2))go to 75
      m=m/2
      k=k-1
      go to 72
   75 m=m+1
      if(k.le.0)go to 80
      p=h(k)
      xn=x(16)
      fn=f(16)
      do 78 i=1,8
      f(2*i)=v(i,k)
   78 x(2*i)=z(i,k)
      go to 30


   80 r=r+c
      if(er.eq.0d0)return
      ar=dabs(r)


   82 q=ar+er
      if(q.ne.ar)return
      er=2d0*er
      go to 82
      end

CDECK  ID>, SIMPS.
      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpxx(a,b,np,ep,func,res)
      implicit real*8 (a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

CDECK  ID>, SIMPDO.
      subroutine simpdo(ain,afi,ep1,ii1,bin,bfi,ep2,ii2,fun,ai)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      common/simpc/x
      external fun			   ! aku
      a1=ain				   ! aku
      b1=afi				   ! aku
      h1=(b1-a1)/ii1			   ! aku
      reps1=ep1 			   ! aku
      aeps1=1d-18			   ! aku
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      call simpxx(bin,bfi,ii2,ep2,fun,functx) ! aku
*aku	  f(1)=funct(x)/3.			! aku
      f(1)=functx/3.				! aku
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
*aku   11 f(k)=funct(x)/3.
   11 continue
      call simpxx(bin,bfi,ii2,ep2,fun,functx) ! aku
*aku	  f(k)=funct(x)/3.			! aku
      f(k)=functx/3.				! aku
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

CDECK  ID>, STRFP2.
************** strfp2 *******************************************

      subroutine strfp2(xi,t,f,ifu)
      implicit real*8(a-h,o-z)
      dimension f(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
c  ich=0  for QED born cross section
c  ich=1  for electroweak (without pure QED) born cross section
c  ich=2  for complete electroweak born cross section
      do 4 ik=1,4
4     f(ik)=0d0
      end


CDECK  ID>, SIGMAB.
****************** sigmab *************************************
      double precision function sigmab(z1,z2)
      implicit real*8(a-h,o-z)
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
       yst=1.-(1.-ys)/(z1*z2)
       xst=xs*ys/yst/z2
       st=s*z1
       sigmab=ys/(yst*z1*z2**2)*boursc(xst,yst,st)

      end

CDECK  ID>, BOURSC.
************** boursc **********************************************

       double precision function boursc(xs,ys,s)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/p/pi,pi2,alfa,i1(8),i2(8)
       dimension tb(8),f(8)
       q2=s*xs*ys
       call strfp2(xs,q2,f,0)
       dmu=amp2*xs/s
       sxm2=1./ys-1.-dmu
      tb(1)=xs*ys
      tb(2)=sxm2
      tb(3)=-xs*ys*(1. + 2.*sxm2)
      tb(4)=4.*dmu*xs
      bour=0.
      do isf=isf1,isf2
c	 if(isf.eq.1.or.isf.eq.2)polst=un
c	 if(isf.eq.3.or.isf.eq.4)polst=pn
c	 if(isf.ge.5)polst=qn/3.
c	 bour=bour+polst*tb(isf)*f(isf)
	 bour=bour+tb(isf)*f(isf)
      enddo
      boursc=bour* 4.*alfa**2*pi/xs/q2*barn
      end


