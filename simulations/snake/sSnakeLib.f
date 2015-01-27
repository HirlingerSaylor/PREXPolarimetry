				subroutine descin
c**************************************************************
c                                                             *
c                          d e s c i n                        *
c                                                             *
c**************************************************************
c				read the description file
c				fill lhe commons /relfra/ , /fidef/ and /alter/
c					x_mag_mon: particle magnetic momentum in unit of ehbar/2m
c					e_over_m: particle cyclotron pulsation in rad s-1 T-1
c					twice_spin: 1. for spin 1/2 particle, 2. for spin 1 particle...
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2,rep*1
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				character*80 line*80,ch*1
				logical cyl,lfil
				character*80 fifi
				character*80 file(2)
				real factor(2)
				common /fidef/method(maxre),itype(maxre),indfi(maxre)
     @								,fact(maxre),cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				common/alter/jflag,prec,hstmax,hstmin,sinc
				logical lsynchro,lspin
				common/logic/lsynchro,lspin
				common/particle/z_particle,xm
				common/spin/x_mag_mon,e_over_m,twice_spin,xmref
				dimension rel(3)

c				min,max and increment of the step change
c					and mass of the particle (in GeV) for average synchrotron
c					radiation losses:
				read(2,*)hstmin,hstmax,sinc
				read(2,*)lspin,lsynchro
				read(2,*)z_particle,xm
				read(2,*)x_mag_mon,e_over_m,twice_spin,xmref
				indre=0
				iepz=0

c**************************new region:*************************
c				define the magnetic field in this region
200			continue
				read(2,'(a)',end=300)line
				if(line(1:4).eq.'end*')goto300
				indre=indre+1
				if(indre.gt.maxre)stop'snake:max.number of regions surpassed'
				title(indre)=line
c				write(6,*)title(indre)

c				name of this new region :
				read(2,'(a)')rname(indre)

c				position option of new region
				read(2,'(a)')rnamer(indre)

c				global coord. of local box origin:
				read(2,*)(xyz0(indre,i),i=1,3)

c				rotation of local box
				read(2,*)zangd,xangd,yangd
				zang(indre)=zangd*dtor
				xang(indre)=xangd*dtor
				yang(indre)=yangd*dtor
c				write(6,*)'    ',rname(indre),' refer.= ',rnamer(indre)
				if(rnamer(indre).ne.'absolute')then
					call itrans
				endif

c				size of box
				read(2,*)(frbox(1,indre,i),i=1,3)
				read(2,*)(frbox(2,indre,i),i=1,3)
				do 1 i=1,3
					if(frbox(1,indre,i).ge.frbox(2,indre,i))then
						write(6,*)' descin:zero or negative box volume'
						stop
					endif
1				continue

c				method of transport: 1 -> no field, 2 -> matrix, 3 -> rung-khuta
				read(2,*)method(indre)

c				type of matrix/field:
c					1 -> 1st order,constant
c					2 -> 2nd order,analytic
c					3 -> read on a file
				read(2,*)itype(indre)

c				index of the matrix/field in this type:
				read(2,*)indfi(indre)
				if((itype(indre).eq.3).and.(indfi(indre).eq.7))then
					read(2,*)file(1),factor(1)
					read(2,*)file(2),factor(2)
					read(2,*)fifi(indre)
					call mapadd(file,factor,fifi(indre))
					go to 747
				endif 

c				field file name:
				read(2,'(a)')line
				if(line(1:1).eq.'#')then
					fifi(indre)=line
				else
					i1=1
					i2=i1
					call nextcoma(i2,line)
					read(line(i1:i2-1),'(a)')fifi(indre)
				endif

c				cartesian or cylindrical coors.syst for box and map?
747			read(2,'(a)')rep
				cyl(indre)=rep.eq.'c'
				if(cyl(indre))then
					if(frbox(1,indre,1).le.0.)then
						write(6,*)' descin: zero or negative minimum',
     @								' radius for a cylindrical region'
						stop
					endif
					frbox(1,indre,2)=frbox(1,indre,2)*dtor
					frbox(2,indre,2)=frbox(2,indre,2)*dtor
				endif

c				multiplicative factor of the field(result in tesla) if method=3
c					or reference momentum if method=2:
				read(2,*)fact(indre)

c				number of additional parameters
				read(2,*)nad
				if(nad.gt.maxdat)stop ' descin: nb.of additional data >maxdat'

c				additional data
				do  8 i=1,maxdat
					adata(i,indre)=0.
					if(i.le.nad)read(2,*)adata(i,indre)
8				continue

c				define the end-planes of this mew region(yep refers to relat.coor.
				do 421 iep=1,maxepr
421			cerfiln(iep,indre)='none'
				read(2,'(a)')epsw(indre)
				if( epsw(indre).eq.'loop')then
					read(2,*)yepmin(indre),yepmax(indre),yepstp(indre)
					nep(indre)=((yepmax(indre)-yepmin(indre))/yepstp(indre))+1
					if(nep(indre).gt.maxepr)stop
     @			' descin:too many end-planes in this region !!!'
				else
					nep(indre)=0
					do 201 i=1,maxepr+1
						read(2,'(a)')line
						if(line(1:3).eq.'eol')goto202
						nep(indre)=i
						i1=1
						i2=i1
						call nextcoma(i2,line)
						read(line(i1:i2-1),'(a)')ch
						i1=i2+1
						i2=i1
						call nextcoma(i2,line)
						read(line(i1:i2-1),*)yept(i,indre)
						i1=i2+1
						i2=i1
						call nextcoma(i2,line)
						read(line(i1:i2-1),*)tept(i,indre)
						i1=i2+1
						i2=i1
						call nextcoma(i2,line)
						read(line(i1:i2-1),'(a)')cerfiln(i,indre)
						do icol=1,6
							i1=i2+1
							i2=i1
							call nextcoma(i2,line)
							read(line(i1:i2-1),*)colli(i,indre,icol)
						enddo
						if(ch.eq.'x')then
							kept(i,indre)=1
						elseif(ch.eq.'y')then
							kept(i,indre)=2
						else
							kept(i,indre)=3
						endif
201				tept(i,indre)=tept(i,indre)*dtor
					write(6,*)' descin:too many end-planes in this region !!!'
					stop
202				if(nep(indre).le.0)then
						write(6,*)' descin:no end plane in this region !!!'
						stop
					endif
					yepmin(indre)=yept(1,indre)
					yepmax(indre)=yept(nep(indre),indre)
				endif

c				store the abs. coor. of the submits of the free-box:
				ir=indre
				it=0
				do 5 ix=1,2
				do 5 iy=1,2
				do 5 iz=1,2
					if(cyl(ir))then
						frx=frbox(ix,ir,1)*cos(frbox(iy,ir,2))
						fry=frbox(ix,ir,1)*sin(frbox(iy,ir,2))
						frz=frbox(iz,ir,3)
					else
						frx=frbox(ix,ir,1)
						fry=frbox(iy,ir,2)
						frz=frbox(iz,ir,3)
					endif
					it=it+1
5				call plorot(frx,fry,frz
     @			,boxp(it,ir,1),boxp(it,ir,2),boxp(it,ir,3),xb,yb,zb)

c				region frames:
				nrs(ir)=0

c				region frames not yet implemented in cylindrical regions:
				if(cyl(ir))goto101
				do iq=1,3
					iqp=mod(iq,3)+1
					iqm=mod(iq+1,3)+1
					if(frbox(1,ir,iqp)*frbox(2,ir,iqp).lt.0.
     @			.and.frbox(1,ir,iqm)*frbox(2,ir,iqm).lt.0.
     @			.and.frbox(2,ir,iq).gt.0.)then
						nrs(ir)=nrs(ir)+1

c						starting point:
						rel(iq)=max(0.,frbox(1,ir,iq))
						rel(iqp)=0.
						rel(iqm)=0.
						call plorot(rel(1),rel(2),rel(3)
     @				,regp(1,nrs(ir),ir,1)
     @				,regp(1,nrs(ir),ir,2)
     @				,regp(1,nrs(ir),ir,3)
     @				,xxb,yyb,zzb)

c						ending point:
						rel(iq)=frbox(2,ir,iq)
						rel(iqp)=0.
						rel(iqm)=0.
						call plorot(rel(1),rel(2),rel(3)
     @				,regp(2,nrs(ir),ir,1)
     @				,regp(2,nrs(ir),ir,2)
     @				,regp(2,nrs(ir),ir,3)
     @				,xxb,yyb,zzb)
					endif
				enddo

c				end planes:
				do 100 iepr=1,nep(ir)
					iepz=iepz+1
					if(iepz.gt.maxep)stop'descin:too many end planes!'
					neps(iepz)=0

c					plot of end planes not yet implemented in cylindrical regions:
					if(cyl(ir))goto 100
					if(epsw(ir).eq.'loop')then

c						case 'loop':
						yep=yepmin(ir)+yepstp(ir)*(iepr-1)
						do iq=1,3,2
							iqp=mod(iq,3)+1
							iqm=mod(iq+1,3)+1

c							starting point:
							rel(1)=0.
							rel(3)=0.
							rel(iq)=frbox(1,ir,iq)
							rel(2)=yep
							if(rel(iqp).gt.frbox(1,ir,iqp)
     @					.and.rel(iqp).lt.frbox(2,ir,iqp))then
								neps(iepz)=neps(iepz)+1
								call plorot(rel(1),rel(2),rel(3)
     @						,epp(1,neps(iepz),iepz,1)
     @						,epp(1,neps(iepz),iepz,2)
     @						,epp(1,neps(iepz),iepz,3)
     @						,xxb,yyb,zzb)

c								ending point:
								rel(iq)=frbox(2,ir,iq)
								call plorot(rel(1),rel(2),rel(3)
     @						,epp(2,neps(iepz),iepz,1)
     @						,epp(2,neps(iepz),iepz,2)
     @						,epp(2,neps(iepz),iepz,3)
     @						,xxb,yyb,zzb)
							endif
						enddo

c					case 'list':
					else

c						1/segment along "iqp":
						yep=yept(iepr,ir)
						iq=kept(iepr,ir)
						iqp=mod(iq,3)+1
						iqm=mod(iq+1,3)+1
						tep=tept(iepr,ir)

c						starting point:
						rel(iq)=yep
						rel(iqp)=frbox(1,ir,iqp)
						rel(iqm)=0.
						if(rel(iq).gt.frbox(1,ir,iq)
     @				.and.rel(iq).lt.frbox(2,ir,iq)
     @				.and.rel(iqm).gt.frbox(1,ir,iqm)
     @				.and.rel(iqm).lt.frbox(2,ir,iqm))then

							neps(iepz)=neps(iepz)+1
							call plorot(rel(1),rel(2),rel(3)
     @				,epp(1,neps(iepz),iepz,1)
     @				,epp(1,neps(iepz),iepz,2)
     @				,epp(1,neps(iepz),iepz,3)
     @				,xxb,yyb,zzb)

c							ending point:
							rel(iq)=yep
							rel(iqp)=frbox(2,ir,iqp)
							rel(iqm)=0.
							call plorot(rel(1),rel(2),rel(3)
     @					,epp(2,neps(iepz),iepz,1)
     @					,epp(2,neps(iepz),iepz,2)
     @					,epp(2,neps(iepz),iepz,3)
     @					,xxb,yyb,zzb)
						endif

c						2/segment along "iqm/iq":
						sit=sin(tep)
						if(sit.ne.0.)then
							xl1=(frbox(1,ir,iq)-yep)/sit
							xl2=(frbox(2,ir,iq)-yep)/sit
						else
							xl1=-1.e30
							xl2=1.e30
						endif
						cot=cos(tep)
						if(cot.ne.0.)then
							xl3=frbox(1,ir,iqm)/cot
							xl4=frbox(2,ir,iqm)/cot
						else
							xl1=1.e30
							xl2=1.e30
						endif
						neps(iepz)=neps(iepz)+1

c						starting point:
						xl=max(xl1,xl2)
						xlm=max(xl3,xl4)
						xlmax=min(xl,xlm)
						rel(iq)=yep+xlmax*sit
						rel(iqp)=0.
						rel(iqm)=xlmax*cot
						call plorot(rel(1),rel(2),rel(3)
     @				,epp(1,neps(iepz),iepz,1)
     @				,epp(1,neps(iepz),iepz,2)
     @				,epp(1,neps(iepz),iepz,3)
     @				,xxb,yyb,zzb)

c						ending point:
						xl=min(xl1,xl2)
						xlm=min(xl3,xl4)
						xlmin=max(xl,xlm)
						rel(iq)=yep+xlmin*sit
						rel(iqp)=0.
						rel(iqm)=xlmin*cot
						call plorot(rel(1),rel(2),rel(3)
     @				,epp(2,neps(iepz),iepz,1)
     @				,epp(2,neps(iepz),iepz,2)
     @				,epp(2,neps(iepz),iepz,3)
     @				,xxb,yyb,zzb)
					endif
100			continue
101			continue
				goto200
300			continue
				write(6,*)
				write(6,*)
				return
				end



				subroutine nextcoma(i,line)
c*************************************************************************
c                                                                        *
c                       n e x t c o m a                                  *
c                                                                        *
c*************************************************************************
c				Updates i to the location of the first occurence in 'line' 
c					of a coma at right of the ith. character 
				character line * (*)
				j=i
				do 1 i=j+1,80
					if(line(i:i).eq.',')return
1				continue
				i=i+1
				return
				end



				subroutine plotin(i,x,y,z,bx,by,bz,px,py,pz,flag)
c*************************************************************************
c                                                                        *
c                       p l o t i n                                      *
c                                                                        *
c*************************************************************************
c				add a plot point to the ith trajectory wile keeping
c					the old trajectories
c				if(flag) : create a new plot trajectory and fill
c					its first point
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				logical flag
				if(i+nplvold.gt.maxplv)return
				iplv=i+nplvold
				if(flag)then
					nplv=iplv
					nplp(iplv)=0
				endif
				if(nplp(iplv).lt.maxplp)then
					nplp(iplv)=nplp(iplv)+1
					iplp=nplp(iplv)
					plott(iplp,iplv,1)=x
					plott(iplp,iplv,2)=y
					plott(iplp,iplv,3)=z
					plott(iplp,iplv,4)=bx
					plott(iplp,iplv,5)=by
					plott(iplp,iplv,6)=bz
					plott(iplp,iplv,7)=px
					plott(iplp,iplv,8)=py
					plott(iplp,iplv,9)=pz
				endif
				return
				end



				subroutine box
c*************************************************************************
c                                                                        *
c                       b o x                                            *
c                                                                        *
c*************************************************************************
c				pre-plot the free boxes
c				input:
c					frbox(maxi or mini,region index,x y or z):relat.coord. of the free-box
c						ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c						ih:       "       "           "           "      horizontal "
c						iv:       "       "           "           "      vertical   "
c						kp= 1:observer is at ip>0.
c						-1:   "     "  "  ip<0.
c				output:
c					plobox(submit index,horiz. or vertic.,region index):abs.coord
c						of submit of free-box(special submit list for 3d projection)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				common /comlight/cl(3)
				dimension ls(5,6)
				dimension nh(3),nv(3),np(3)
				data ls/1,3,4,2,1
     @					,5,6,8,7,5
     @					,1,2,6,5,1
     @					,3,7,8,4,3
     @					,1,5,7,3,1
     @					,2,4,8,6,2/
				ul2=0.
				do 5 i=1,3
					nh(i)=0.
					nv(i)=0.
5				ul2=ul2+cl(i)*cl(i)
				ul=sqrt(ul2)
				nh(ih)=1
				nv(iv)=1
				np(1)=nh(2)*nv(3)-nh(3)*nv(2)
				np(2)=nh(3)*nv(1)-nh(1)*nv(3)
				np(3)=nh(1)*nv(2)-nh(2)*nv(1)
				do 6 i=1,3
					if(np(i).ne.0)then
						ip=i
						kp=np(i)
					endif
6				continue
				do 4 ir=1,nure
					indre=ir
					nface(ir)=0 
					do 8 if=1,6

c						put this face in the screen (h*v*p) coord. syst.:
						hl1=boxp(ls(1,if),indre,ih)
						vl1=boxp(ls(1,if),indre,iv)
						pl1=boxp(ls(1,if),indre,ip)*float(kp)
						hl2=boxp(ls(2,if),indre,ih)
						vl2=boxp(ls(2,if),indre,iv)
						pl2=boxp(ls(2,if),indre,ip)*float(kp)
						hl3=boxp(ls(3,if),indre,ih)
						vl3=boxp(ls(3,if),indre,iv)
						pl3=boxp(ls(3,if),indre,ip)*float(kp)

c						is this face seen direct or inverse ?
						dh2=hl2-hl1
						dv2=vl2-vl1
						dp2=pl2-pl1
						dh3=hl3-hl1
						dv3=vl3-vl1
						dp3=pl3-pl1
						sh=dv2*dp3-dp2*dv3
						sv=dp2*dh3-dh2*dp3
						sp=dh2*dv3-dv2*dh3

						if(sp.lt.0.)then

c							this face has to be plotted :
							nface(ir)=nface(ir)+1
							do 7 is=1,5
								plobox(is,nface(ir),1,indre)
     @						=boxp(ls(is,if),indre,ih)
7							plobox(is,nface(ir),2,indre)
     @					=boxp(ls(is,if),indre,iv)
							arg=sh*sh+sv*sv+sp*sp
							if(arg.ne.0.)then
								coef=1./(sqrt(arg)*ul)
							else
								coef=1.
							endif
							clight(nface(ir),indre)=coef*(sh*cl(ih)
     @					+sv*cl(iv)+sp*cl(ip)*float(kp))
							if(nface(ir).ge.4)goto 4
						endif
8					continue
4				continue
				return
				end



				subroutine catch(ah,av,nreg)
c*************************************************************************
c                                                                        *
c                       c a t c h                                        *
c                                                                        *
c*************************************************************************
c				catch the free boxes who plotted a submit closest to cursor
c				input:
c					ah,av(horizontal and vertical cursor lcation 
c						in physical gks coordinates)
c				output:
c				nreg(number of the region whose free box has been caugth)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				logical flag
				d2min=1.e10
				do 1 i=1,nure
					do 2 j=1,5
					do 2 k=1,nface(i)
						d2=(ah-plobox(j,k,1,i))**2+(av-plobox(j,k,2,i))**2
						if(d2.lt.d2min)then
							d2min=d2
							nreg=i
						endif
2					continue
1				continue

c				save the number of clicked boxes and their list
c					supress redundency:
				flag=.false.
				do 3 ic=1,nclic
					if(lclic(ic).eq.nreg)then
						nclic=nclic-1
						flag=.true.
					endif
					if(flag)then
						if(ic.le.nclic)lclic(ic)=lclic(ic+1)
					endif
3				continue
				nclic=nclic+1
				lclic(nclic)=nreg
				return
				end



				subroutine unclip
c*************************************************************************
c                                                                        *
c                       u n c l i p                                      *
c                                                                        *
c*************************************************************************
c				unclip the plot vectors :
				include 'sSnake.inc'
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				do 1 i=1,nplv
					iplpi(i)=1
1				iplpf(i)=nplp(i)
				return
				end



				subroutine clip(sup,iqb,qb)
c*************************************************************************
c                                                                        *
c                           c l i p                                      *
c                                                                        *
c*************************************************************************
c				clip the plot vectors :
				include 'sSnake.inc'
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				do 1 i=1,nplv
					js=1
					do 2 j=1,nplp(i)
						if((qb-plott(j,i,iqb))*sup.gt.0.)goto3
2					js=j
3					iplpi(i)=js
					js=nplp(i)
					do 4 j=nplp(i),iplpi(i),-1
						if((qb-plott(j,i,iqb))*sup.gt.0.)goto 5
4					js=j
5					iplpf(i)=js
1				continue			
				return
				end



				subroutine pretra(*)
c*************************************************************************
c                                                                        *
c                       p r e t r a                                      *
c                                                                        *
c*************************************************************************
c				find the initial window
c				input:
c					nplv<=maxplv:number of trajectories to be plotted
c					nplp(traj.index)<=maxplp:number of submits in the plot-line
c					plott(submit index,traj.index,x y z bx by bz
c						px py pz):absol.coord.of submit,
c					absol. field components and
c					absol. spin components at this point
c						ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c						ih:       "       "           "           "      horizontal "
c						iv:       "       "           "           "      vertical   "
				include 'sSnake.inc'
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /initiale/ coxi,shxi,coyi,shyi
				common /courante/ cox,shx,coy,shy,zmarge
				common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax
				character*10 menu(6)
				dimension XZ(8),YZ(8)
				logical hsvu,hsvz
				dimension work(maxplp)
				parameter (fracti=10.)
				character formh*15,formv*15
				common/ctext/formh,formv
				data menu/'initial','zoom','unzoom','hard copy'
     @						,'paint','exit'/
				data xz/.1,0.,.5,0.,21.,20.5,21.,20.9/
				data yz/.5,0.,.1,0.,29.7,29.6,29.7,29.2/

c				find the plot-window according to ih and iv:
				hsvu=unit(iv).eq.unit(ih)
				xmin= 1.e10
				xmax=-1.e10
				ymin= 1.e10
				ymax=-1.e10
				if(lbox)then

c					then scale also on the free-boxes (11 submits * nure regions)
					do 4 ire=1,nure
					do 4 is=1,5
					do 4 if=1,nface(ire)
						if(plobox(is,if,1,ire).lt.xmin)xmin=plobox(is,if,1,ire)
						if(plobox(is,if,1,ire).gt.xmax)xmax=plobox(is,if,1,ire)
						if(plobox(is,if,2,ire).lt.ymin)ymin=plobox(is,if,2,ire)
4						if(plobox(is,if,2,ire).gt.ymax)ymax=plobox(is,if,2,ire)
				endif

c				then scale on the plot vectors:
				do 1 indve=1,nplv
					np=nplp(indve)
					if(ih.eq.10)then
						call plotl(indve,work)
					elseif(ih.eq.11.or.ih.eq.12.or.ih.eq.15)then
						call plotm(indve,ih,work)
					elseif(ih.eq.13)then
						call plota(indve,work)
					elseif(ih.eq.14)then
						call plotra(indve,work)
					else
						call plotg(indve,ih,work)
					endif
					do 2 n=1,np
						if(work(n).lt.xmin)xmin=work(n)
2						if(work(n).gt.xmax)xmax=work(n)
					if(iv.eq.10)then
						call plotl(indve,work)
					elseif(iv.eq.11.or.iv.eq.12.or.iv.eq.15)then
						call plotm(indve,iv,work)
					elseif(iv.eq.13)then
						call plota(indve,work)
					elseif(iv.eq.14)then
						call plotra(indve,work)
					else
						call plotg(indve,iv,work)
					endif
					do 3 n=1,np
						if(work(n).lt.ymin)ymin=work(n)
3					if(work(n).gt.ymax)ymax=work(n)
1				continue

c				add 1% margins:
				hm=1.e-2*(xmax-xmin)
				xmin=xmin-hm
				xmax=xmax+hm
				vm=1.e-2*(ymax-ymin)
				ymin=ymin-vm
				ymax=ymax+vm
				hsvz=hsvu

c				open gks and graphic screen
				xmini=xmin
				xmaxi=xmax
				ymini=ymin
				ymaxi=ymax
				iwr=1
				zmarge=2.
				call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz
     @				,cox,shx,coy,shy)
				coxi=cox
				shxi=shx
				coyi=coy
				shyi=shy
29			continue
c				CALL IGSSE(6,1)
				call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz
     @		,cox,shx,coy,shy)

c				higz: the size of the physical window on the wk is fixed by the ith. line
c					of the file 'higz_windows.dat' (if it exists ) , where 'i' is the second
c					argument of igsse. in that line data are: h(def=0),v(def=0),dh(def=600) and
c					dv(def=600). Unit is pixel. (h,v) is the position of the upper left corner,
c					(dh,dv) is the size of the window. for sun h<1052, v<913. Inside this 
c					window, higz plots on the biggest viewport it can draw, keeping its dimensions
c					h' and v' in the same ratio as d2 and d3 .(Keeping also the aspect ratio)
				if(xmin.eq.xmax.or.ymin.eq.ymax)then
c				call igend
c				CALL IGSA (1)
				write(6,*)' pretra:can''t scale'
				return1
				endif
30			continue

c				CALL ICLRWK(0,1)
c				CALL IGRNG(21.,29.7)

c				plot global coord
c				call iselnt(2)
c				call isclip(1)
c				call isln(3)

c				if(.not.bare)then
c					do 32 i=1,nq
c						call ipl(2,h1(1,i),v1(1,i))
c32				continue    
c				endif

c				call isln(1)
c				call iselnt(0)
				if(.not.bare)call grad
c				call iselnt(2)

				alp=(xmax1-xmin1)/(xnmax-xnmin)
				xp=xmin1+(30.-xnmin)*alp
				xm=xmin1-(30.+xnmin)*alp
				alp=(ymax1-ymin1)/(ynmax-ynmin)
				yp=ymin1+(30.-ynmin)*alp
				ym=ymin1-(30.+ynmin)*alp
				call traj(xp,xm,yp,ym)

c				plot in cm:
c				call iselnt(1)
c				call ipl(3,xz,yz)
c				call ipl(3,xz(6),yz(6))
				if(iwr.eq.2)then
					iwr=1
c					call igend
					close(10)
					goto 29
				endif

c				plot menu in cm:
c				call iselnt(1)
c				call ischh(.35)
c				call istxci(1)
				dx=20.9/6.
				dz=dx/10.
				y1=28.5+.1
				y2=29.7-dz
				x1=.1
				do i=1,6
					x2=x1+dx-dz
c					call igpave(x1,x2,y1,y2,dz,0.,0.,'TR')
c					call itx(x1+dz,y1+dz/2.,menu(i))
					x1=x1+dx
				enddo
c				call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
31			continue   
33			xw=xndc*29.7
				yw=yndc*29.7
				if(yw.ge.y1.and.yw.le.y2)then	
					x1=.1
					do i=1,6
						x2=x1+dx-dz
						if(xw.gt.x1.and.xw.le.x2)then
							goto(101,102,103,104,105,106),i
						endif
						x1=x1+dx
					enddo
					goto 31
				endif
				goto 31

c				initial:
101			write(6,*)menu(1)
				hsvz=hsvu
				xmin=xmini
				xmax=xmaxi
				ymin=ymini
				ymax=ymaxi
				cox=coxi
				shx=shxi
				coy=coyi
				shy=shyi
				call map1(xmini,ymini,xmaxi,ymaxi,10,21.,28.5
     @						,zmarge,hsvz,coxi,shxi,coyi,shyi)
				goto 30

c				zoom:
102			write(6,*)menu(2)
c				call igloc(1,nt,ibn,xndc1,yndc1,xwc,ywc)
				x1=xndc1*29.7
				y1=yndc1*29.7
				call cross(x1,y1)
c				call igloc(1,nt,ibn,xndc2,yndc2,xwc,ywc)
				x2=xndc2*29.7
				y2=yndc2*29.7
				if(abs(x2-x1).lt..1.and.abs(y2-y1).lt..1)then
					hsvz=.false.
c					call igloc(1,nt,ibn,xndc2,yndc2,xwc,ywc)
				else
					hsvz=hsvu
				endif
				xmin=(min(xndc1,xndc2)-zmarge/29.7)/cox+shx
				ymin=(min(yndc1,yndc2)-zmarge/29.7)/coy+shy
				xmax=(max(xndc1,xndc2)-zmarge/29.7)/cox+shx
				ymax=(max(yndc1,yndc2)-zmarge/29.7)/coy+shy
				call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz
     @				,cox,shx,coy,shy)
				goto 30

c				unzoom:
103			write(6,*)menu(3)
c				call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
				x1=(xndc-zmarge/29.7)/cox+shx
				y1=(yndc-zmarge/29.7)/coy+shy
c				call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
				x2=(xndc-zmarge/29.7)/cox+shx
				y2=(yndc-zmarge/29.7)/coy+shy
				alp=(xmax-xmin)/(max(x1,x2)-min(x1,x2))
				bet=xmin-alp*min(x1,x2)
				xmin=alp*xmin+bet
				xmax=alp*xmax+bet
				alp=(ymax-ymin)/(max(y1,y2)-min(y1,y2))
				bet=ymin-alp*min(y1,y2)
				ymin=alp*ymin+bet
				ymax=alp*ymax+bet
				call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz
     @			,cox,shx,coy,shy)
				goto 30

c				hard copy:
104			write(6,*)menu(4)
				iwr=2
c				call igend
				open(10,file='higz.ps',status='unknown')
c				call iopks(6)

c				PostScript metafile  A4 Portrait (-111):
c				call iopwk(1,10,-111)
c				call iacwk(1)
				call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz
     @			,cox,shx,coy,shy)
				go to 30

c				paint:
105			write(6,*)menu(5)
c				call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
				if(yndc*29.7.ge.y1)goto 33
				ah=(xndc-zmarge/29.7)/cox+shx
				av=(yndc-zmarge/29.7)/coy+shy
				call catch(ah,av,nreg)
c				call iselnt(2)
				call paint(nreg)
				goto 105

c				exit:
106			write(6,*)menu(6)
c				call igend
				return
				END



				subroutine map1(xmin,ymin,xmax,ymax,n,xsize1,ysize1,zmarge,hsvz
     @						,cox,shx,coy,shy)
c*************************************************************************
c                                                                        *
c                       m a p 1                                          *
c                                                                        *
c*************************************************************************
				include 'sSnake.inc'
				logical hsvz
				parameter (fracti=10.)
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				character formh*15,formv*15
				real a
				common/ctext/formh,formv
				common/text/hldep,bh,vldep,bv,iexph,iexpv,ith,itv
				common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax
				external scale
				xnmin=zmarge/29.7
				xnmax=(xsize1-zmarge)/29.7
				ynmin=zmarge/29.7
				ynmax=(ysize1-zmarge)/29.7
				hsvzn=(xnmax-xnmin)/(ynmax-ynmin)
				if(hsvz)then
					dx=xmax-xmin
					dy=ymax-ymin
					if(dx/dy.gt.hsvzn)then
						ymoy=(ymax+ymin)*.5
						dys2=.5*dx/hsvzn
						xmin1=xmin
						xmax1=xmax
						ymin1=ymoy-dys2
						ymax1=ymoy+dys2
					else
						xmoy=(xmax+xmin)*.5
						dxs2=.5*dy*hsvzn
						xmin1=xmoy-dxs2
						xmax1=xmoy+dxs2
						ymin1=ymin
						ymax1=ymax
					endif
				else
					xmin1=xmin
					xmax1=xmax
					ymin1=ymin
					ymax1=ymax
				endif
c				call iswn(2,xmin1,xmax1,ymin1,ymax1)
c				call isvp(2,xnmin,xnmax,ynmin,ynmax)
				cox=(xnmax-xnmin)/(xmax1-xmin1)
				shx=xmin1
				coy=(ynmax-ynmin)/(ymax1-ymin1)
				shy=ymin1

c				compute the scale and store the tic marks:
				dtic=(xmax1-xmin1)/fracti
				nq=0
				call scale(xmin1,xmax1,dtic,a,bh,formh,ith,iexph)
				hldep=a+bh
				do 3 hl=hldep,xmax1,bh
					if(nq.ge.maxq)goto5
					nq=nq+1
					h1(1,nq)=hl
					v1(1,nq)=ymin1
					h1(2,nq)=hl
3				v1(2,nq)=ymax1
				dtic=(ymax1-ymin1)*hsvzn/fracti
				call scale(ymin1,ymax1,dtic,a,bv,formv,itv,iexpv)
				vldep=a+bv
				do 4 vl=vldep,ymax1,bv
					if(nq.ge.maxq)goto 5
					nq=nq+1
					h1(1,nq)=xmin1
					v1(1,nq)=vl
					h1(2,nq)=xmax1
4				v1(2,nq)=vl
5				return
				end



				subroutine paint(nreg)
c*************************************************************************
c                                                                        *
c                       p a i n t                                        *
c                                                                        *
c*************************************************************************
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension icoll(3)

c				call iselnt(2)
c				set: style=solid
c				call isln(1)
c				call isfais(1)

c				paint the 3 faces by filling them with background or foreground
c					color depending on the lightning:
				do 6 if=1,nface(nreg)
					if(clight(if,nreg).lt.0.)then
						icols=0
						icoll(if)=1
					else
						icols=1
						icoll(if)=0
					endif
c					call isfaci(icols)
6				continue
c				call ifa(5,plobox(1,if,1,nreg),plobox(1,if,2,nreg))

c				re-draw the ridges of the 3 faces (erased by gfa) :
c				do 7 if=1,nface(nreg)
c				call isplci(icoll(if))
c7				continue

c				call ipl(5,plobox(1,if,1,nreg),plobox(1,if,2,nreg))
c				call isplci(1)
c				call igterm
				return
				end



				subroutine cross(h,v)
c*************************************************************************
c                                                                        *
c                       c r o s s                                        *
c                                                                        *
c*************************************************************************
c				plot a cross at the position h*v (normalized coordinates)
				dimension hch(2),hcv(2),vch(2),vcv(2)
				hch(1)=h
				hch(2)=h
				hcv(1)=0.
				hcv(2)=21.
				vch(1)=0.
				vch(2)=28.5
				vcv(1)=v
				vcv(2)=v
c				call isln(2)
c				call ipl(2,hch(1),vch(1))
c				call ipl(2,hcv(1),vcv(1))
				return
				end



				subroutine plotl(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t l                                        *
c                                                                        *
c*************************************************************************
c				for the trajectory of index=indve :
c					integrate the track-length from the data of plott
c					and put it in work(maxplp)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu 
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension work(maxplp)
				work(1)=0.
				do 1 i=2,nplp(indve)
					dl2=0.
					do 2 j=1,3
2					dl2=dl2+(plott(i,indve,j)-plott(i-1,indve,j))**2
1				work(i)=work(i-1)+sqrt(dl2)
				return
				end



				subroutine plota(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t a                                        *
c                                                                        *
c*************************************************************************
c				for the trajectory of index=indve :
c					compute the angular component (t) of cylindrical coordinates (r*t*z)
c					and put it in work(maxplp)
c				axis of the cylindrical coord.: z
c				origin of the angles on axis x:  (t=arctan(y/x))
c				warning!: t is given in deg. , the starting point is taken
c					between -180. and 180. but the other points of the trajectory
c					can have their t outside this range .
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension work(maxplp)
				scal(xn,yn,xa,ya)=xn*xa+yn*ya
				vect(xn,yn,xa,ya)=xn*ya-yn*xa
				x1=plott(1,indve,1)
				y1=plott(1,indve,2)

c				define starting angle in  ( -pi , pi ) :
				if(x1.eq.0..and.y1.eq.0.)then
					work(1)=0.
				else
					work(1)=atan2(y1,x1)*rtod
				endif

c				next code allows the angle to go outside ( -pi , pi ) :
				do 1 i=2,nplp(indve)
					x=plott(i,indve,1)
					y=plott(i,indve,2)
					if((x1.eq.0..and.y1.eq.0.).or.(x.eq.0..and.y.eq.0.))then
						da=0.
					else
						da=atan2(vect(x1,y1,x,y),scal(x1,y1,x,y))
					endif
					work(i)=work(i-1)+da*rtod
					x1=x
1				y1=y
				return
				end



				subroutine plotra(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t r a                                      *
c                                                                        *
c*************************************************************************
c				for the trajectory of index=indve :
c					compute the angular component (t) by calling plota,
c					the radius (r) and put their product in work(maxplp)
c				axis of the cylindrical coord.: z
c				origin of the angles on axis x:    (t=arctan(y/x))
c				warning!: t is given in deg. , the starting point is taken
c					between -180. and 180. but the other points of the trajectory
c					can have their t outside this range .
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension work(maxplp)
				call plota(indve,work)
				do 1 i=1,nplp(indve)
					b2=0.
					do 2 j=1,2
2					b2=b2+plott(i,indve,j)**2
1				work(i)=work(i)*dtor*sqrt(b2)
				return
				end



				subroutine plotg(indve,iq,work)
c*************************************************************************
c                                                                        *
c                       p l o t g                                        *
c                                                                        *
c*************************************************************************
c				for the trajectory of index=indve :
c					copy plott(i,indve,iq) in work(i)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension work(maxplp)
				do 1 i=1,nplp(indve)
1				work(i)=plott(i,indve,iq)
				return
				end
				subroutine plotm(indve,iq,work)
c*************************************************************************
c                                                                        *
c                       p l o t m                                        *
c                                                                        *
c*************************************************************************
c				for the trajectory of index=indve :
c					compute the module of the field or the radius from the data of plott
c					and put it in work(maxplp)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				dimension work(maxplp)
				if(iq.eq.11)then
					j1=4
					j2=6
				else if(iq.eq.12)then
					j1=1
					j2=2
				else if(iq.eq.15)then
					j1=7
					j2=9
				else
					stop'plotm: invalid "iq" argument'
				endif
				do 1 i=1,nplp(indve)
					b2=0.
					do2 j=j1,j2
2         b2=b2+plott(i,indve,j)**2
1				work(i)=sqrt(b2)
				return
				end



				subroutine traj(xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       t r a j                                          *
c                                                                        *
c*************************************************************************
c				plot trajectories,tic marks ( and free-boxes if necessary)
c				input:
c					nure:number of region(s)  defined
c					nplv<=maxplv:number of trajectories to be plotted
c					nplp(traj.index)<=maxplp:number of submits in the plot-line
c					plott(submit index,traj.index,x y z bx by bz px py pz):absol.coord.of submit,
c					absol. field components and
c					absol. spin components at this point
c					nq<=maxq:number of tic marks to be plotted
c					h1(submit index,tic mark index):absol.horiz.coord of the submit
c						of the tic mark
c					v1(submit index,tic mark index):absol.vertic.coord of the submit
c						of the tic mark
c					ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c					ih:       "       "           "           "      horizontal "
c					iv:       "       "           "           "      vertical   "
c					plobox(submit index,horiz. or vertic.,region index):abs.coord
c						of submit of free-box(special submit list for 3d projection)
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				dimension workh(maxplp),workv(maxplp)

c				trace vectors:
c				call isln(1)
				do 1 i=1,nplv
					if(ih.eq.10)then
						call plotl(i,workh)
					elseif(ih.eq.11.or.ih.eq.12.or.ih.eq.15)then
						call plotm(i,ih,workh)
					elseif(ih.eq.13)then
						call plota(i,workh)
					elseif(ih.eq.14)then
						call plotra(i,workh)
					else
						call plotg(i,ih,workh)
					endif
					if(iv.eq.10)then
						call plotl(i,workv)
					elseif(iv.eq.11.or.iv.eq.12.or.iv.eq.15)then
						call plotm(i,iv,workv)
					elseif(iv.eq.13)then
						call plota(i,workv)
					elseif(iv.eq.14)then
						call plotra(i,workv)
					else
						call plotg(i,iv,workv)
					endif
					iplp=iplpf(i)-iplpi(i)+1

c					if(iplp.gt.1)then
c						call ipl(iplp,workh,workv)

c					HIGZ limit 1000 points:
					iplpc=iplp
					iplpini=iplpi(i)
200				if(iplpc.gt.1000)then
c						call ipl1(1000,workh(iplpini),workv(iplpini)
c    @				,xp,xm,yp,ym)
						iplpini=iplpini+999
						iplpc=iplpc-999
						goto 200
					endif

c					call ipl1(iplpc,workh(iplpini),workv(iplpini),xp,xm,yp,ym)
c					call igterm
c					endif
1				continue

c				trace user plot:
c				call isln(1)
c				if(lbox.and.luser)then
c					do i=1,nulu
c						call ipl1(nuseu(i),plotus(1,i,ih),plotus(1,i,iv)
c     @				,xp,xm,yp,ym)
c					end do
c				endif

c				trace free-boxes:
c				first: not clickd boxes
c				call isln(1)
				if(lbox.and..not.bare)then
c					do 3 i=1,nure
c						do 4 ic=1,nclic
c4						if(i.eq.lclic(ic))goto 3
c						do 6 if=1,nface(i)	
c							call ipl1(5,plobox(1,if,1,i),plobox(1,if,2,i)
c     @					,xp,xm,yp,ym)
c6						continue
c3					continue

c					second: clicked boxes,in chronological order
					do 5 i=1,nclic
						nreg=lclic(i)
						call paint(nreg)
5					continue

c					plot region frames:
c					call isplci(2)
c					do ir=1,nure
c						do irs=1,nrs(ir)
c							call ipl1(2,regp(1,irs,ir,ih),regp(1,irs,ir,iv)
c     @					,xp,xm,yp,ym)
c						enddo
c					enddo
c					call isplci(1)

c					plot end-planes:
c					call isplci(4)
c					do iep=1,nuep
c						do ieps=1,neps(iep)
c							call ipl1(2,epp(1,ieps,iep,ih),epp(1,ieps,iep,iv)
c     @					,xp,xm,yp,ym)
c						enddo
c					enddo

c					call isplci(1)
				endif
				return
				end



				subroutine grad
c*************************************************************************
c                                                                        *
c                       g r a d                                          *
c                                                                        *
c*************************************************************************
c				plot numerical values of cordinates corresponding to the tic marks
c				for simplicity the plot is done in gks's normalized coordinates
c				input:
c					/text/:data to plot strings
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				character formh*15,formv*15,string*20
				common/ctext/formh,formv
				common/text/hldep,bh,vldep,bv,iexph,iexpv,ith,itv
				common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax

c				call ischh(1./120.)
c				call ischup(-1.,0.)
c				call istxci(1)

				coef=10.**(-iexph)
				alp=(xnmax-xnmin)/(xmax1-xmin1)
				do 13 hl=hldep,xmax1,bh
					hlco=hl*coef
					write(string,formh)hlco
					hnp=hn
					hn=xnmin+hl*alp-xmin1*alp
c					call itx(hn,ynmin,string)
13			continue
c				call ischup(0.,1.)
				coef=10.**(-iexpv)
				alp=(ynmax-ynmin)/(ymax1-ymin1)
				do 16 vl=vldep,ymax1,bv
					vlco=vl*coef
					write(string,formv)vlco
					vnp=vn
					vn=ynmin+vl*alp-ymin1*alp
c					call itx(xnmin,vn,string)
16			continue
c				call ischh(1./60.)
				vn=(vnp*5.+vn*3.)/8.
c				call itx(xnmin,vn,name(iv)//'-'//unit(iv))
c				call ischup(-1.,0.)
				hn=(hnp*3.+hn*5.)/8.
c				call itx(hn,ynmin,name(ih)//'-'//unit(ih))
c				call ischup(0.,1.)
c				call ischh(1./3.)
				return
				end



				subroutine ipl1(n,x,y,xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       i p l 1                                          *
c                                                                        *
c*************************************************************************
c				La routine ipl de Higz genere un trace incorrect lorsque la polyline
c					a l'un de ses sommets situe a grande distance de la fenetre effectivement
c					plottee (environ 38 en coordonnee normalisee pour la version 93c de cernlib).
c					ipl1(n,x,y,xp,xm,yp,ym) se substitue a ipl(n,x,y) en assurant un
c					ecretage sur le rectangle (xp,xm,yp,ym) (coord. world) suppose 
c					correspondre a x=+-30, y=+-30 (coord. normalisees).
c				La polyline initiale est decoupee en sous-polylines dont seuls le point
c					origine et le point extremite sortent eventuellement du cadre. Ces 2 points
c					sont alors rabatus sur le cadre par interpolation lineaire, vers le 
c					point suivant pour le premier, vers le point precedent pour le dernier.
c					Les parties de polylines entierement a l'exterieur du cadre sont converties
c					en points de plot (polyline a 2 sommets coincidents) et ne genent par
c					le plot.
				parameter (maxpoly=1000)
				dimension x(*),y(*),xt(maxpoly),yt(maxpoly)
				logical previous_in
				if(n.gt.maxpoly)then
				write(6,*)'ipl1: polyline size limited to ',maxpoly,'!'
				stop
				endif
				do ii=1,n
					if(x(ii).lt.xp.and.x(ii).gt.xm
     @			.and.y(ii).lt.yp.and.y(ii).gt.ym)then
						if(ii.eq.1)then
							if=1
							xt(if)=x(ii)
							yt(if)=y(ii)
						elseif(previous_in)then
							if=if+1
							xt(if)=x(ii)
							yt(if)=y(ii)
c							if(ii.eq.n)then
c								call ipl(if,xt,yt)
c								write(10,*)if,xt,yt
c							endif
						else
							if=1
							call rabat(ii-1,ii,if,x,y,xt,yt
     @					,xp,xm,yp,ym)
							if=if+1
							xt(if)=x(ii)
							yt(if)=y(ii)
c							if(ii.eq.n)then
c								call ipl(if,xt,yt)
c							endif
						endif
						previous_in=.true.
					else
						if(ii.ne.1.and.previous_in)then
							if=if+1
							call rabat(ii,ii-1,if,x,y,xt,yt
     @					,xp,xm,yp,ym)
c							call ipl(if,xt,yt)
						elseif(ii.ne.1)then
							den=x(ii)-x(ii-1)
							xl1=2.
							if(den.ne.0.)xl1=(xp-x(ii-1))/den
							y1=xl1*y(ii)+(1.-xl1)*y(ii-1)
							xl2=2.
							if(den.ne.0.)xl2=(xm-x(ii-1))/den
							y2=xl2*y(ii)+(1.-xl2)*y(ii-1)
							den=y(ii)-y(ii-1)
							xl3=2.
							if(den.ne.0.)xl3=(yp-y(ii-1))/den
							x1=xl3*x(ii)+(1.-xl3)*x(ii-1)
							if((xl1.gt.0..and.xl1.lt.1.
     @					.and.y1.gt.ym.and.y1.lt.yp)
     @					.or.(xl2.gt.0..and.xl2.lt.1.
     @					.and.y2.gt.ym.and.y2.lt.yp)
     @					.or.(xl3.gt.0..and.xl3.lt.1.
     @					.and.x1.gt.xm.and.x1.lt.xp))then
								if=1
								call rabat(ii-1,ii,if,x,y,xt,yt
     @						,xp,xm,yp,ym)
								if=if+1
								call rabat(ii,ii-1,if,x,y,xt,yt
     @						,xp,xm,yp,ym)
c								call ipl(if,xt,yt)
							endif
						endif
						previous_in=.false.
					endif
				enddo
				return
				end



				subroutine ipl1old(n,x,y,xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       i p l 1                                          *
c                                                                        *
c*************************************************************************
c				La routine ipl de Higz genere un trace incorrect lorsque la polyline
c					a l'un de ses sommets situe a grande distance de la fenetre effectivement
c					plottee (environ 38 en coordonnee normalisee pour la version 93c de cernlib).
c					ipl1(n,x,y,xp,xm,yp,ym) se substitue a ipl(n,x,y) en assurant un
c					ecretage sur le rectangle (xp,xm,yp,ym) (coord. world) suppose 
c					correspondre a x=+-30, y=+-30 (coord. normalisees).
c				La polyline initiale est decoupee en sous-polylines dont seuls le point
c					origine et le point extremite sortent eventuellement du cadre. Ces 2 points
c					sont alors rabatus sur le cadre par interpolation lineaire, vers le 
c					point suivant pour le premier, vers le point precedent pour le dernier.
c				Les parties de polylines entierement a l'exterieur du cadre sont converties
c					en points de plot (polyline a 2 sommets coincidents) et ne genent par
c					le plot.
				parameter (maxpoly=1000)
				dimension x(*),y(*),xt(maxpoly),yt(maxpoly)
				if(n.gt.maxpoly)then
					write(6,*)'ipl1: polyline size limited to ',maxpoly,'!'
					stop
				endif
				idep=0
2				if(x(idep+1).gt.xp.or.x(idep+1).lt.xm
     @		.or.y(idep+1).gt.yp.or.y(idep+1).lt.ym)then
					call rabat(idep+1,idep+2,1,x,y,xt,yt,xp,xm,yp,ym)
				else
					xt(1)=x(idep+1)
					yt(1)=y(idep+1)
				endif
				do 1 i=idep+2,n
					ic=i
					if(x(i).gt.xp.or.x(i).lt.xm
     @			.or.y(i).gt.yp.or.y(i).lt.ym)then
						call rabat(i,i-1,i-idep
     @				,x,y,xt,yt,xp,xm,yp,ym)
c						call ipl(i-idep,xt,yt)
						if(i.eq.n)return
						idep=i-1
						goto 2
					else
						xt(i-idep)=x(i)
						yt(i-idep)=y(i)
					endif
1				continue
c				call ipl(ic-idep,xt,yt)
				return
				end



				subroutine rabat(i,ip,j,x,y,xt,yt,xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       r a b a t                                        *
c                                                                        *
c*************************************************************************
c 
c				copie dans (xt,yt) a l'adresse j le point (x,y), adresse i, ecrete par
c					le rectangle (xp,xm,yp,ym). L'interpolation est faite de (x,y)(i) vers
c					(x,y)(ip) avec ip=i+-1. Le point (x,y)(ip) est suppose situe de sorte 
c					que le segment i-->ip intercepte le rectangle (xp,xm,yp,ym)
				dimension x(*),y(*),xt(*),yt(*)
				xc=x(i)
				yc=y(i)
				if(xc.gt.xp)then
					yc=(y(ip)*(xc-xp)+yc*(xp-x(ip)))/(xc-x(ip))
					xc=xp
				elseif(xc.lt.xm)then
					yc=(y(ip)*(xc-xm)+yc*(xm-x(ip)))/(xc-x(ip))
					xc=xm
				endif
				if(yc.gt.yp)then
					xc=(x(ip)*(yc-yp)+xc*(yp-y(ip)))/(yc-y(ip))
					yc=yp
				elseif(yc.lt.ym)then
					xc=(x(ip)*(yc-ym)+xc*(ym-y(ip)))/(yc-y(ip))
					yc=ym
				endif
				xt(j)=xc
				yt(j)=yc
				return
				end



				subroutine scale(ba,bb,fn,a,fnsr,form,it,iexp)
c*************************************************************************
c                                                                        *
c                       s c a l e                                        *
c                                                                        *
c*************************************************************************
c				find the most economic way to mark from ba to bb
c					with an increment aproximately=fn
c				input:
c					ba,bb and fn
c				output:
c					a(near ba) and fnsr(near fn)
c					form:most economical write format
c					it:number of chatacter(s) in form
c					iexp:value of the decimal exponent (if any)
				character form*15,formm*8,forme*7,formf*17
				b1=ba
				b2=bb
				if(b1.eq.b2)then

c					in this case:a and b arbitrary:
					a=b1
					b=1.
					return
				endif
				call c125(abs(fn),b,imin)
				fnsr=sign(b,b2-b1)
				ap=abs(b1)
				ap=ap-amod(ap,b)
				a=sign(ap,b1)

c				find first and last usefull digit:
				if(a.ne.0.)then
					i1max=ifix(alog10(abs(a))+.0001)
					if(abs(a).lt.1.)i1max=i1max-1
				else
					i1max=-1
				endif
				if(b2.ne.0.)then
					i2max=ifix(alog10(abs(b2))+.0001)
					if(abs(b2).lt.1.)i2max=i2max-1
				else
					i2max=-1
				endif
				imax=max(i1max,i2max)
				imoy=(imin+imax)/2+1
				iexp=3*(imoy/3)
				if(imin.le.2.and.imax.ge.-3)iexp=0
				imin=imin-iexp
				imax=imax-iexp
				ifmin=min(imin,0)
				ifmax=max(imax,-1)
				it=ifmax-ifmin+2
				if(b1.lt.0..or.b2.lt.0.)it=it+1
				if=-ifmin
				write(formm,100)it,if
100			format('(f',i2,'.',i2)
				forme=')'
				if(iexp.ne.0)then
					nexp=ifix(alog10(abs(float(iexp))+.0001))+1
					if(iexp.lt.0)nexp=nexp+1
					write(formf,200)nexp
200				format('('',''''e'',i',i1,','''''')'')')
					write(forme,formf)iexp
					it=it+1+nexp
				endif
				form=formm//forme
				return
				end



				subroutine c125(x,xp,imin)
c*************************************************************************
c                                                                        *
c                       c 1 2 5                                          *
c                                                                        *
c*************************************************************************
c				for x>0. find xp as close as possible from x
c					such that: xp=(1. , 2. , 2.5  or 5.)*10.**n
c				input: x
c				output:
c					xp
c					imin:distance in digit between the last stginficant digit
c						and the decimal point( >0 at right)
				dimension fl(5),xl(5),id(5)
				data fl/1.,2.,2.5,5.,10./
				data id/0 ,0 ,-1 ,0 , 1 /
				do 18 i=1,5
					xl(i)=alog10(fl(i))
18			continue
				xlc=alog10(x)
				ixl=int(xlc+.0001)
				if(x.lt.1.)ixl=ixl-1
				dixl=10.**ixl
				xl0=xlc-ixl
				do 21 i=2,5
					if(xl0.lt.xl(i))then
						xp=fl(i)*dixl
						goto1
					endif
21			continue
1				xm=fl(i-1)*dixl
				imin=id(i)+ixl
				if(xp/x.gt.x/xm)then
					xp=xm
					imin=id(i-1)+ixl
				endif
				return
				end



				subroutine itrans
c*************************************************************************
c                                                                        *
c                       i t r a n s                                      *
c                                                                        *
c*************************************************************************
				include 'sSnake.inc'
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				call plorot(0.,0.,0.,x0,y0,z0,xb,yb,zb)
				call plorot(1.,0.,0.,xb,yb,zb,xx,xy,xz)
				call plorot(0.,1.,0.,xb,yb,zb,yx,yy,yz)
				call plorot(0.,0.,1.,xb,yb,zb,zx,zy,zz)
				irsave=indre
				do 3 i=1,indre-1
					if(rnamer(indre).eq.rname(i))then
						indrer(indre)=i
						goto4
					endif
3				continue
				write(6,*)' can''t find this region :',rnamer(indre)
				stop
4				indre=indrer(indre)
				call plorot(x0,y0,z0,x0,y0,z0,xb,yb,zb)
				call plorot(xx,xy,xz,xb,yb,zb,xx,xy,xz)
				call plorot(yx,yy,yz,xb,yb,zb,yx,yy,yz)
				call plorot(zx,zy,zz,xb,yb,zb,zx,zy,zz)
				indre=irsave
				xyz0(indre,1)=x0
				xyz0(indre,2)=y0
				xyz0(indre,3)=z0
				xang(indre)=asin(yz)
				if((xz.eq.0.and.zz.eq.0.).or.(yx.eq.0..and.yy.eq.0.))then

c					polar ambiguity:
					yang(indre)=0.
					zang(indre)=atan2(xy,xx)
				else
					yang(indre)=atan2(-xz,zz)
					zang(indre)=atan2(-yx,yy)
				endif
				zad=zang(indre)*rtod
				xad=xang(indre)*rtod
				yad=yang(indre)*rtod
				return
				end



				subroutine rot
c*************************************************************************
c                                                                        *
c                       r o t                                            *
c                                                                        *
c*************************************************************************
c				at each region convert absolute coord. in relative one
c					to fill the current relative coord. buffer
c				input:
c					nuve:number of input vector(s) defined
c					vea(iv,iep,ic):description of the particules at each end-plane
c					iv<=nuve<=maxve:vector index
c					iep<=nuep<=maxep:end-plane index
c					ic= 1:x coord. in absolute frame
c							2:y coord
c							3:z coord
c							4:x component of the normalized mommentum in absol.frame
c							5:y momentum
c							6:z momentum
c							7:module of the mommentum
c							8:trace length measured from the "spring point"
c							9:life-flag(-1.:buried,0.:just dead,1.:alive)
c							10 to 15:the same as 1 to 6 but in relative frame coordinates
c							16:x component of the spin
c							17:y spin
c							18:z spin
c							19 to 21: the same as 16 to 18 but in relative frame coord.
c					indre=current region   index
c					indep=current end-plane index
c					xyz0(region index,x y or z):absolute coord. of the rel.origine
c					zang(region index),xang(region index) and yang(region index):
c						rotation angles to define the relat.frame in the absol.frame
c					adata(data index < maxdat,region index):additional data to describe a region
c				output:
c					ver(iv,ic):description of the particules in the current region
c					iv<=nuve<=maxve:vector index
c					ic= 1:x coord.relative to the current region(in relative frame)
c							2:y            "
c							3:z            "
c							4:x component of the normalized mommentum(in relat.frame)
c							5:y            "
c							6:z            "
c							7:module of the mommentum
c							8:trace length measured from the "spring point"
c							9:life-flag(-1.:buried,0.:just dead,1.:alive)
c							10:x component of the spin (in relat.frame)
c							11:y          "
c							12:z          "
c				output for subroutine rot and plorot:
c					x0,y0,z0,stx,ctx,stz,ctz,sty,cty=rotation matrix elements
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /crot/x0,y0,z0,stx,ctx,stz,ctz,sty,cty
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				common/regin/veri(maxve,maxre,12)
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				logical cyl,lfil
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @					,cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				x0=xyz0(indre,1)
				y0=xyz0(indre,2)
				z0=xyz0(indre,3)
				stz=sin(zang(indre))
				ctz=cos(zang(indre))
				stx=sin(xang(indre))
				ctx=cos(xang(indre))
				sty=sin(yang(indre))
				cty=cos(yang(indre))
				do 1 i=1,nuve

c				if this vector is dead or buried, no rotation:
					if(vea(i,indep,9).le.0.)then
						ver(i,9)=vea(i,indep,9)
						do 3 iq=1,12
3						veri(i,indre,iq)=0.
						veri(i,indre,9)=vea(i,indep,9)
						goto1
					endif

c					no field in the new region:
c					plot these points to zero the field:
c					next not usefull for the first region: this plot has been done 
c						by inject - 92 
					if(method(indre).le.2)then
						if(indre.ne.1)call plotin(i,vea(i,indep,1)
     @												,vea(i,indep,2)
     @												,vea(i,indep,3)
     @												,0.,0.,0.,vea(i,indep,16)
     @												,vea(i,indep,17)
     @												,vea(i,indep,18),.false.)
					endif

c					first:translation (the new origine is the point "0"):
					x=vea(i,indep,1)-x0
					y=vea(i,indep,2)-y0
					z=vea(i,indep,3)-z0
					cx=vea(i,indep,4)
					cy=vea(i,indep,5)
					cz=vea(i,indep,6)
					px=vea(i,indep,16)
					py=vea(i,indep,17)
					pz=vea(i,indep,18)

c					second:rotations if necessary:
c						in this case:rotation/z axis
					if(zang(indre).ne.0.)then
						xt=x*ctz+y*stz
						y=y*ctz-x*stz
						x=xt
						cxt=cx*ctz+cy*stz
						cy=cy*ctz-cx*stz
						cx=cxt
						pxt=px*ctz+py*stz
						py=py*ctz-px*stz
						px=pxt
					endif

c					in this case:rotation/x axis
					if(xang(indre).ne.0.)then
						yt=y*ctx+z*stx
						z=z*ctx-y*stx
						y=yt
						cyt=cy*ctx+cz*stx
						cz=cz*ctx-cy*stx
						cy=cyt
						pyt=py*ctz+pz*stz
						pz=pz*ctz-py*stz
						py=pyt
					endif

c					in this case:rotation/y axis
					if(yang(indre).ne.0.)then
						zt=z*cty+x*sty
						x=x*cty-z*sty
						z=zt
						czt=cz*cty+cx*sty
						cx=cx*cty-cz*sty
						cz=czt
						pzt=pz*ctz+px*stz
						px=px*ctz-pz*stz
						pz=pzt
					endif
					ver(i,1)=x
					ver(i,2)=y
					ver(i,3)=z
					ver(i,4)=cx
					ver(i,5)=cy
					ver(i,6)=cz
					ver(i,7)=vea(i,indep,7)
					ver(i,8)=vea(i,indep,8)
					ver(i,9)=vea(i,indep,9)
					ver(i,10)=px
					ver(i,11)=py
					ver(i,12)=pz
					do 2 iq=1,12
2					veri(i,indre,iq)=ver(i,iq)
1				continue
				return
				end



				subroutine unrot
c*************************************************************************
c                                                                        *
c                       u n r o t                                        *
c                                                                        *
c*************************************************************************
c				at each end-plane convert relative coord. in absolute one
c					to fill the absolute coord. buffer and the plot buffer
c				input from subroutine rot:
c					x0,y0,z0,stx,ctx,stz,ctz,sty,cty=rotation matrix elements
c				input:
c					nuve:number of input vector(s) defined
c					indre=current region   index
c					indep=   "    end-plane "
c					ver(iv,ic):description of the particules in the current region
c					iv<=nuve<=maxve:vector index
c					ic= 1:x coord.relative to the current region(in relative frame)
c							2:y            "
c							3:z            "
c							4:x component of the normalized mommentum(in relat.frame)
c							5:y            "
c							6:z            "
c							7:module of the mommentum
c							8:trace length measured from the "spring point"
c							9:life-flag(-1.:buried,0.:just dead,1.:alive)
c							10:x component of the spin (in relat.frame)
c							11:y          "
c							12:z          "
c					xyz0(region index,x y or z):absolute coord. of the rel.origine
c					zang(region index),xang(region index) and yang(region index):
c						rotation angles to define the relat.frame in the absol.frame
c					adata(data index < maxdat,region index):additional data to
c						describe a region
c				output:
c					vea(iv,iep,ic):description of the particules at each end-plane
c					iv<=nuve<=maxve:vector index
c					iep<=nuep<=maxep:end-plane index
c					ic= 1:x coord. in absolute frame
c							2:y            "
c							3:z            "
c							4:x component of the normalized mommentum in absol.frame
c							5:y            "
c							6:z            "
c							7:module of the mommentum
c							8:trace length measured from the "spring point"
c							9:life-flag(-1.:buried,0.:just dead,1.:alive)
c							10 to 15:the same as 1 to 6 but in relative frame coordinates
c							16:x component of the spin
c							17:y          "
c							18:z          "
c							19 to 21: the same as 16 to 18 but in relative frame coord.
c					ilive=   "    number of alive vector(s)
c					nplv<=maxplv:number of trajectories to be plotted
c					nplp(traj.index)<=maxplp:number of submits in the plot-line
c					plott(submit index,traj.index,x y z bx by bz px py pz):
c						absol.coord.of submit, absol. field components and
c						absol. spin components at this point
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /crot/x0,y0,z0,stx,ctx,stz,ctz,sty,cty
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				logical cyl,lfil
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @								,cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				indep=indep+1
				if(indep.gt.maxep)stop'unrot:max num. of end-planes surpas.'
				ilive=0
				do 2 i=1,nuve

c					if the vector is buried:no rotation
					if(ver(i,9).le.-1.)then
						do 5 j=1,21
5						vea(i,indep,j)=0.
						vea(i,indep,9)=-1.
						goto 2
					endif
					x=ver(i,1)
					y=ver(i,2)
					z=ver(i,3)
					cx=ver(i,4)
					cy=ver(i,5)
					cz=ver(i,6)
					px=ver(i,10)
					py=ver(i,11)
					pz=ver(i,12)

c					first:rotations if necessary:
					if(yang(indre).ne.0.)then

c						in this case:rotation/y axis
						zt=z*cty-x*sty
						x=x*cty+z*sty
						z=zt
						czt=cz*cty-cx*sty
						cx=cx*cty+cz*sty
						cz=czt
						pzt=pz*cty-px*sty
						px=px*cty+pz*sty
						pz=pzt
					endif

c					in this case:rotation/x axis
					if(xang(indre).ne.0.)then
						yt=y*ctx-z*stx
						z=z*ctx+y*stx
						y=yt
						cyt=cy*ctx-cz*stx
						cz=cz*ctx+cy*stx
						cy=cyt
						pyt=py*ctx-pz*stx
						pz=pz*ctx+py*stx
						py=pyt
					endif

c					in this case:rotation/z axis
					if(zang(indre).ne.0.)then
						xt=x*ctz-y*stz
						y=y*ctz+x*stz
						x=xt
						cxt=cx*ctz-cy*stz
						cy=cy*ctz+cx*stz
						cx=cxt
						pxt=px*ctz-py*stz
						py=py*ctz+px*stz
						px=pxt
					endif

c					second:translation:
					vea(i,indep,1)=x+x0
					vea(i,indep,2)=y+y0
					vea(i,indep,3)=z+z0
					vea(i,indep,4)=cx
					vea(i,indep,5)=cy
					vea(i,indep,6)=cz
					vea(i,indep,7)=ver(i,7)
					vea(i,indep,8)=ver(i,8)
					vea(i,indep,9)=ver(i,9)
					if(ver(i,9).gt.0.)ilive=ilive+1
					vea(i,indep,10)=ver(i,1)
					vea(i,indep,11)=ver(i,2)
					vea(i,indep,12)=ver(i,3)
					vea(i,indep,13)=ver(i,4)
					vea(i,indep,14)=ver(i,5)
					vea(i,indep,15)=ver(i,6)
					vea(i,indep,16)=px
					vea(i,indep,17)=py
					vea(i,indep,18)=pz
					vea(i,indep,19)=ver(i,10)
					vea(i,indep,20)=ver(i,11)
					vea(i,indep,21)=ver(i,12)

c					plot these points to end properly the previous plot:
c						(the last field seen in rungk (at landing time) 
c						has been saved in bplot(j,i))
					if(i.le.maxplv)then
						call plotin(i,vea(i,indep,1)
     @					,vea(i,indep,2)
     @					,vea(i,indep,3)
     @					,bplot(1,i)
     @					,bplot(2,i)
     @					,bplot(3,i)
     @					,vea(i,indep,16)
     @					,vea(i,indep,17)
     @					,vea(i,indep,18)
     @					,.false.)

c						reset bplot:
						do4j=1,3
4						bplot(j,i)=0.
					endif
2				continue
				return
				end



				subroutine plorot(xr,yr,zr,xa,ya,za,xb,yb,zb)
c*************************************************************************
c                                                                        *
c                       p l o r o t                                      *
c                                                                        *
c*************************************************************************
c				convert the relative coortdinates of a point (xr,yr,zr)
c					in the absolutes one (xa,ya,za)
c				input:
c					xr,yr,zr
c					xyz0(region index,x y or z):absolute coord. of the rel.origine
c					zang(region index),xang(region index) and yang(region index):
c					rotation angles to define the relat.frame in the absol.frame
c				output:
c					xa,ya,za
				include 'sSnake.inc'
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				data indrep/0/
				if(indre.ne.indrep)then
					indrep=indre
					x0=xyz0(indre,1)
					y0=xyz0(indre,2)
					z0=xyz0(indre,3)
					stz=sin(zang(indre))
					ctz=cos(zang(indre))
					stx=sin(xang(indre))
					ctx=cos(xang(indre))
					sty=sin(yang(indre))
					cty=cos(yang(indre))
				endif
				x=xr
				y=yr
				z=zr

c				in this case:rotation/y axis
				if(yang(indre).ne.0.)then
					zt=z*cty-x*sty
					x=x*cty+z*sty
					z=zt
				endif

c				in this case:rotation/x axis
				if(xang(indre).ne.0.)then
					yt=y*ctx-z*stx
					z=z*ctx+y*stx
					y=yt
				endif

c				in this case:rotation/z axis
				if(zang(indre).ne.0.)then
					xt=x*ctz-y*stz
					y=y*ctz+x*stz
					x=xt
				endif
				xa=x+x0
				ya=y+y0
				za=z+z0
				xb=x
				yb=y
				zb=z
				return
				end



				subroutine unplorot(xa,ya,za,xr,yr,zr,xb,yb,zb)
c*************************************************************************
c                                                                        *
c                       u n p l o r o t                                  *
c                                                                        *
c*************************************************************************
c				convert the absolute coortdinates of a point (xr,yr,zr)
c					into the relative one (xa,ya,za)
c				input:
c					xa,ya,za
c					xyz0(region index,x y or z):absolute coord. of the rel.origine
c					zang(region index),xang(region index) and yang(region index):
c					rotation angles to define the relat.frame in the absol.frame
c				output:
c					xr,yr,zr
				include 'sSnake.inc'
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				data indrep/0/
				if(indre.ne.indrep)then
					indrep=indre
					x0=xyz0(indre,1)
					y0=xyz0(indre,2)
					z0=xyz0(indre,3)
					stz=sin(zang(indre))
					ctz=cos(zang(indre))
					stx=sin(xang(indre))
					ctx=cos(xang(indre))
					sty=sin(yang(indre))
					cty=cos(yang(indre))
				endif
				x=xa-x0
				y=ya-y0
				z=za-z0

c				in this case:rotation/z axis
				if(zang(indre).ne.0.)then
					xt=x*ctz+y*stz
					y=y*ctz-x*stz
					x=xt
				endif

c				in this case:rotation/x axis
				if(xang(indre).ne.0.)then
					yt=y*ctx+z*stx
					z=z*ctx-y*stx
					y=yt
				endif

c				in this case:rotation/y axis
				if(yang(indre).ne.0.)then
					zt=z*cty+x*sty
					x=x*cty-z*sty
					z=zt
				endif
				xr=x
				yr=y
				zr=z
				x=xa
				y=ya
				z=za

c				in this case:rotation/z axis
				if(zang(indre).ne.0.)then
					xt=x*ctz+y*stz
					y=y*ctz-x*stz
					x=xt
				endif 

c				in this case:rotation/x axis
				if(xang(indre).ne.0.)then
					yt=y*ctx+z*stx
					z=z*ctx-y*stx
					y=yt
				endif

c				in this case:rotation/y axis
				if(yang(indre).ne.0.)then
					zt=z*cty+x*sty
					x=x*cty-z*sty
					z=zt
				endif
				xb=x
				yb=y
				zb=z
				return
				end



				subroutine arot(xu,yu,zu,xv,yv,zv,ltr)
c*************************************************************************
c                                                                        *
c                       a r o t                                          *
c                                                                        *
c*************************************************************************
c				move the plot data point (xu,yu,zu)
c					by translation (no translation if ltr=.false.)
c					and rotation / z , x and y axis.
c				input:
c					xu,yu,zu
c					xtp,ytp,ztp=translation vector
c					zap,xap,yap=set of angles of rotation in radian
c				output:
c					xv,yv,zv
				include 'sSnake.inc'
				logical ltr
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				if(ltr)then
					x=xu+xtp
					y=yu+ytp
					z=zu+ztp
				else
					x=xu
					y=yu
					z=zu
				endif

c				in this case:rotation/z axis
				if(zap.ne.0.)then
					stz=sin(zap)
					ctz=cos(zap)
					xt=x*ctz-y*stz
					y=y*ctz+x*stz
					x=xt
				endif

c				in this case:rotation/x axis
				if(xap.ne.0.)then
					stx=sin(xap)
					ctx=cos(xap)
					yt=y*ctx-z*stx
					z=z*ctx+y*stx
					y=yt
				endif

c				in this case:rotation/y axis
				if(yap.ne.0.)then
					sty=sin(yap)
					cty=cos(yap)
					zt=z*cty-x*sty
					x=x*cty+z*sty
					z=zt
				endif
				xv=x
				yv=y
				zv=z
				return
				end



				subroutine add_field
c*************************************************************************
c                                                                        *
c                       a d d _ f i e l d                                *
c                                                                        *
c*************************************************************************
				include 'sSnake.inc'
				common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				character*50 form
				character*1 sym
				common /chfield/sym(3),form
				common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @							,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin
     @							,zmax,zstep,jmem
				logical cyl,lfil
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre)
     @								,fact(maxre),cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				dimension indre_arg(10)
				character nomreg*8,sym_cour*1
				dimension nomreg(10),sym_cour(3),xyz(3),f(3)
				dimension b1(maxx,maxy,maxz,3),ylis_cour(maxy)
				logical lmap,ladd,eject

c				extract the name and index of the argument regions:
				i1=2
				i2=i1
				call nextcoma(i2,fifi(indre))
				read(fifi(indre)(i1:i2-1),*)ibarg
				if(ibarg.gt.10)stop'add_field: too many regions'
				do ib=1,ibarg
					i1=i2+1
					i2=i1
					call nextcoma(i2,fifi(indre))
					read(fifi(indre)(i1:i2-1),'(a)')nomreg(ib)
					do 5 i=1,nure
						if(nomreg(ib).eq.rname(i))then
							isave=i
							goto6
						endif
5					continue
					write(6,*)' can''t find this region :',nomreg(ib)
					stop
6					indre_arg(ib)=isave
					if(fifi(indre_arg(ib))(1:1).eq.'#')then
						write(6,*)'add_field#1: region ',nomreg(ib)
     @ 				,'is already an addition type region!'
						stop
					endif
					if(method(indre_arg(ib)).ne.3)then
						write(6,*)'add_field: region ',nomreg(ib)
     @ 				,'is not a field distribution type region!'
						stop
					endif
				enddo
				i1=i2+1
				i2=i1
				call nextcoma(i2,fifi(indre))

c				extract the symmmetry word and save the current region data:
				read(fifi(indre)(i1:i2-1),'(3a1)')sym_cour
				do l=1,3
					if(sym_cour(l).ne.'s'.and.sym_cour(l).ne.'a'
     @			.and.sym_cour(l).ne.'n'
     @			.and.sym_cour(l).ne.'S'.and.sym_cour(l).ne.'A'
     @			.and.sym_cour(l).ne.'N')
     @			stop'add_field: invalid symmetry word'
				enddo
				if(cyl(indre))
     @		stop'add_field: cylindrical map not yet implemented'
				indre_cour=indre
				ndx_cour=1
				xmin_cour=0.
				xmax_cour=0.
				xstep_cour=0.
				ndy_cour=1
				ymin_cour=0.
				ymax_cour=0.
				ystep_cour=0.
				ndz_cour=1
				zmin_cour=0.
				zmax_cour=0.
				zstep_cour=0.
				goto(801,802,803,804,805,806)indfi(indre)

c				1 dim. Poisson type map
801			ndy_cour=nint(adata(1,indre))
				ymin_cour=adata(2,indre)
				ymax_cour=adata(3,indre)
				ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
				goto 100

c				2 dim. measured type map, uniform y steps
802			ndx_cour=nint(adata(1,indre))
				xmin_cour=adata(2,indre)
				xmax_cour=adata(3,indre)
				ndy_cour=nint(adata(4,indre))
				ymin_cour=adata(5,indre)
				ymax_cour=adata(6,indre)
				xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
				ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
				goto 100

c				2 dim. measured type map, variable y steps
803			ndx_cour=nint(adata(1,indre))
				xmin_cour=adata(2,indre)
				xmax_cour=adata(3,indre)
				ndy_cour=nint(adata(4,indre))
				do j=1,ndy_cour
					ylis_cour(j)=adata(j+4,indre)
				enddo
				xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
				goto 100

c				2 dim. Poisson type map
804			ndy_cour=nint(adata(1,indre))
				ymin_cour=adata(2,indre)
				ymax_cour=adata(3,indre)
				ndz_cour=nint(adata(4,indre))
				zmin_cour=adata(5,indre)
				zmax_cour=adata(6,indre)
				ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
				zstep_cour=(zmax_cour-zmin_cour)/real(ndz_cour-1)
				goto 100

c				3 dim. Poisson type map
805			continue
806			ndx_cour=nint(adata(1,indre))
				xmin_cour=adata(2,indre)
				xmax_cour=adata(3,indre)
				ndy_cour=nint(adata(4,indre))
				ymin_cour=adata(5,indre)
				ymax_cour=adata(6,indre)
				ndz_cour=nint(adata(7,indre))
				zmin_cour=adata(8,indre)
				zmax_cour=adata(9,indre)
				xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
				ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
				zstep_cour=(zmax_cour-zmin_cour)/real(ndz_cour-1)
				goto 100

100			if(ndx_cour.gt.maxx)
     @		stop'ifield:max. field points in x surpassed'
				if(ndy_cour.gt.maxy)
     @		stop'ifield:max. field points in y surpassed'
				if(ndz_cour.gt.maxz)
     @		stop'ifield:max. field points in z surpassed'

c				build the field map of the current region by adding the field
c					distribution of the two argument regions
				do i=1,ndx_cour
					do j=1,ndy_cour
						do k=1,ndz_cour
							do l=1,3
								b1(i,j,k,l)=0.
							enddo
						enddo
					enddo
				enddo

c				in case of diag. in field:
				indve=1
				do ib=1,ibarg
					indre=indre_arg(ib)
					lmap=.true.
					call ifield(lmap,ladd)
					if(ladd)then
						write(6,*)'add_field#2: region ',nomreg(ib)
     @				,'is already an addition type region!'
						stop
					endif
					do i=1,ndx_cour
						x_cour=xmin_cour+xstep_cour*real(i-1)
						do j=1,ndy_cour
							y_cour=ymin_cour+ystep_cour*real(j-1)
							do k=1,ndz_cour
								z_cour=zmin_cour+zstep_cour*real(k-1)
								indre=indre_cour
								call plorot(x_cour,y_cour,z_cour
     @						,x_abs,y_abs,z_abs,xb,yb,zb)
								indre=indre_arg(ib)
								call unplorot(x_abs,y_abs,z_abs
     @						,xyz(1),xyz(2),xyz(3),xb,yb,zb)

c								don't eject if outside the free box:
								call field(eject,xyz,f
     @						,.false.,.false.,.false.,.false.)

c								no contribution to the final field in case of eject:
								if(.not.eject)then
									call plorot(f(1),f(2),f(3)
     @							,xb,yb,zb,bx_abs,by_abs,bz_abs)
									indre=indre_cour
									call unplorot(bx_abs,by_abs,bz_abs
     @							,xb,yb,zb,bx_cour,by_cour,bz_cour)
									b1(i,j,k,1)=b1(i,j,k,1)+bx_cour
									b1(i,j,k,2)=b1(i,j,k,2)+by_cour
									b1(i,j,k,3)=b1(i,j,k,3)+bz_cour
								endif
							enddo
						enddo
					enddo
				enddo

c				initialize the current region field map:
				indre=indre_cour
				do l=1,3
					sym(l)=sym_cour(l)
				enddo
				ndx=ndx_cour
				xmin=xmin_cour
				xmax=xmax_cour
				xstep=xstep_cour
				ndy=ndy_cour
				ymin=ymin_cour
				ymax=ymax_cour
				ystep=ystep_cour
				ndz=ndz_cour
				zmin=zmin_cour
				zmax=zmax_cour
				zstep=zstep_cour
				jmem=2
				do i=1,ndx
					do j=1,ndy
						if(indfi(indre).eq.3)ylis(j)=ylis_cour(j)
						do k=1,ndz
							do l=1,3
								b(i,j,k,l)=b1(i,j,k,l)
							enddo
						enddo
					enddo
				enddo

c				1 dim. Poisson type map: compute the gradient
				if(indfi(indre).eq.1)then
					if(ndy.lt.2)stop'add_field: not enough y in a index=1 map'
					db(1)=(b(1,2,1,3)-b(1,1,1,3))/ystep
					do j=2,ndy-1
						db(j)=(b(1,j+1,1,3)-b(1,j-1,1,3))/(2*ystep)
					enddo
					db(ndy)=(b(1,ndy,1,3)-b(1,ndy-1,1,3))/ystep
				endif
				return
				end



				subroutine ifield(lmap,ladd)
c*************************************************************************
c                                                                        *
c                       i f i e l d                                      *
c                                                                        *
c*************************************************************************
c				initialize the field in the new region for the subroutine send and/or field
c					(input and output list depend on methode,type and field index)
c				parameters:
c					maxx:maximum number of x values for the field
c					maxy:   "      "    "  y  "      "   "   "
c					maxz:   "      "    "  z  "      "   "   "
c					maxdat: "      "    "  additional data to describe a region
c				input:
c					indre=current region   index
c					method(region index):methode to use to carry particules in this region
c					itype(region index):type in this method
c					indfi(region index):index in this type and in this method
c					fifi(region index):field file name
c					adata(data index < maxdat,region index):additional data to describe a region
c				output for send: a(6,6):1st order transport type matrix
c				output for subroutine field:
c					indfi(region index): over-read on the field-file
c					sym(x y or z)='s' for symmetry,='a' for antisymmetry,else='n'
c					ndx,xmin,xmax:control param. for x-mesh
c					ndy,ymin,ymax:       "           y  "
c					ndz,zmin,zmax:       "           z  "
c					form:read format for field and/or y values
c					b(xindex,yindex,zindex,bx by or bz):components of the magn.field
c					db(yindex)=y component of the gradient of the field
c					ylis(yindex)=y value for yindex in case of irregular mesh
c				In case of field map to be read on a disk file:
c				lmap=.t. : read the header, the field map and the user plot (if any)
c				lmap=.f. : read only the header because the map and the user plot will be 
c					overwriten by fil.
c   for meth=3 , type=2 , index=1 (dipole with grad. ):
c   for meth=3 , type=2 , index=16+17 (arbitrary field numeric.deriv.)
c   for meth=3 , type=2 , index=2 to 14 and 18 (raytrace library ):
c   for meth=3 , type=2 , index=15 (clam):
				include 'sSnake.inc'
				logical lbox,luser,bare
				character unit*3,name*2
				common /cplot/unit(15),name(15)
				common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @							,nplvold,h1(2,maxq),v1(2,maxq),nq
     @							,ip,ih,iv,kp,bplot(3,maxplv)
     @							,plobox(5,4,2,maxre),nface(maxre),lbox
     @							,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @							,iplpi(maxplv),iplpf(maxplv)
     @							,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @							,luser,nclic,lclic(maxre),clight(4,maxre)
     @							,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @							,nrs(maxre)
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				character*50 form
				character*1 sym
				common /chfield/sym(3),form
				common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @							,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin
     @							,zmax,zstep,jmem
				common/dipgra/rm,alp,bet
				common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @							,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @							,tbound,dst
				common/raytra/fact1(5),qrad1,fact2(5),qrad2,dl
				common/diprt/dipdata(75)
				common/clam1/d0,dx,dy,bdref,ta
				logical cyl,lfil
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre)
     @								,fact(maxre),cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				dimension ad(3),bd(3),cd(3),plotusc(3)
				logical l,lmap,ladd
				common/matrix/a(6,6),dty,s(3,3)
				ladd=.false.
				go to(101,102,103),method(indre)

c*********************no field**************************
101			return

c**********************matrix:**************************
102			goto(301,302),itype(indre)

c**********************1st order matrix*
c				fill the current matrix array acordind to the index value:
301			goto(201,202),indfi(indre)

c				matrix for 1000.mm long drif,parallele faces:
201			dty=1000.
				goto 5

c				matrix for 500.mm long drif,parallele faces:
202			dty=500.
				goto 5
5				continue

c				matrix for drift:
				do 1 i=1,6
				do 1 j=1,6
					a(i,j)=0.
1					s(i,j)=0.

c				x:
				a(1,1)=1.
				s(1,1)=1.
				a(2,1)=dty

c				xp:
				a(2,2)=1.
				s(2,2)=1.

c				z:
				a(3,3)=1.
				s(3,3)=1.
				a(4,3)=dty

c				zp:
				a(4,4)=1.

c				l:
				a(5,5)=1.

c				d:
				a(6,6)=1.

				return

c******************2nd order matrix*
302			continue
				stop'ifield:not yet implemented'

c************************rung-khuta************************
103			continue

c				initialize the subroutine "field"
				goto(501,502,503),itype(indre)

c				itype=1:homogenous field
501			goto(401,402,403),indfi(indre)

c				indfi=1:pure dipolar field given in data:
401			fc(1)=adata(1,indre)
				fc(2)=adata(2,indre)
				fc(3)=adata(3,indre)
				return

c				indfi=2:pure quadrupolar field given in data:
c					( fc(1)>0. for the field of a qpole having x as axis,focusing positive
c						particles in y )
c					( fc(2)>0. for the field of a qpole having y as axis,focusing positive
c						particles in z )
c					( fc(3)>0. for the field of a qpole having x as axis,focusing positive
c						particles in x )
402			fc(1)=adata(1,indre)
				fc(2)=adata(2,indre)
				fc(3)=adata(3,indre)
				return

c				indfi=3: dipolar field along z, qpole+sextupole having y as axis:
c					( fc(1)=dipolar field, fc(2)>0.=qpolar gradiant focusing in z,
c					fc(3)>0.=sextupolar grad. defocusing in x, no effect in x,
c					fc(4)>0.=sextupolar grad. focusing in z, no effect in x )
403			fc(1)=adata(1,indre)
				fc(2)=adata(2,indre)
				fc(3)=adata(3,indre)
				fc(4)=adata(4,indre)
				return

c				itype=2:analytical field
502			goto(901,902,903,904,905,906,907,908,909,910,911,912,913,914
     @		,915,916,917,918,919,920),indfi(indre)

c				indfi=1;field with gradient:
901			rm=adata(1,indre)
				alp=adata(2,indre)
				bet=adata(3,indre)
				return

c				indfi=2:unit gradient (1 t/mm ) along x and z (qpole axis = y )
c					nothing to do
902			return

c				indfi=3:quadrupole
903			fact1(1)=adata(1,indre)
				qrad1=adata(2,indre)
				return

c				indfi=4:dipole entrance fringing field:
904			fact1(1)=adata(1,indre)
				qrad1=adata(2,indre)
				return

c				indfi=5:dipole exit fringing field:
905			fact1(1)=adata(1,indre)
				qrad1=adata(2,indre)
				return

c				indfi=6:uniform dipole :
906			fact1(1)=adata(1,indre)
				fact1(2)=adata(2,indre)
				fact1(3)=adata(3,indre)
				return

c				indfi=7:multipole entrance fringing field:
907			fact1(1)=adata(1,indre)
				fact1(2)=-adata(2,indre)
				fact1(3)=adata(3,indre)
				fact1(4)=-adata(4,indre)
				fact1(5)=adata(5,indre)
				qrad1=adata(6,indre)
				return

c				indfi=8:uniform multipole :
908			fact1(1)=adata(1,indre)
				fact1(2)=-adata(2,indre)
				fact1(3)=adata(3,indre)
				fact1(4)=-adata(4,indre)
				fact1(5)=adata(5,indre)
				qrad1=adata(6,indre)
				return

c				indfi=9:multipole exit fringing field:
909			fact1(1)=adata(1,indre)
				fact1(2)=-adata(2,indre)
				fact1(3)=adata(3,indre)
				fact1(4)=-adata(4,indre)
				fact1(5)=adata(5,indre)
				qrad1=adata(6,indre)
				return

c				indfi=10:qq overlaping fields:
910			fact1(1)=adata(1,indre)
				fact1(2)=-adata(2,indre)
				fact1(3)=adata(3,indre)
				fact1(4)=-adata(4,indre)
				fact1(5)=adata(5,indre)
				qrad1=adata(6,indre)
				fact2(1)=adata(7,indre)
				fact2(2)=-adata(8,indre)
				fact2(3)=adata(9,indre)
				fact2(4)=-adata(10,indre)
				fact2(5)=adata(11,indre)
				qrad2=adata(12,indre)
				dl=adata(13,indre)
				return

c				indfi=11:qd overlaping fields:
911			fact1(1)=adata(1,indre)
				fact1(2)=-adata(2,indre)
				fact1(3)=adata(3,indre)
				fact1(4)=-adata(4,indre)
				fact1(5)=adata(5,indre)
				qrad1=adata(6,indre)
				fact2(1)=adata(7,indre)
				fact2(2)=adata(8,indre)
				qrad2=adata(9,indre)
				dl=adata(10,indre)
				return

c				indfi=12:dd overlaping fields:
912			fact1(1)=adata(1,indre)
				fact1(2)=adata(2,indre)
				qrad1=adata(3,indre)
				fact2(1)=adata(4,indre)
				fact2(2)=adata(5,indre)
				qrad2=adata(6,indre)
				dl=adata(7,indre)
				return

c				indfi=13:dq overlaping fields (old version, PV 4/9/97):
913			fact1(1)=adata(1,indre)
				fact1(2)=adata(2,indre)
				qrad1=adata(3,indre)
				fact2(1)=adata(4,indre)
				fact2(2)=-adata(5,indre)
				fact2(3)=adata(6,indre)
				fact2(4)=-adata(7,indre)
				fact2(5)=adata(8,indre)
				qrad2=adata(9,indre)
				dl=adata(10,indre)
				return

c				indfi=14:dq overlaping fields (new version, PV 4/9/97):
914			fact1(1)=adata(1,indre)
				fact1(2)=adata(2,indre)
				qrad1=adata(3,indre)
				fact2(1)=adata(4,indre)
				fact2(2)=-adata(5,indre)
				fact2(3)=adata(6,indre)
				fact2(4)=-adata(7,indre)
				fact2(5)=adata(8,indre)
				qrad2=adata(9,indre)
				dl=adata(10,indre)
				open(3,file=fifi(indre),status='old')
				read(3,107)(dipdata(j),j=1,5),(dipdata(j),j=11,22 )
     @		,( dipdata(j ) , j=25,64)
				close(3)
				return

c				indfi=15:clam central field:
915			bref=adata(1,indre)
				dref=adata(2,indre)
				bdref=bref*dref
				beta=adata(3,indre)*dtor
				alpha=adata(4,indre)
				ta=tan(alpha)
				dx=cos(angt)
				dy=sin(angt)
				do 23 i=1,3
					ad(i)=xyz0(indre,i)
					bd(i)=0.
23			continue
				cd(1)=-sin(beta)
				cd(2)= cos(beta)
				cd(3)=0.
				call dist(ad,bd,cd,d02,d0)
				return

c				indfi=16;field with gradient , numeric. derivation:
916			rms=adata(1,indre)
				alps=adata(2,indre)
				bets=adata(3,indre)
				xc(1)=adata(4,indre)
				yc(1)=adata(5,indre)
				r0(1)=adata(6,indre)
				e0(1)=adata(7,indre)
				s0(1)=adata(8,indre)
				s1(1)=adata(9,indre)
				s2(1)=adata(10,indre)
				s3(1)=adata(11,indre)
				xc(2)=adata(12,indre)
				yc(2)=adata(13,indre)
				r0(2)=adata(14,indre)
				e0(2)=adata(15,indre)
				s0(2)=adata(16,indre)
				s1(2)=adata(17,indre)
				s2(2)=adata(18,indre)
				s3(2)=adata(19,indre)
				tbound=adata(20,indre)*dtor
				dst=adata(21,indre)
				luser=.true.
				return

c				indfi=17;dipole field from raytrace, old fashion:
917			open(3,file=fifi(indre),status='old')
				read (3,107) ( dipdata( j ) , j=1,5 ), (dipdata( j ), j=11,22)
     @		,( dipdata( j ) , j=25,64)
107			format(5f10.5/ 5f10.5/3f10.5/4f10.5/ 4f10.5/ 6f10.5/ 6f10.5
     @				/6f10.5/ 4f10.5/ 7f10.5/ 7f10.5)
				close(3)
				return

c				indfi=18; dipole field from raytrace, new fashion:
918			open(3,file=fifi(indre),status='old')
				read (3,107) ( dipdata( j ) , j=1,5 ), (dipdata( j ), j=11,22)
     @		,( dipdata( j ) , j=25,64)
				close(3)
				return

c				indfi=19:helmotz coil
919			fact1(1)=adata(1,indre)
				qrad1=adata(2,indre)
				return

c				indfi=20: Solenoid
920			fact1(1)=adata(1,indre)
				fact1(2)=adata(2,indre)
				fact1(3)=adata(3,indre)
				return

c				itype=3:read the file:fifi(indre)
503			continue

c				built a new map by adding two existing field distributions?
				if(fifi(indre)(1:1).eq.'#')then
					ladd=.true.
					return
				endif
				open(3,file=fifi(indre),status='old')
				write(6,'(2a)')'reading map file: ',fifi(indre)

c				overread the interpolation method type to be used:
				read(3,*)indfi(indre)

c				read the symetries of the field:
				read(3,'(3a1)')sym

c				next reads depend on indfi:
				goto(801,802,803,804,805,806),indfi(indre)

c				indfi=1:1 dim. poisson type table(bz(0,y,0) and dbz/dy(0,y,0))
c					assume:bz(x,y,z)=bz(0,y,z) and bz(x,y,-z)=bz(x,y,z)
801			read(3,*)ndy,ymin,ymax
				if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					ymin=ymin*dtor
					ymax=ymax*dtor
				endif
				ystep=(ymax-ymin)/real(ndy-1)
				q1=ymin+ystep
				q2=ymax
				f1=frbox(1,indre,2)
				a1=abs(f1)
				f2=frbox(2,indre,2)
				a2=abs(f2)
				l=sym(2).eq.'N'.or.sym(2).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in y'

c				read the format:
				read(3,'(a)')form

c				read the field and gradian:
				if(lfil)then
					ndy=ndy+2
					ymin=ymin-ystep
					ymax=ymax+ystep
				endif
				ndx=1
				xmin=0.
				xmax=0.
				xstep=0.
				ndz=1
				zmin=0.
				zmax=0.
				zstep=0.
				if(.not.lmap)then
					close(3)
					return
				endif
				read(3,form)(b(1,i,1,3),db(i),i=1,ndy)

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					do 811 i=1,ndy
811				db(i)=db(i)*rtod
				endif
				goto 100

c				indfi=2:2 dim. measured type table(bz(x,y,0))constant y increment
c					assume: bz(x,y,-z)=bz(x,y,z)
802			read(3,*)ndx,xmin,xmax
				if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
				xstep=(xmax-xmin)/real(ndx-1)
				q1=xmin+xstep
				q2=xmax
				f1=frbox(1,indre,1)
				a1=abs(f1)
				f2=frbox(2,indre,1)
				a2=abs(f2)
				l=sym(1).eq.'N'.or.sym(1).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in x'
				read(3,*)ndy,ymin,ymax
				if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					ymin=ymin*dtor
					ymax=ymax*dtor
				endif
				ystep=(ymax-ymin)/real(ndy-1)
				q1=ymin+ystep
				q2=ymax
				f1=frbox(1,indre,2)
				a1=abs(f1)
				f2=frbox(2,indre,2)
				a2=abs(f2)
				l=sym(2).eq.'N'.or.sym(2).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in y'
				ndz=1
				zmin=0.
				zmax=0.
				zstep=0.

c				read the format:
				read(3,'(a)')form

c				read the field
				if(.not.lmap)then
					close(3)
					return
				endif
				read(3,form)((b(i,j,1,3),i=1,ndx),j=1,ndy)
				goto 100

c				indfi=3:2 dim. measured type table(bz(x,y,0))variable y increment
c					assume: bz(x,y,-z)=bz(x,y,z)
803				read(3,*)ndx,xmin,xmax
				if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
				xstep=(xmax-xmin)/real(ndx-1)
				q1=xmin+xstep
				q2=xmax
				f1=frbox(1,indre,1)
				a1=abs(f1)
				f2=frbox(2,indre,1)
				a2=abs(f2)
				l=sym(1).eq.'N'.or.sym(1).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in x'
				read(3,*)ndy
				if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
				ndz=1
				zmin=0.
				zmax=0.
				zstep=0.

c				read the format:
				read(3,'(a)')form

c				read the field
				read(3,form)(ylis(j),j=1,ndy)

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					do 813 j=1,ndy 
813				ylis(j)=ylis(j)*dtor
				endif
				if(.not.lmap)then
					close(3)
					return
				endif
				read(3,form)((b(i,j,1,3),i=1,ndx),j=1,ndy)
				q1=ylis(2)
				q2=ylis(ndy)
				f1=frbox(1,indre,2)
				a1=abs(f1)
				f2=frbox(2,indre,2)
				a2=abs(f2)
				l=sym(2).eq.'N'.or.sym(2).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in y'
				jmem=2
				goto 100

c				indfi=4:2 dim. poisson type table(by(0,y,z) and bz(0,y,z))
c					assume: bz(x,y,z)=bz(0,y,z)
804			read(3,*)ndy,ymin,ymax
				if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					ymin=ymin*dtor
					ymax=ymax*dtor
				endif
				ystep=(ymax-ymin)/real(ndy-1)
				q1=ymin+ystep
				q2=ymax
				f1=frbox(1,indre,2)
				a1=abs(f1)
				f2=frbox(2,indre,2)
				a2=abs(f2)
				l=sym(2).eq.'N'.or.sym(2).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in y'
				read(3,*)ndz,zmin,zmax
				if(ndz.gt.maxz)stop'ifield:max. field points in z surpassed'
				zstep=(zmax-zmin)/real(ndz-1)
				q1=zmin+zstep
				q2=zmax
				f1=frbox(1,indre,3)
				a1=abs(f1)
				f2=frbox(2,indre,3)
				a2=abs(f2)
				l=sym(3).eq.'N'.or.sym(3).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in z'
				ndx=1
				xmin=0.
				xmax=0.
				xstep=0.

c				read the format:
				read(3,'(a)')form

c				read the field
				if(.not.lmap)then
					close(3)
					return
				endif
				read(3,form)((b(1,j,k,2),b(1,j,k,3),j=1,ndy),k=1,ndz)
				goto 100

c				indfi=5:3 dim table
805			continue

c				index 6: same initialization as index 5:
806			read(3,*)ndx,xmin,xmax
				if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
				xstep=(xmax-xmin)/real(ndx-1)
				if(indfi(indre).eq.5)then
					q1=xmin+xstep
				else
					q1=xmin
				endif
				q2=xmax
				f1=frbox(1,indre,1)
				a1=abs(f1)
				f2=frbox(2,indre,1)
				a2=abs(f2)
				l=sym(1).eq.'N'.or.sym(1).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in x'
				read(3,*)ndy,ymin,ymax
				if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'

c				if cylindrical , convert degree to radian :
				if(cyl(indre))then
					ymin=ymin*dtor
					ymax=ymax*dtor
				endif
				ystep=(ymax-ymin)/real(ndy-1)
				if(indfi(indre).eq.5)then
					q1=ymin+ystep
				else
					q1=ymin
				endif
				q2=ymax
				f1=frbox(1,indre,2)
				a1=abs(f1)
				f2=frbox(2,indre,2)
				a2=abs(f2)
				l=sym(2).eq.'N'.or.sym(2).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in y'
				read(3,*)ndz,zmin,zmax
				if(ndz.gt.maxz)stop'ifield:max. field points in z surpassed'
				zstep=(zmax-zmin)/real(ndz-1)
				if(indfi(indre).eq.5)then
					q1=zmin+zstep
				else
					q1=zmin
				endif
				q2=zmax
				f1=frbox(1,indre,3)
				a1=abs(f1)
				f2=frbox(2,indre,3)
				a2=abs(f2)
				l=sym(3).eq.'N'.or.sym(3).eq.'n'
				if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l
     @		.and.(a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @		write(6,*)' warning : field < free-box in z'

c				read the format:
				read(3,'(a)')form
				if(.not.lmap)then
					close(3)
					return
				endif

c				read the field
				if(form(1:3).eq.'bin')then
					read(3)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3),i=1,ndx),j=1,ndy)
     @			,k=1,ndz)
				else
					read(3,form)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3),i=1,ndx),j=1
     @			,ndy),k=1,ndz)
				endif
				goto100

c				read the user plot asuming that it is in relat. coord.,
c					convert it in absol. coord. and add it to previous data in /plot/:
100			read(3,*,end=200)inulu
				do 300 i=1,inulu
					nulu=nulu+1
					if(nulu.gt.nulumax)then
						nulu=nulumax
						write(6,*)' warning ifield: too many user line'
						goto 200
					endif
					luser=.true.
					read(3,*)nuseuc
					nuseu(nulu)=nuseuc
					if(nuseu(nulu).gt.nuseumax)nuseu(nulu)=nuseumax
					do 400 is=1,nuseuc
						read(3,'(3e12.5)')plotusc
						if(is.le.nuseumax)then
							call plorot(plotusc(1),plotusc(2),plotusc(3)
     @					,plotus(is,nulu,1),plotus(is,nulu,2),plotus(is,nulu,3)
     @					,xb,yb,zb)
						endif
400				continue
300			continue
200			close(3)
				write(6,'(a)')'done'
				return
				end



				subroutine field(eject,xyz,f,logx,logy,logz,lfreeb)
c*************************************************************************
c                                                                        *
c                       f i e l d                                        *
c                                                                        *
c*************************************************************************
c				find the 3 components of the field f(3) at the point xyz(3)
c				field is called by rungk  .true. output value for eject causes
c					the death of the particle.
c				input:
c					xyz(x y or z):relative coord. of the point
c				input for subroutine ifield:
c					fact(region index):multipl.factor of the field or ref.mommentum
c					indfi(region index): over-read on the field-file(if any)
c					sym(x y or z)='s' for symmetry,='a' for antisymmetry,else='n'
c					ndx,xmin,xmax:control param. for x-mesh (if any)
c					ndy,ymin,ymax:       "           y  "      "
c					ndz,zmin,zmax:       "           z  "      "
c					form:read format for field (b),gradient (db) and y values (ylis)
c					b(xindex,yindex,zindex,bx by or bz):components of the magn.field
c					db(yindex)=y component of the gradient of the field
c					ylis(yindex)=y value for yindex in case of irregular mesh
c				output:
c					f(x y or z):field component along x y or z
c					eject:value .true. cause the death of the particule
				include 'sSnake.inc'
				character diag*100,diagt*100
				common/cdiag/diag(maxve),diagt
				common/ndiag/iepdead(maxve),iredead(maxve)
				character mdiag*3,qdiag(6)*1
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ver(maxve,12),xyz0(maxre,3)
     @							,frbox(2,maxre,3)
     @							,zang(maxre),xang(maxre),yang(maxre)
     @							,indve,indre,indep,ilive,adata(maxdat,maxre)
     @							,indrer(maxre),tept(maxepr,maxre)
     @							,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @							,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @							,colli(maxepr,maxre,6)
				character*50 form
				character*1 sym
				common /chfield/sym(3),form
				common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @		,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin,zmax,zstep
     @		,jmem
				logical cyl,lfil
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @					,cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				dimension f(3),xyz(3),xyzc(3)
				logical eject,logx,logy,logz,ok,lfreeb
				data qdiag/'x','y','z','r','t','z'/

c				convert to cylindrical coord. ?
				if(cyl(indre))then

c					save the cartesian coord. in xyzc :
					do 13 i=1,3
13				xyzc(i)=xyz(i)
					xyz(1)=sqrt(xyzc(1)*xyzc(1)+xyzc(2)*xyzc(2))
					xyz(2)=atan2(xyzc(2),xyzc(1))
				endif

c				in the free-box?
				eject=.false.

c				no free box test in case of field addition:
				if(lfreeb)then
					do 11 k=1,3
						if(xyz(k).lt.frbox(1,indre,k))then
							idiag=k
							mdiag='min'
							eject=.true.
						endif
						if(xyz(k).gt.frbox(2,indre,k))then
							idiag=k
							mdiag='max'
							eject=.true.
						endif
11				continue
					if(cyl(indre))idiag=idiag+3

c					built the diag:
					if(eject)then
						write(diagt,'(4a)')'exits free box in :'
     @				,qdiag(idiag),mdiag,' during step by step ray tracing'
						return
					endif
				endif

c				get the field acording to its type:
				goto(601,602,603)itype(indre)

c				itype=1:homogenous field
601			goto(611,612,613)indfi(indre)
611			f(1)=fc(1)
				f(2)=fc(2)
				f(3)=fc(3)
				goto 1000
612			f(1)=fc(2)*xyz(3)+fc(3)*xyz(2)
				f(2)=fc(1)*xyz(3)+fc(3)*xyz(1)
				f(3)=fc(2)*xyz(1)+fc(1)*xyz(2)
				goto 1000
613			f(1)=fc(2)*xyz(3)+fc(3)*xyz(1)*xyz(3)
     @			+fc(4)*(xyz(3)**2-xyz(1)**2)*.5
				f(2)=0.
				f(3)=fc(1)+fc(2)*xyz(1)+fc(3)*(xyz(1)**2-xyz(3)**2)*.5
     @			+fc(4)*xyz(1)*xyz(3)
				goto 1000

c				itype=2:analytical field
602			call analyf(eject,xyz,f,indfi(indre),indve)
				goto 1000

c				itype=3:possible symetry and interpolation
603			continue

c				1/range reduction due to symetry
				x=xyz(1)
				y=xyz(2)
				z=xyz(3)
				if(sym(1).ne.'N'.and.sym(1).ne.'n'.and.x.lt.0.)x=-xyz(1)
				if(sym(2).ne.'N'.and.sym(2).ne.'n'.and.y.lt.0.)y=-xyz(2)
				if(sym(3).ne.'N'.and.sym(3).ne.'n'.and.z.lt.0.)z=-xyz(3)

c				2/interpolation
				goto(701,702,703,704,705,706)indfi(indre)

c				indfi=1:1 dim. poisson type table(bz(0,y,0) and dbz/dy(0,y,0))
c				assume:bz(x,y,z)=bz(0,y,z) and bz(x,y,-z)=bz(x,y,z)
701			continue
				stop'field:not yet implemented'

c				indfi=2:2 dim. measured type table(bz(x,y,0))constant y increment
c					assume: bz(x,y,-z)=bz(x,y,z)
702			continue

c				finds the mesh
				xr=(x-xmin)/xstep
				ixm=int(xr)
				if(logx)then
					if(mod(ixm,2).eq.0)ixm=ixm-1
					ix=ixm+2
					ixp=ix+2
					p=.5*(xr-real(ix-1))
				else
					ix=ixm+1
					ixp=ix+1
					p=xr-real(ixm)
				endif
				eject=eject.or.ix.lt.0
				eject=eject.or.ix.gt.ndx
				yr=(y-ymin)/ystep
				iym=int(yr)
				if(logy)then
					if(mod(iym,2).eq.0)iym=iym-1
					iy=iym+2
					iyp=iy+2
					q=.5*(yr-real(iy-1))
				else
					iy=iym+1
					iyp=iy+1
					q=yr-real(iym)
				endif
				eject=eject.or.iy.lt.0
				eject=eject.or.iy.gt.ndy

c				build the diag:
				if(eject)then
					write(diagt,'(5(a,i3))')'exits field map :x index='
     @			,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy,']'
					return
				endif
				f00= b(ix,iy,1,3)
				f0m= b(ix,iym,1,3)
				f0p= b(ix,iyp,1,3)
				fm0= b(ixm,iy,1,3)
				fp0= b(ixp,iy,1,3)
				fpp= b(ixp,iyp,1,3)
				ok=f00.ne.1.e30.and.f0m.ne.1.e30.and.f0p.ne.1.e30
				ok=ok.and.fm0.ne.1.e30.and.fp0.ne.1.e30.and.fpp.ne.1.e30
				if(ok)then
					ixs=1
					iys=1
					goto322
				elseif(logx.or.logy)then
					eject=.true.
					return
				endif

c				finds the 6-point formula center
				do 331 jj=1,4
					jx=jj/3
					jy=1-mod(jj,2)
					iwx=1+2*jx
					iwy=1+2*jy
					do 332 iax=ix-jx,ix-jx+iwx,iwx
					do 332 iay=iy-jy,iy-jy+iwy,iwy
						ixcen=iax
						iycen=iay
						if((iax-1).lt.1.or.(iax+1).gt.ndx)goto332
						if((iay-1).lt.1.or.(iay+1).gt.ndy)goto332
						do 303 lx=ixcen-1,ixcen+1
						do 303 ly=iycen-1,iycen+1
							if(b(lx,ly,1,3).eq.1.e30)goto332
303					continue
						goto311


332				continue
331			continue
304			eject=.true.

c				build the diag:
				write(diagt,'(a,5(a,i3))')
     @		'can''t find 6 non zero field points'
     @		,':x index=',ix,'[1,',ndx,'] y index=',iy,'[1,',ndy,']'
				return

c				checks the 6th point
311			ixd=ixcen-ix
				iyd=iycen-iy
				ix0=isign(1,-ixd)
				iy0=isign(1,-iyd)
				do 312 jx=ix0,-ix0,-2*ix0
				do 312 jy=iy0,-iy0,-2*iy0
					ixs=jx
					iys=jy
					if((ixcen+jx).lt.1.or.(ixcen+jx).gt.ndx)goto312
					if((iycen+jy).lt.1.or.(iycen+jy).gt.ndy)goto312
					if(b(ixcen+jx,iycen+jy,1,3).eq.1.e30)goto312
					goto 321
312			continue
				goto 304

c				2-dim 2nd order inter/extra-polation for bz (x,y,0)
321			p=(p-real(ixd))*real(ixs)
				q=(q-real(iyd))*real(iys)
				f00= b(ix+ixd     ,iy+iyd     ,1,3)
				f0m= b(ix+ixd     ,iy+iyd-iys ,1,3)
				f0p= b(ix+ixd     ,iy+iyd+iys ,1,3)
				fm0= b(ix+ixd-ixs ,iy+iyd     ,1,3)
				fp0= b(ix+ixd+ixs ,iy+iyd     ,1,3)
				fpp= b(ix+ixd+ixs ,iy+iyd+iys ,1,3)
322			b30= f0m*q*(q-1.) + fm0*p*(p-1.) + fp0*p*(p-2.*q+1.)
     @			+ f0p*q*(q-2.*p+1.)
				b30=.5*b30 + f00*( 1. +p*q - p**2 - q**2) + fpp*p*q

c				1st & 2nd field derivatives  in the (x,y) plane
				ffp1 =  fm0 + fp0 -2.*f00
				ffp2 = -fp0 - f0p + f00 + fpp
				ffq1 =  f0m + f0p -2.*f00
				grx1 = (p*ffp1+q*ffp2+.5*(-fm0+fp0))*real(ixs)/xstep
				gry1 = (p*ffp2+q*ffq1+.5*(-f0m+f0p))*real(iys)/ystep
				grx2 = ffp1/xstep/xstep
				gry2 = ffq1/ystep/ystep
				if(cyl(indre))then
					gry1=gry1/xyz(1)
					xlaplace=grx2+grx1/xyz(1)+gry2/(xyz(1)*xyz(1))
				else
					xlaplace=grx2+gry2
				endif

c				taylor's  serie for z.ne.0.
				f(1) = z*grx1
				f(2) = z*gry1
				f(3) = b30 -.5*z*z*(xlaplace)
				goto 12

c				indfi=3:2 dim. measured type table(bz(x,y,0))variable y increment
c					assume: bz(x,y,-z)=bz(x,y,z)
703			continue
				if(cyl(indre))stop' not yet impl. with cylindr. coord. option'

c				find the mesh:
				rim=(x-xmin)/xstep
				im=int(rim)
				ex=rim-real(im)
				eject=eject.or.im.lt.1
				i=im+1
				ip=i+1
				eject=eject.or.i.gt.(ndx-1)

c				use again the previous value of j?:
				j=jmem
				jm=j-1
				jp=j+1
				if(y.lt.ylis(j))then

c					search backward:
					do 3 j=jm,2,-1
						if(ylis(j).lt.y)then
							jp=j+1
							goto 4
						endif
3					continue
					eject=.true.
					bid=0.
				endif
				j1=jp+1
				if(y.gt.ylis(jp))then

c					search forward:
					do 6 jp=j1,ndy
						j=jp-1
						if(ylis(jp).gt.y)goto4
6					continue
					eject=.true.
				endif
4				if(eject)then
					jmem=2
					write(diagt,'(5(a,i3))')
     @			'exits field map :x index='
     @			,i,'[1,',ndx,'] y index=',j,'[1,',ndy,']'
					return
				endif
				jm=j-1
				yl=ylis(j)
				ylp=ylis(jp)
				ey=(y-yl)/(ylp-yl)
				py=ylp-yl
				ph=yl-ylis(jm)
				px=xstep
				zh=z

c				2dim. interpolation/extrapolation at 2nd order:
				c0=b(i,j,1,3)
				c1=b(im,j,1,3)
				c2=b(ip,j,1,3)
				c3=b(i,jm,1,3)
				c4=b(i,jp,1,3)
				c5=((c4*ph**2-c3*py**2)/(ph+py)-(ph-py)*c0)/ph
				c6=c2+c1-2.*c0
				c7=(c4*ph+c3*py)/(ph+py)-c0
				c8=b(ip,jp,1,3)+c0-c2-c4
				f(1)=(0.5*(c2-c1)+ex*c6+ey*c8)/px*zh
				f(2)=((c5+ex*c8)/py+ey*2.*c7/ph)*zh
				f(3)=c0+0.5*ex*(c2-c1)+ey*c5+0.5*(ex**2-(zh/px)**2)*c6+ex*ey*c8
     @			+(ey**2*py/ph-zh**2/(py*ph))*c7
				goto 12

c				indfi=4:2 dim. poisson type table(bz(0,y,z) and by(0,y,z))
c					assume: bz(x,y,z)=bz(0,y,z)
704			stop'field:not yet implemented'

c				indfi=5:3 dim. table
705			continue

c				find the cube:
				xr=(x-xmin)/xstep
				ixm=int(xr)
				if(logx)then
					if(mod(ixm,2).eq.0)ixm=ixm-1
					ix=ixm+2
					ixp=ix+2
					p=.5*(xr-real(ix-1))
				else
					ix=ixm+1
					ixp=ix+1
					p=xr-real(ixm)
				endif
				eject=eject.or.ixm.lt.1
				eject=eject.or.ixp.gt.ndx
				yr=(y-ymin)/ystep
				iym=int(yr)
				if(logy)then
					if(mod(iym,2).eq.0)iym=iym-1
					iy=iym+2
					iyp=iy+2
					q=.5*(yr-real(iy-1))
				else
					iy=iym+1
					iyp=iy+1
					q=yr-real(iym)
				endif   
				eject=eject.or.iym.lt.1
				eject=eject.or.iyp.gt.ndy
				zr=(z-zmin)/zstep
				izm=int(zr)
				if(logz)then
					if(mod(izm,2).eq.0)izm=izm-1
					iz=izm+2
					izp=iz+2
					r=.5*(zr-real(iz-1))
				else
					iz=izm+1
					izp=iz+1
					r=zr-real(izm)
				endif
				eject=eject.or.izm.lt.1
				eject=eject.or.izp.gt.ndz
				if(eject)then
					write(diagt,'(7(a,i3))')
     @			'exits field map :x index='
     @			,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     @			,'] z index=',iz,'[1,',ndz,']'
					return
				endif

c				test the geometry:
				if(eject)then
					return
				endif

c				3 dim. 2nd order interpolation inside the cube for each
c					component of the field.
c				use 11 values of the field:
				do 7 i=1,3
					f000=b(ix,iy,iz,i)
					eject=eject.or.(f000.gt.50.)
					f00p=b(ix,iy,izp,i)
					eject=eject.or.(f00p.gt.50.)
					f00m=b(ix,iy,izm,i)
					eject=eject.or.(f00m.gt.50.)
					f0p0=b(ix,iyp,iz,i)
					eject=eject.or.(f0p0.gt.50.)
					f0m0=b(ix,iym,iz,i)
					eject=eject.or.(f0m0.gt.50.)
					fp00=b(ixp,iy,iz,i)
					eject=eject.or.(fp00.gt.50.)
					fm00=b(ixm,iy,iz,i)
					eject=eject.or.(fm00.gt.50.)
					f0pp=b(ix,iyp,izp,i)
					eject=eject.or.(f0pp.gt.50.)
					fp0p=b(ixp,iy,izp,i)
					eject=eject.or.(fp0p.gt.50.)
					fpp0=b(ixp,iyp,iz,i)
					eject=eject.or.(fpp0.gt.50.)
					fppp=b(ixp,iyp,izp,i)
					eject=eject.or.(fppp.gt.50.)
					if(eject)then
						write(diagt,'(a,7(a,i3))')'try to use undefined '
     @				,'field values :x index='
     @				,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     @				,'] z index=',iz,'[1,',ndz,']'
						return
					endif
					c=2.*f000
					cp=fp00-c+fm00
					cq=f0p0-c+f0m0
					cr=f00p-c+f00m
					dp=f000-fp00+fpp0-f0p0
					dq=f000-f0p0+f0pp-f00p
					dr=f000-f00p+fp0p-fp00
					e=-dp-f0pp+f00p-fp0p+fppp

c					compute the taylor's serie:
					pq=p*q
					pqr=pq*r
					qr=q*r
					pr=p*r
					f(i)=f000+.5*(p*(fp00-fm00)+p*p*cp
     @				+q*(f0p0-f0m0)+q*q*cq
     @				+r*(f00p-f00m)+r*r*cr)
     @				+pq*dp+qr*dq+pr*dr+pqr*e
7				continue
				goto 12

c				indfi=6:3 dim. table, no extra mesh needed for interpolation
706			continue

c				find the cube:
				xr=(x-xmin)/xstep
				xrlim=real(ndx-1)*.5
				if(xr.ge.xrlim)then
					if(logx)then
						ix=int(xr/2.)*2+1
						ixm=ix-2
						ixp=ix+2
						p=.5*(xr-real(ix-1))
						eject=eject.or.ixm.lt.1
						eject=eject.or.ixp.gt.ndx
					else
						ixm=int(xr)
						ix=ixm+1
						ixp=ix+1
						p=xr-real(ixm)
						eject=eject.or.ixm.lt.1
						eject=eject.or.ixp.gt.ndx
					endif
				else
					if(logx)then
						ix=int(xr/2.)*2+3
						ixm=ix+2
						ixp=ix-2
						p=-.5*(xr-real(ix-1))
						eject=eject.or.ixm.gt.ndx
						eject=eject.or.ixp.lt.1
					else
						ix=int(xr)+2
						ixp=ix-1
						ixm=ix+1
						p=real(ix-1)-xr
						eject=eject.or.ixm.gt.ndx
						eject=eject.or.ixp.lt.1
					endif
				endif

				yr=(y-ymin)/ystep
				yrlim=real(ndy-1)*.5
				if(yr.ge.yrlim)then
					if(logy)then
						iy=int(yr/2.)*2+1
						iym=iy-2
						iyp=iy+2
						q=.5*(yr-real(iy-1))
						eject=eject.or.iym.lt.1
						eject=eject.or.iyp.gt.ndy
					else
						iym=int(yr)
						iy=iym+1
						iyp=iy+1
						q=yr-real(iym)
						eject=eject.or.iym.lt.1
						eject=eject.or.iyp.gt.ndy
					endif
				else
					if(logy)then
						iy=int(yr/2.)*2+3
						iym=iy+2
						iyp=iy-2
						q=-.5*(yr-real(iy-1))
						eject=eject.or.iym.gt.ndy
						eject=eject.or.iyp.lt.1
					else
						iy=int(yr)+2
						iyp=iy-1
						iym=iy+1
						q=real(iy-1)-yr
						eject=eject.or.iym.gt.ndy
						eject=eject.or.iyp.lt.1
					endif
				endif

				zr=(z-zmin)/zstep
				zrlim=real(ndz-1)*.5
				if(zr.ge.zrlim)then
					if(logz)then
						iz=int(zr/2.)*2+1
						izm=iz-2
						izp=iz+2
						r=.5*(zr-real(iz-1))
						eject=eject.or.izm.lt.1
						eject=eject.or.izp.gt.ndz
					else
						izm=int(zr)
						iz=izm+1
						izp=iz+1
						r=zr-real(izm)
						eject=eject.or.izm.lt.1
						eject=eject.or.izp.gt.ndz
					endif
				else
					if(logz)then
						iz=int(zr/2.)*2+3
						izm=iz+2
						izp=iz-2
						r=-.5*(zr-real(iz-1))
						eject=eject.or.izm.gt.ndz
						eject=eject.or.izp.lt.1
					else
						iz=int(zr)+2
						izp=iz-1
						izm=iz+1
						r=real(iz-1)-zr
						eject=eject.or.izm.gt.ndz
						eject=eject.or.izp.lt.1
					endif
				endif 
				if(eject)then
					write(diagt,'(7(a,i3))')
     @			'exits field map :x index='
     @			,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     @			,'] z index=',iz,'[1,',ndz,']'
					return
				endif

c				test the geometry:
				if(eject)then
					return
				endif

c				3 dim. 2nd order interpolation inside the cube for each
c					component of the field.
c				use 11 values of the field:
				do  i=1,3
					f000=b(ix,iy,iz,i)
					eject=eject.or.(abs(f000).gt.50.)
					f00p=b(ix,iy,izp,i)
					eject=eject.or.(abs(f00p).gt.50.)
					f00m=b(ix,iy,izm,i)
					eject=eject.or.(abs(f00m).gt.50.)
					f0p0=b(ix,iyp,iz,i)
					eject=eject.or.(abs(f0p0).gt.50.)
					f0m0=b(ix,iym,iz,i)
					eject=eject.or.(abs(f0m0).gt.50.)
					fp00=b(ixp,iy,iz,i)
					eject=eject.or.(abs(fp00).gt.50.)
					fm00=b(ixm,iy,iz,i)
					eject=eject.or.(abs(fm00).gt.50.)
					f0pp=b(ix,iyp,izp,i)
					eject=eject.or.(abs(f0pp).gt.50.)
					fp0p=b(ixp,iy,izp,i)
					eject=eject.or.(abs(fp0p).gt.50.)
					fpp0=b(ixp,iyp,iz,i)
					eject=eject.or.(abs(fpp0).gt.50.)
					fppp=b(ixp,iyp,izp,i)
					eject=eject.or.(abs(fppp).gt.50.)
					if(eject)then
						write(diagt,'(a,7(a,i3))')'try to use undefined '
     @				,'field values :x index='
     @				,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     @				,'] z index=',iz,'[1,',ndz,']'
						return
					endif
					c=2.*f000
					cp=fp00-c+fm00
					cq=f0p0-c+f0m0
					cr=f00p-c+f00m
					dp=f000-fp00+fpp0-f0p0
					dq=f000-f0p0+f0pp-f00p
					dr=f000-f00p+fp0p-fp00
					e=-dp-f0pp+f00p-fp0p+fppp

c					compute the taylor's serie:
					pq=p*q
					pqr=pq*r
					qr=q*r
					pr=p*r
					f(i)=f000+.5*(p*(fp00-fm00)+p*p*cp
     @				+q*(f0p0-f0m0)+q*q*cq
     @				+r*(f00p-f00m)+r*r*cr)
     @				+pq*dp+qr*dq+pr*dr+pqr*e
				enddo
				goto 12

c				3/range restitution due to symetry
12			if(x.ne.xyz(1))then
					if(sym(1).eq.'S'.or.sym(1).eq.'s')then
						f(2)=-f(2)
						f(3)=-f(3)
					else
						f(1)=-f(1)
					endif
				endif
				if(y.ne.xyz(2))then
					if(sym(2).eq.'S'.or.sym(2).eq.'s')then
						f(3)=-f(3)
						f(1)=-f(1)
					else
						f(2)=-f(2)
					endif
				endif
				if(z.ne.xyz(3))then
					if(sym(3).eq.'S'.or.sym(3).eq.'s')then
						f(1)=-f(1)
						f(2)=-f(2)
					else
						f(3)=-f(3)
					endif
				endif

c				4/normalize the field and restore cartesian coord.
c					and field components :
1000		continue
				if(eject)then
					return
				endif
				if(cyl(indre))then
					cst=cos(xyz(2))
					sst=sin(xyz(2))
					bx=f(1)*cst-f(2)*sst
					by=f(1)*sst+f(2)*cst
					f(1)=bx
					f(2)=by
					do 14 i=1,3
						xyz(i)=xyzc(i)
14				continue
				endif
				do 2 i=1,3
					f(i)=f(i)*fact(indre)
2				continue
				return
				end



				subroutine dist(a,b,c,d2,d)
c*************************************************************************
c                                                                        *
c                       d i s t                                          *
c                                                                        *
c*************************************************************************
c				compute d=distance between a point a and a straight line
c					of equation b + d*c (a,b and c = vectors , c assumed of unit length
				dimension a(3),b(3),c(3),bma(3)
				bmasc=0.
				do 29 i=1,3
					bmac=b(i)-a(i)
					bma(i)=bmac
					bmasc=bmasc+bmac*c(i)
29			continue
				d2=0.
				do 32 i=1,3
					dc=bma(i)-bmasc*c(i)
					d2=d2+dc*dc
32			continue
				d=sqrt(d2)
				return
				end
        


				subroutine mapadd(file,factor,outfile)        
c				subroutine that takes 2 3-d maps for SNAKE (output of mapwrt_3d)
c					and adds them together having multiplied each by a scaling
c					factor 
c				file(1)=name of 1st map file
c				factor(1)=multiplicative factor for 1st map
c				file(2)=name of 2nd map file
c				factor(2)=multiplicative factor for 2nd map
c				outfile=name of output map file
				Parameter (nn=900000)
				integer n(2),ndx(2),ndy(2),ndz(2),index
				character*1 typ(2,3)
				real factor(2)
				real xmax(2),ymax(2),zmax(2)
				real xmin(2),ymin(2),zmin(2),b(3,nn)
				character*80 file(2),outfile
				character*10 frmt(2)
				write(6,*)file(1),file(2)
				write(6,*)factor(1),factor(2)
				open(15,file=file(1),status='old',form='formatted')  !,err=97)
				open(16,file=file(2),status='old',form='formatted',err=98)
				do i=1,2
					k=i+14
					read(k,'(i2)')n(i)
					read(k,'(3a1)') (typ(i,j),j=1,3)
					read(k,*) ndx(i), xmin(i), xmax(i)
					read(k,*) ndy(i), ymin(i), ymax(i)
					read(k,*) ndz(i), zmin(i), zmax(i)
					read(k,*) frmt(i)
					write(6,'(i2)')n(i)
					write(6,'(3a1)') (typ(i,j),j=1,3)
					write(6,*) ndx(i), xmin(i), xmax(i)
					write(6,*) ndy(i), ymin(i), ymax(i)
					write(6,*) ndz(i), zmin(i), zmax(i)
					write(6,*) frmt(i)
				enddo

c				check that files match
				if((n(1).ne.n(2))
     @		.or.(ndx(1).ne.ndx(2))
     @		.or.(ndy(1).ne.ndy(2))
     @		.or.(ndz(1).ne.ndz(2))
     @		.or.(xmax(1).ne.xmax(2))
     @		.or.(ymax(1).ne.ymax(2))
     @		.or.(zmax(1).ne.zmax(2))
     @		.or.(xmin(1).ne.xmin(2))
     @		.or.(ymin(1).ne.ymin(2))
     @		.or.(zmin(1).ne.zmin(2))
     @		.or.(typ(1,1).ne.typ(2,1))
     @		.or.(typ(1,2).ne.typ(2,2))
     @		.or.(typ(1,3).ne.typ(2,3))) then
					go to 96
				else
					write(6,*) 'map files match'
				endif

c			open output file
				open(3,file=outfile,status='unknown',form='formatted',err=95)
				write(3,'(i2)')n(1)
				write(3,'(3a1)') (typ(1,j),j=1,3)
				write(3,*) ndx(1), xmin(1), xmax(1)
				write(3,*) ndy(1), ymin(1), ymax(1)
				write(3,*) ndz(1), zmin(1), zmax(1)
				write(3,*) frmt(1) 

				index=ndx(1)*ndy(1)*ndz(1)*3
				do i=1,2
					k=i+14
					read(k,(frmt(i)))(b(i,j),j=1,index)
				enddo
				close(15)
				close(16)
				do i=1,index
					b(3,i)=(b(1,i)*factor(1))+(b(2,i)*factor(2))
				enddo
				write(3,(frmt(1)))(b(3,i),i=1,index)
				close(3)    

				return 
95			stop 'error opening outfile'
96			stop 'map files do not match'
98			stop 'open error on file 2'   
				end 



        subroutine userout(in1,in2)
c*************************************************************************
c                                                                        *
c                       u s e r o u t                                    *
c                                                                        *
c*************************************************************************
        include 'sSnake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
        iw=0
        do 521 i=1,nuve
          alive=vea(i,in1,9)
          if(alive.le.0.)goto521
          alive=vea(i,in2,9)
          if(alive.le.0.)goto521
          iw=iw+1
          tl1=vea(i,in1,8)
          tl2=vea(i,in2,8)

c         built a set of independent variables:
          xr1=vea(i,in1,10)
          axr1=asin(vea(i,in1,13))
          zr1=vea(i,in1,12)
          azr1=asin(vea(i,in1,15))
          pxr1=(vea(i,in1,19))
          pyr1=(vea(i,in1,20))
          pzr1=(vea(i,in1,21))
          xr2=vea(i,in2,10)
          axr2=asin(vea(i,in2,13))
          zr2=vea(i,in2,12)
          azr2=asin(vea(i,in2,15))
          pxr2=(vea(i,in2,19))
          pyr2=(vea(i,in2,20))
          pzr2=(vea(i,in2,21))
          dl=tl2-tl1
          psq=vea(i,in1,7)

c         temporary: convert to transport's conventions (cm , mr proj. , %):
          psqref=vea(1,0,7)
          psq=((psq-psqref)/psqref)*100.
          cy1r=vea(1,in1,14)
          axr1r=atan2(vea(1,in1,13),cy1r)
          cy1=vea(i,in1,14)
          xr1=(xr1-vea(1,in1,10))/10.
          axr1=(atan2(vea(i,in1,13),cy1)-axr1r)*1000.
          zr1=(zr1-vea(1,in1,12))/10.
          azr1=atan2(vea(i,in1,15),cy1)*1000.
          cy2=vea(i,in2,14)
          xr2=(xr2-vea(1,in2,10))/10.
          cy2r=vea(1,in2,14)
          axr2r=atan2(vea(1,in2,13),cy2r)
          axr2=(atan2(vea(i,in2,13),cy2)-axr2r)*1000.
          zr2=(zr2-vea(1,in2,12))/10.
          azr2=atan2(vea(i,in2,15),cy2)*1000.
          tl1r=vea(1,in1,8)
          tl2r=vea(1,in2,8)
          dlr=tl2r-tl1r
          dl=-(dl-dlr)/10.

c         temporary end
523       write(4,'(10e14.6)')xr1,axr1,zr1,azr1,xr2,axr2,zr2,azr2,dl
     @          ,psq
521     continue
        return
        end



        subroutine spring
c*************************************************************************
c                                                                        *
c                       s p r i n g                                      *
c                                                                        *
c*************************************************************************
c
c                create the initial vector(s) set (end-plane index=0)
c     use 6 nested loops whose parameters (initial,maximum and increment value)
c         are given by the user
c       output:
c             nuve:number of input vector(s) defined
c             vea(iv,0,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord. in absolute frame
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum in absol.frame
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c            10 to 15:the same as 1 to 6 but in relative frame coordinates
c      16:x component of the spin
c      17:y          "
c      18:z          "
c      19 to 21: the same as 16 to 18 but in relative frame coord.
c             indre=current region   index
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(1,traj.index,x y or z):absol.coord.of submit
c
        include 'sSnake.inc'
        character rep*1,filtra*20,fitr*20,fibl*20,line*60
        logical lrand,ldisk,lref
        logical lbox,luser,bare
        integer ul
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @              ,iplpi(maxplv),iplpf(maxplv)
     @              ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @              ,luser,nclic,lclic(maxre),clight(4,maxre)
     @              ,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @              ,nrs(maxre)
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                ,frbox(2,maxre,3)
     @                ,zang(maxre),xang(maxre),yang(maxre)
     @                ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                ,indrer(maxre),tept(maxepr,maxre)
     @                ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                ,colli(maxepr,maxre,6)
				common ntraj
        data fitr/'t.dat'/
        data fibl/' '/
400     nuve=0
300     format(3f13.6/)
330     format(3f13.6)

c       data from terminal or disk ?
c        write(6,'(a)')'$ data from t(erminal) or d(isk) (def.=t):'
c        read(5, '(a)')rep
c        ldisk=rep.eq.'d'.or.rep.eq.'D'
				ldisk=.true.
        if(ldisk)then
c          write(6,'(''$ trajectory-file name(def='',a,'')'')')fitr
c          read(5,'(a)')filtra
c          if(filtra.ne.fibl)fitr=filtra
c          write(6,*)fitr
          ul=7  
          open(ul,file=fitr,status='unknown')
        else
          ul=5
        endif

c       loop random or individual ?
c        write(6,'(a)')'$ l(oops) r(andom) or i(ndividual) (def.=l):'
        read(ul, '(a)')rep
c        if(ldisk)write(6,*)rep
        if(rep.eq.' ')rep='l'
        lrand=rep.eq.'R'.or.rep.eq.'r'
        if(lrand)then
          write(6,'(a)')'$ number of trajectories :'
          read(ul,*)istat
          if(ldisk)write(6,*)istat
          open(8,file='sn-rand.dat',status='old')
          read(8,*)inidef
          close(8)
          write(6,'(a,i15,a)')
     @          '$ initial value of random generator (def from disk='
     @          ,inidef,'):'
          read(ul,'(i15)')iniran
          if(iniran.eq.0)iniran=inidef
          write(6,*)iniran
          call ranset(iniran)
        endif

c       case rep.e.'i' : loop or random :
        if(rep.ne.'I'.and.rep.ne.'i')then

c         input x min, max, and step          
          write(6,'(''$ x(def=0.,0.,1.mm):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)xmin,xmax,xstep
          else
            xmin=0.
            xmax=0.
            xstep=0.
          endif
          if(xstep.eq.0.)xstep=1.
          if(((xmax-xmin)*xstep).lt.0.)xmax=xmin
          write(6,330)xmin,xmax,xstep

c         input y min, max, and step          
          write(6,'(''$ y(def=0.,0.,1.mm):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)ymin,ymax,ystep
          else
            ymin=0.
            ymax=0.
            ystep=0.
          endif
          if(ystep.eq.0.)ystep=1.
          if(((ymax-ymin)*ystep).lt.0.)ymax=ymin
          write(6,330)ymin,ymax,ystep

c         input z min, max, and step
          write(6,'(''$ z(def=0.,0.,1.mm):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)zmin,zmax,zstep
          else
            zmin=0.
            zmax=0.
            zstep=0.
          endif
          if(zstep.eq.0.)zstep=1.
          if(((zmax-zmin)*zstep).lt.0.)zmax=zmin
          write(6,330)zmin,zmax,zstep

c         input theta x min, max, and step          
          write(6,'(''$ theta x (proj.angle)''
     @          ,''(def=0.,0.,1.rad):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)txmin,txmax,txstep
          else
            txmin=0.
            txmax=0.
            txstep=0.
          endif
          if(txstep.eq.0.)txstep=1.
          if(((txmax-txmin)*txstep).lt.0.)txmax=txmin
          write(6,330)txmin,txmax,txstep

c         input theta z min, max, and step          
          write(6,'(''$ theta z (space angle)''
     @          ,''(def=0.,0.,1.rad):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)tzmin,tzmax,tzstep
          else
            tzmin=0.
            tzmax=0.
            tzstep=0.
          endif
          if(tzstep.eq.0.)tzstep=1.
          if(((tzmax-tzmin)*tzstep).lt.0.)tzmax=tzmin
          write(6,330)tzmin,tzmax,tzstep

c         input p/q min, max, and step
          write(6,'(''$ p/q(def=1.,1.,1.gev/c):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)psqmin,psqmax,psqstep
          else
            psqmin=1.
            psqmax=1.
            psqstep=1.
          endif
          if(psqstep.eq.0.)psqstep=1.
          if(((psqmax-psqmin)*psqstep).lt.0.)psqmax=psqmin
          write(6,330)psqmin,psqmax,psqstep

c         input spin x min, max, and step
          write(6,'(''$ spin x''
     @          ,''(def=0.,0.,1.):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)spxmin,spxmax,spxstep
          else
            spxmin=0.
            spxmax=0.
            spxstep=0.
          endif
          if(spxstep.eq.0.)spxstep=1.
          if(((spxmax-spxmin)*spxstep).lt.0.)spxmax=spxmin
          write(6,330)spxmin,spxmax,spxstep

c         input spin y min, max, and step
          write(6,'(''$ spin y''
     @          ,''(def=0.,0.,1.):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)spymin,spymax,spystep
          else
            spymin=0.
            spymax=0.
            spystep=0.
          endif
          if(spystep.eq.0.)spystep=1.
          if(((spymax-spymin)*spystep).lt.0.)spymax=spymin
          write(6,330)spymin,spymax,spystep

c         input spin z min, ma, and step
          write(6,'(''$ spin z''
     @          ,''(def=0.,0.,1.):min,max,step?'')')
          read(ul,'(a)')line
          if(line(1:1).ne.' ')then
            read(line,*)spzmin,spzmax,spzstep
          else
            spzmin=0.
            spzmax=0.
            spzstep=0.
          endif
          if(spzstep.eq.0.)spzstep=1.
          if(((spzmax-spzmin)*spzstep).lt.0.)spzmax=spzmin
          write(6,330)spzmin,spzmax,spzstep

c         insert ref. traj. ?
          write(6,'(a)')'$ insert ref. traj. ? (y or n , def.=y)'
          read(ul,'(a)')rep
          if(ldisk)write(6,*)rep
          lref=rep.ne.'N'.and.rep.ne.'n'

c         ok ?
          write(6,'(a)')'$ o.k. ? (y or n , def.=y)'
          read(5,'(a)')rep
          if(rep.eq.'N'.or.rep.eq.'n')goto400

c         if random , then compute the loop param. 
c           to obtain istat*1*1*1*1*1=istat 
          if(lrand)then
            ixmax=istat
            iymax=1
            izmax=1
            itxmax=1
            itzmax=1
            ipsqmax=1
            ispxmax=1
            ispymax=1
            ispzmax=1

c         else compute the loop parameter :
          else
            ixmax=(xmax-xmin)/xstep+1.5
            iymax=(ymax-ymin)/ystep+1.5
            izmax=(zmax-zmin)/zstep+1.5
            itxmax=(txmax-txmin)/txstep+1.5
            itzmax=(tzmax-tzmin)/tzstep+1.5
            ipsqmax=(psqmax-psqmin)/psqstep+1.5
            ispxmax=(spxmax-spxmin)/spxstep+1.5
            ispymax=(spymax-spymin)/spystep+1.5
            ispzmax=(spzmax-spzmin)/spzstep+1.5
          endif

c         creation of the input vectors:
c           insert reference trajectory in 1st position ?
          if(lref)then
            call inject(0.,0.,0.,0.,1.,0.,(psqmin+psqmax)*.5
     @          ,0.,1.,0.)
          endif

c         loop on the 9 indep. components(x,y,z,tx,tz,p,px,py,pz):
c         x:
          do 1 ix=1,ixmax
            x=(ix-1)*xstep+xmin

c           y:
            do 2 iy=1,iymax
              y=(iy-1)*ystep+ymin

c             z:
              do 3 iz=1,izmax
                z=(iz-1)*zstep+zmin

c               theta x:
                do 4 itx=1,itxmax
                  tx=(itx-1)*txstep+txmin

c                 theta z:
                  do 5 itz=1,itzmax
                    tz=(itz-1)*tzstep+tzmin

c                   p:
                    do 6 ipsq=1,ipsqmax
                      psq=(ipsq-1)*psqstep+psqmin

c                     spin x:
                      do 7 ispx=1,ispxmax
                        px=(ispx-1)*spxstep+spxmin

c                       spin y:
                        do 8 ispy=1,ispymax
                          py=(ispy-1)*spystep+spymin

c                         spin z:
                          do 9 ispz=1,ispzmax
                            pz=(ispz-1)*spzstep+spzmin

                            if(lrand)then
                              x=aphasa(xmin,xmax)
                              y=aphasa(ymin,ymax)
                              z=aphasa(zmin,zmax)
                              tx=aphasa(txmin,txmax)
                              tz=aphasa(tzmin,tzmax)
                              psq=aphasa(psqmin,psqmax)
                              px=aphasa(spxmin,spxmax)
                              py=aphasa(spymin,spymax)
                              pz=aphasa(spzmin,spzmax)
                            endif
                            cz=sin(tz)
                            cx=sin(tx)*cos(tz)
                            cy2=1.-cx*cx-cz*cz
                            if(cy2.lt.0.)stop' spring:bad angles'
                            cy=sqrt(cy2)
                            if(nuve.ge.maxve)then
                              write(6,*)' room for only',maxve
     @                                  ,' vector(s)'
                              write(6,*)'  first ignored vector:'
                              write(6,100)x,y,z,tx,tz,psq,px,py,pz
100                           format(' x=',f10.3,' y=',f10.3,' z='
     @                              ,f10.3,' tx=',f10.6,' tz=',f10.6
     @                              ,' p=',f8.3,' sx=',f8.3,' sy='
     @                              ,f8.3,' sz=',f8.3)
                              goto10
                            endif
                            call inject(x,y,z,cx,cy,cz,psq,px,py,pz)
9                         continue
8                       continue
7                     continue
6                   continue
5                 continue
4               continue
3             continue
2           continue
1         continue
10        continue
          if(lrand)then

            call ranget(iniran)
            write(6,'(a,i15/a)')' value of the random generator :'
     @            ,iniran,' ( saved on disk file ''snake-rand.dat'')'
            open(8,file='sn-rand.dat',status='old')
            write(8,*)iniran
            close(8)
          endif

        else
c       case rep.eq.'i' : individual
c          write(6,'(a)')'$ number of trajectories :'
          read(ul,*)istat
					ntraj=istat
c          if(ldisk)write(6,*)istat
          do 108,ix=1,istat
            if(.not. ldisk)write(6,*)' type in x, y, z, tx, tz, psq'
     @                                ,'px, py, pz    for ray #',ix
            read(ul,*)x,y,z,tx,tz,psq
c            if(ldisk)write(6,100)x,y,z,tx,tz,psq
            cz=sin(tz)
            cx=sin(tx)*cos(tz)
            cy2=1.-cx*cx-cz*cz
            if(cy2.lt.0.)stop' spring:bad angles'
            cy=sqrt(cy2)
            if(nuve.ge.maxve)then
              write(6,*)' room for only',maxve,' vector(s)'
              write(6,*)'  first ignored vector:'
              write(6,100)x,y,z,tx,tz,psq,px,py,pz
              goto107
            endif
            call inject(x,y,z,cx,cy,cz,psq,px,py,pz)
108       continue
        endif

107     if(ldisk)close(ul)
        return
        end



        subroutine inject(x,y,z,cx,cy,cz,psq,px,py,pz)
c*************************************************************************
c                                                                        *
c                       i n j e c t                                      *
c                                                                        *
c*************************************************************************
        include 'sSnake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                ,frbox(2,maxre,3)
     @                ,zang(maxre),xang(maxre),yang(maxre)
     @                ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                ,indrer(maxre),tept(maxepr,maxre)
     @                ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                ,colli(maxepr,maxre,6)
        nuve=nuve+1
        indve=nuve
        ilive=nuve

c       fill the abs. coodr. buff. with this vector,region index=0
        indre=0
        vea(indve,0,1)=x
        vea(indve,0,2)=y
        vea(indve,0,3)=z
        vea(indve,0,4)=cx
        vea(indve,0,5)=cy
        vea(indve,0,6)=cz
        vea(indve,0,7)=psq
        vea(indve,0,8)=0.
        vea(indve,0,9)=1.
        vea(indve,0,10)=x
        vea(indve,0,11)=y
        vea(indve,0,12)=z
        vea(indve,0,13)=cx
        vea(indve,0,14)=cy
        vea(indve,0,15)=cz
        vea(indve,0,16)=px
        vea(indve,0,17)=py
        vea(indve,0,18)=pz
        vea(indve,0,19)=px
        vea(indve,0,20)=py
        vea(indve,0,21)=pz

c       fill the plot buffer with absolute data if room:
        call plotin(indve,x,y,z,0.,0.,0.,px,py,pz,.true.)
        return
        end



        subroutine send(kep,yep,tep,iep)
c*************************************************************************
c                                                                        *
c                       s e n d                                          *
c                                                                        *
c*************************************************************************
c                send the vectors from one end-plane to the next one
c       according to method.control the death of the vectors
c       input:
c             nuve:number of input vector(s) defined
c             method(region index):methode to use to carry particules
c                in this region
c             itype(region index):type in this method
c             fact(region index):multipl.factor of the field or ref.mommentum
c       input/output:
c             ver(iv,ic):description of the particules in the current region
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord.relative to the current region(in relative frame)
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum(in relat.frame)
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c      10:x component of the spin (in relat.frame)
c      11:y          "
c      12:z          "
c

        include 'sSnake.inc'
        parameter (maxmir=10)
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @              ,iplpi(maxplv),iplpf(maxplv)
     @              ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @              ,luser,nclic,lclic(maxre),clight(4,maxre)
     @              ,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @              ,nrs(maxre)
        common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     @        ,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
        common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
        character diag*100,diagt*100
        common/cdiag/diag(maxve),diagt
        common/ndiag/iepdead(maxve),iredead(maxve)
        character mdiag*3,qdiag(6)*1
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                ,frbox(2,maxre,3)
     @                ,zang(maxre),xang(maxre),yang(maxre)
     @                ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                ,indrer(maxre),tept(maxepr,maxre)
     @                ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                ,colli(maxepr,maxre,6)

        logical cyl,lfil
        character*80 fifi
        common /fidef/method(maxre),itype(maxre),indfi(maxre)
     @                ,fact(maxre),cyl(maxre),lfil
        common /cfidef/ fifi(maxre)
        logical eject,lcer
        common/alter/jflag,prec,hstmax,hstmin,sinc
        common/params/pin(10),fin,ctep,step,idir,epsil,pout(10)
     @                ,h,ifail,pal
        common/matrix/a(6,6),dy,s(3,3)
        dimension q(6)
        equivalence (x,q(1)),(y,q(2)),(z,q(3)),(xp,q(4)),(yp,q(5))
     @              ,(zp,q(6))
        dimension xyz(3),pm(5,3)
        data qdiag/'x','y','z','r','t','z'/


c       data for rungk:

c       end-plane orthogonal to axis y (1:x,2:y,3:z):
        idir=kep

c       accuracy of the final value of the coord. idir:(mm)
        epsil=.001

c       variable step (jflag=1):
        jflag=1

c       step change:(adapted for spectro.900:r=1800.,half gap=60.)
c         threshold value of relative field change for step change:
        prec=.02
        ctep=cos(tep)
        step=sin(tep)
        fin=yep*ctep

c       cerenkov mirrors?:
        lcer=cerfiln(iep,indre).ne.'none'
        if(.not.lcer)goto 197
        open(9,file=cerfiln(iep,indre),status='old')
        read(9,*)nmir,rmirin,rmirout,dzmir,dxmir

c       next temp:
        if(nmir.gt.maxmir)then
          nmir=maxmir
          write(6,*)'too many mirrors,max=',maxmir
        endif
        read(9,*)axp,axm,azp,azm

c       next temp:
        write(6,*)axp,axm,azp,azm
        axp=axp*dtor
        axm=axm*dtor
        azp=azp*dtor
        azm=azm*dtor
        do 198 imir=1,nmir
          read(9,*)(qmir(iq,imir),iq=1,3),tiltz(imir),tiltx(imir)

c         next temp:
          tiltz(imir)=tiltz(imir)*dtor
          tiltx(imir)=tiltx(imir)*dtor

c         I call "middle" of a mirror the barycenter of its 4 corners
c           distance (qmir="middle" of the mirror)
c           _(cv=center of the sphere):
          xm=sqrt(rmirin**2-.25*(dxmir**2+dzmir**2))

c         sv and tv are two unit vectors, perpendicular between them. 
c           s is normal to the planes limiting the mirror in z'',
c           t is normal to the planes limiting the mirrors in x''
c           .(x'',z'') coincidate with (xrel.,zrel.) if tiltx=tiltz=0. 
c           To go from rel. frame to '' frame one rotates first /z, 
c           second /x
          cosz=cos(tiltz(imir))
          sinz=sin(tiltz(imir))
          cosx=cos(tiltx(imir))
          sinx=sin(tiltx(imir))
          tv(1,imir)=cosz
          tv(2,imir)=sinz*cosx
          tv(3,imir)=sinz*sinx
          sv(1,imir)=0.
          sv(2,imir)=-sinx
          sv(3,imir)=cosx

c         unit vector between "middle" and center:
          cv(1,imir)=tv(2,imir)*sv(3,imir)-tv(3,imir)*sv(2,imir)
          cv(2,imir)=tv(3,imir)*sv(1,imir)-tv(1,imir)*sv(3,imir)
          cv(3,imir)=tv(1,imir)*sv(2,imir)-tv(2,imir)*sv(1,imir)

          do 932 iq=1,3

c           compute the position ( cv ) of the center of this mirror:
            cv(iq,imir)=qmir(iq,imir)+cv(iq,imir)*xm

c           plot this mirror as a rectangle whose submits are those of 
c             the inner (reflecting) surface of the mirror:
            pm(1,iq)=qmir(iq,imir)+.5*(tv(iq,imir)*dxmir+sv(iq,imir)
     @				*dzmir)
            pm(2,iq)=qmir(iq,imir)+.5*(tv(iq,imir)*dxmir-sv(iq,imir)
     @				*dzmir)
            pm(3,iq)=qmir(iq,imir)+.5*(-tv(iq,imir)*dxmir-sv(iq,imir)
     @				*dzmir)
            pm(4,iq)=qmir(iq,imir)+.5*(-tv(iq,imir)*dxmir+sv(iq,imir)
     @				*dzmir)
            pm(5,iq)=pm(1,iq)
932       continue  
          nulu=nulu+1
          if(nulu.gt.nulumax)then
            nulu=nulumax
            write(6,*)' warning send: too many user line'
            goto198
          endif
          luser=.true.
          nuseu(nulu)=5
          if(nuseu(nulu).gt.nuseumax)nuseu(nulu)=nuseumax
          do 400 is=1,nuseu(nulu)
400         call plorot(pm(is,1),pm(is,2),pm(is,3)
     @          ,plotus(is,nulu,1),plotus(is,nulu,2),plotus(is,nulu,3)
     @          ,xb,yb,zb)
198     continue
        close(9)
197     continue

c       loop over the nuve vectors:
        do 1 indve=1,nuve
          i=indve

c         check if trajectory is dead
          if(ver(i,9).le.0.)then
            do 2 j=1,8
              ver(i,j)=0.
2           continue
            do 3 j=10,12
              ver(i,j)=0.
3           continue

c           if it was dead or buried,now it is buried:
            ver(i,9)=-1.
            goto1
          endif

          diag(indve)=' '
          x=ver(i,1)
          y=ver(i,2)
          z=ver(i,3)
          cx=ver(i,4)
          cy=ver(i,5)
          cz=ver(i,6)
          psq=ver(i,7)
          px=ver(i,10)
          py=ver(i,11)
          pz=ver(i,12)

c         input point in the free-box?
          if(cyl(indre))then
            xyz(1)=sqrt(x*x+y*y)
            xyz(2)=atan2(y,x)     
          else
            xyz(1)=x
            xyz(2)=y
          endif
          xyz(3)=z
          eject=.false.
          do 11 k=1,3
            if(xyz(k).lt.frbox(1,indre,k))then
              idiag=k
              mdiag='min'
              eject=.true.
            endif
11        continue
          if(cyl(indre))idiag=idiag+3

c         built the diag:
          if(eject)then
            write(diag(indve),'(3a)')'outside the free box in :'
     @            ,qdiag(idiag),mdiag
            goto 200
          endif

c         cerenkov? If yes, emitt light from the current point with the
c           current direction; the radiator is the previous end-plane:
          if(lcer)then

c           x,y,z: point of emitted cerenkov photon
c           cx,cy,cz: direction of emitted photon
c           x1,y1,z1: point of reflection of photon from mirror "imir"
c           x2,y2,z2: point of interception of photon from "imir1"
            call mirror(x,y,z,cx,cy,cz,x1,y1,z1,x2,y2,z2,imir,imir1)
            call plorot(x,y,z,xa,ya,za,xb,yb,zb)
            call plotin(indve,xa,ya,za,0.,0.,0.,0.,0.,0.,.false.)
            call plorot(x1,y1,z1,xa1,ya1,za1,xb,yb,zb)
            call plotin(indve,xa1,ya1,za1,0.,0.,0.,0.,0.,0.,.false.)
            call plorot(x2,y2,z2,xa2,ya2,za2,xb,yb,zb)
            call plotin(indve,xa2,ya2,za2,0.,0.,0.,0.,0.,0.,.false.)
            call plotin(indve,xa1,ya1,za1,0.,0.,0.,0.,0.,0.,.false.)
            call plotin(indve,xa,ya,za,0.,0.,0.,0.,0.,0.,.false.)
          endif

c         field-less,matrix or rung-khuta type?
          goto(101,102,103),method(indre)

c         no field:
101       continue
          yp=-x*step+y*ctep
          cyp=-cx*step+cy*ctep
          if(cyp.le.0.)goto15
          gly=(fin-yp)/cyp
          if(gly.lt.0.)goto15
          xt=x+gly*cx
          yt=y+gly*cy
          zt=z+gly*cz

c         output point in the free box?
          if(xt.gt.frbox(1,indre,1).and.xt.lt.frbox(2,indre,1).and.
     @      yt.gt.frbox(1,indre,2).and.yt.lt.frbox(2,indre,2).and.
     @      zt.gt.frbox(1,indre,3).and.zt.lt.frbox(2,indre,3))then
            
            x=xt
            y=yt
            z=zt
            pal=gly

c           test for collimator clearance
            if(colli(iep,indre,1).ne.0.)then
              call collitest(x,z,iep,i,eject)
            endif
            goto199
          endif

c         find the escape point:
          glxz=1.e30
          if(xt.le.frbox(1,indre,1))then
            glxt=(frbox(1,indre,1)-x)/cx
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=1
              mdiag='min'
            endif
          elseif(xt.ge.frbox(2,indre,1))then
            glxt=(frbox(2,indre,1)-x)/cx
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=1
              mdiag='max'
            endif
          endif
          if(yt.le.frbox(1,indre,2))then
            glxt=(frbox(1,indre,2)-y)/cy
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=2
              mdiag='min'
            endif
          elseif(yt.ge.frbox(2,indre,2))then
            glxt=(frbox(2,indre,2)-y)/cy
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=2
              mdiag='max'
            endif
          endif
          if(zt.le.frbox(1,indre,3))then
            glxt=(frbox(1,indre,3)-z)/cz
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=3
              mdiag='min'
            endif
          elseif(zt.ge.frbox(2,indre,3))then
            glxt=(frbox(2,indre,3)-z)/cz
            if(glxt.lt.glxz)then
              glxz=glxt
              idiag=3
              mdiag='max'
            endif
          endif
          x=x+glxz*cx
          y=y+glxz*cy
          z=z+glxz*cz
          pal=glxz
          eject=.true.

c         build the diag:
          write(diag(indve),'(4a)')'exits free box in :'
     @          ,qdiag(idiag),mdiag,' in a no field region'
          goto 199

c         can't reach the end-plane:
15        pal=0.
          eject=.true.

c         build the diag:
          write(diag(indve),'(2a)')'cant reach the end plane'
     @          , ' in a no field region'
          goto 200

c         cerenkov mirrors?:
199       continue
          if(cerfiln(iep,indre).eq.'none')goto 200
          goto 200

c         matrix:
102       continue
          tx=asin(cx)
          tz=asin(cz)
          if(fact(indre).eq.0.)stop'send/1stmatrix:ref. momentum=0.'
          del=psq/fact(indre)-1.

c         length of the reference trajectory in the directive file
c           (first auxiliary data) :
          pal=adata(1,indre)
          goto(301,302),itype(indre)

c         1st order matrix:
301       xt=x*a(1,1)+tx*a(2,1)+z*a(3,1)+tz*a(4,1)+del*a(6,1)
          txt=x*a(1,2)+tx*a(2,2)+z*a(3,2)+tz*a(4,2)+del*a(6,2)
          zt=x*a(1,3)+tx*a(2,3)+z*a(3,3)+tz*a(4,3)+del*a(6,3)
          tzt=x*a(1,4)+tx*a(2,4)+z*a(3,4)+tz*a(4,4)+del*a(6,4)
          pl=x*a(1,5)+tx*a(2,5)+z*a(3,5)+tz*a(4,5)+del*a(6,5)
          pal=pal+pl
          cx=sin(txt)
          cz=sin(tzt)
          x=xt
          z=zt
          y=yep
          pxt=px*s(1,1)+py*s(2,1)+pz*s(3,1)
          pyt=px*s(1,2)+py*s(2,2)+pz*s(3,2)
          pzt=px*s(1,3)+py*s(2,3)+pz*s(3,3)
          px=pxt
          py=pyt
          pz=pzt
          cy2=1.-cx*cx-cz*cz
          if(cy2.lt.0.)then
            eject=.true.

c           build the diag:
            write(diag(indve),'(2a)')'crazy direction after 1st order'
     @            , ' matrix'
            goto 200
          endif
          cy=sqrt(cy2)

c         output point in the free-box
          do 12 k=1,3
            if(q(k).lt.frbox(1,indre,k))then
              idiag=k
              mdiag='min'
              eject=.true.
            endif
            if(q(k).gt.frbox(2,indre,k))then
              idiag=k
              mdiag='max'
              eject=.true.
            endif
12        continue

c         build the diag:
          if(eject)then
            write(diag(indve),'(4a)')'exit free box in :'
     s            ,qdiag(idiag),mdiag,' after 1st order matrix'
            goto 200
          endif

c         test for collimator clearance
          if(colli(iep,indre,1).ne.0.) call collitest(x,z,iep,i,eject)
          goto 200

c         2nd order matrix
302       stop'send/2nd matrix:not yet implemented'

c         rung-khuta:
c           fill the standard rungk buffer:
103       continue
          pin(1)=x
          pin(2)=y
          pin(3)=z
          pin(4)=cx
          pin(5)=cy
          pin(6)=cz
          pin(7)=px
          pin(8)=py
          pin(9)=pz
          pin(10)=psq
          h=hstmin
          call rungk(eject,i)
          x=pout(1)
          y=pout(2)
          z=pout(3)
          cx=pout(4)
          cy=pout(5)
          cz=pout(6)
          px=pout(7)
          py=pout(8)
          pz=pout(9)
          psq=pout(10)
          if(eject)goto 200

c         test for collimator clearance
          if(colli(iep,indre,1).ne.0.) call collitest(x,z,iep,i,eject)
200       continue
          ver(i,1)=x
          ver(i,2)=y
          ver(i,3)=z
          ver(i,4)=cx
          ver(i,5)=cy
          ver(i,6)=cz
          ver(i,7)=psq
          ver(i,8)=ver(i,8)+pal
          if(eject)then
            ver(i,9)=0.
            iepdead(i)=indep+1
            iredead(i)=indre
            diag(i)=diagt
          endif
          ver(i,10)=px
          ver(i,11)=py
          ver(i,12)=pz

1       continue
        return
        end



        subroutine collitest(x,z,iep,iv,eject)
c*************************************************************************
c                                                                        *
c                       c o l l i t e s t                                *
c                                                                        *
c*************************************************************************
        include 'sSnake.inc'
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                ,frbox(2,maxre,3)
     @                ,zang(maxre),xang(maxre),yang(maxre)
     @                ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                ,indrer(maxre),tept(maxepr,maxre)
     @                ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                ,colli(maxepr,maxre,6)
        character diag*100,diagt*100
        common/cdiag/diag(maxve),diagt
        common/ndiag/iepdead(maxve),iredead(maxve)
        logical eject
        real k,m1,m2

c       there are three kinds of collimators; trapazoid, ellipse
c         , and blocker; parameters in colli are as follows
c                           trapezoid      ellipse   blocker
c       colli(iep,indre,1)    xmin          xmax      xmin
c       colli(iep,indre,2)    xmax          zmax      xmax
c       colli(iep,indre,3)    zmax1          h        zmin
c       colli(iep,indre,4)    zmax2          k        zmax
c       colli(iep,indre,5)    zmin1          -         -
c       colli(iep,indre,6)    zmin2          -         999.

        if(colli(iep,indre,1).eq.0.) go to 1
        if(colli(iep,indre,2).eq.0.) go to 1

c       test for trapezoid or ellipse
        if((colli(iep,indre,5).eq.0.)
     @    .and.(colli(iep,indre,6).eq.0.))then

c         ellipse
          xmax=colli(iep,indre,1)
          zmax=colli(iep,indre,2)
          h=colli(iep,indre,3)
          k=colli(iep,indre,4)
          ans=((x-h)/xmax)**2+((z-k)/zmax)**2
          if(ans.gt.1.)then
            eject=.true.

c           build the diag:
            write(diag(iv),'(2a,i3,a)')'discarded by collimator'
     @      ,' at end-plane #',iep,' of current region'
          endif
          go to 1

c       trapezoid
        elseif(colli(iep,indre,6).ne.999.
     @		.and.colli(iep,indre,6).ne.888.
     @		.and.colli(iep,indre,6).ne.777.)then
          xmin=colli(iep,indre,1)
          xmax=colli(iep,indre,2)
          zmax1=colli(iep,indre,3)
          zmax2=colli(iep,indre,4)
          zmin1=colli(iep,indre,5)
          zmin2=colli(iep,indre,6)
          if(xmax.eq.xmin)then
            write(6,*)' ************colli error xmax=xmin**********'
            go to 1
          endif
          m1=(zmax2-zmax1)/(xmax-xmin)
          b1=zmax1-(m1*xmin)
          m2=(zmin2-zmin1)/(xmax-xmin)
          b2=zmin1-(m2*xmin)
          if((z.lt.(x*m2)+b2).or.(z.gt.(x*m1)+b1)
     @      .or.(x.lt.xmin).or.(x.gt.xmax))then
            eject=.true.

c           build the diag:
            write(diag(iv),'(2a,i3,a)')'discarded by collimator'
     @            ,' at end-plane #',iep,' of current region'
c            write(6,110)iv,iep
          endif

c				Dual Collimator
				elseif(colli(iep,indre,6).ne.999.
     @		.and.colli(iep,indre,6).ne.777.)then
					xlim=colli(iep,indre,1)
					zmin1=colli(iep,indre,2)
					zmax1=colli(iep,indre,3)
					zmin2=colli(iep,indre,4)
					zmax2=colli(iep,indre,5)
					if((x.gt.xlim).or.(x.lt.-xlim))then
						eject=.true.
c						write(6,110)iv,iep
					elseif((z.lt.zmin1).or.(z.gt.zmax2))then
						eject=.true.
c						write(6,110)iv,iep
					elseif((z.gt.zmax1).and.(z.lt.zmin2))then
						eject=.true.
c						write(6,110)iv,iep
					endif
						
c				Dual Collimator 2
				elseif(colli(iep,indre,6).ne.999.)then
					xmin=colli(iep,indre,1)
					xmax=colli(iep,indre,2)
					zlim1=colli(iep,indre,3)
					zlim2=colli(iep,indre,4)
					if((x.gt.xmax).or.(x.lt.xmin))then
						eject=.true.
					elseif((z.lt.-zlim2).or.(z.gt.zlim2))then
						eject=.true.
					elseif((z.gt.-zlim1).and.(z.lt.zlim1))then
						eject=.true.
					endif

c       blocking rectangular collimator
        else
          xmin=colli(iep,indre,1)
          xmax=colli(iep,indre,2)
          zmin=colli(iep,indre,3)
          zmax=colli(iep,indre,4)
          if(xmax.eq.xmin)then
            write(6,*)' ************colli error xmax=xmin**********'
            go to 1
          endif
          if((z.lt.zmax).and.(z.gt.zmin)
     @      .and.(x.lt.xmax).and.(x.gt.xmin))then
            eject=.true.
            write(diag(iv),'(2a,i3,a)')'discarded by blocker'
     @      ,' at end-plane #',iep,' of current region'
c           write(6,110)iv,iep
          endif
        endif
110     format(1x,'vector #',i4,'   intercepted by collimator #',i3)
1       return
        end



        subroutine mirror(x,y,z,cx,cy,cz,x1,y1,z1,x2,y2,z2,imir,imir1)
c*************************************************************************
c                                                                        *
c                       m i r r o r                                        *
c                                                                        *
c*************************************************************************
        parameter (maxmir=10)
        common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     @    ,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
        common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
        dimension av(3),uv(3),bv(3),vv(3),wv(3)
        dimension bvold(3),vvold(3),fvold(3)
        av(1)=x
        av(2)=y
        av(3)=z
        uv(1)=cx
        uv(2)=cy
        uv(3)=cz
        xlold=1.e30
        iold=0
        imir=0
        imir1=0
        rmir=rmirin
        imirn=0
        sign=1.
        call sphere(av,uv,xlold,iold,imirn,rmir,bvold,vvold,sign)
        if(iold.gt.0)then
          imir=iold
          do 10 i=1,3
            bv(i)=bvold(i)
10          vv(i)=vvold(i)
          x1=bv(1)
          y1=bv(2)
          z1=bv(3)
        else

c         intersection of incident particle:
          x1=x
          y1=y
          z1=z
          x2=x
          y2=y
          z2=z
          return
        endif

c       calculation of particle reflection:
        vsu=0.
        do 4 i=1,3
4         vsu=vsu+vv(i)*uv(i)
        do 5 i=1,3
5         wv(i)=uv(i)-2.*vsu*vv(i)

c       particle reflection from another mirror
        xlold=.5*rmirin
        iold=0
        rmir=rmirin
        imirn=imir
        sign=-1.
        call sphere(bv,wv,xlold,iold,imirn,rmir,fvold,vvold,sign)
        if(iold.gt.0)then
          imir1=iold
          x2=fvold(1)
          y2=fvold(2)
          z2=fvold(3)
        else
          x2=bv(1)+.5*rmirin*wv(1)
          y2=bv(2)+.5*rmirin*wv(2)
          z2=bv(3)+.5*rmirin*wv(3)
        endif
        return
        end



        subroutine sphere(av,uv,xlold,iold,imirn,rmir,bvold,vvold,sign)
c*************************************************************************
c                                                                        *
c                       s p h e r e                                      *
c                                                                        *
c*************************************************************************
c
c  computes the intersection between a ray and a set of "spherical rectangles"
c  called "mirrors".
c  If there is an intersection, and if it occurs after a flight distance
c  smaller than "xlold", updates "xlold", returns the intersection point "bv",
c  the unit vector "vv" (from "bv" to the center "cv" of the concerned mirror),
c  the number "iold" of this mirror. Dont test mirror#"imirn". Dont keep
c  solutions for which xl<0. If "sign"=+1. choose the solution for which 
c  xl is the greatest (snallest for "sign"=-1). Returns 0 for "iold" if
c  no intersection. If sevral intersection, keep that for which xl is 
c  the smallest.
        parameter (maxmir=10)
        common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     s    ,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
        common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
        dimension av(3),uv(3),dg(3),bv(3),dv(3),vv(3),ts(3),vn(3),gv(3)
        dimension bvold(3),vvold(3)

        do 2 jmir=1,nmir
          if(jmir.eq.imirn)goto2

c         unit vector between "middle" and center:
          ts(1)=tv(2,jmir)*sv(3,jmir)-tv(3,jmir)*sv(2,jmir)
          ts(2)=tv(3,jmir)*sv(1,jmir)-tv(1,jmir)*sv(3,jmir)
          ts(3)=tv(1,jmir)*sv(2,jmir)-tv(2,jmir)*sv(1,jmir)

          xl0=0.
          dg2=0.
          do 1 i=1,3
            dg(i)=cv(i,jmir)-av(i)
            xl0=xl0+dg(i)*uv(i)
1           dg2=dg2+dg(i)*dg(i)
          delta=xl0*xl0-dg2+rmir*rmir
          if(delta.le.0.)goto2

          xl=xl0+sign*sqrt(delta)
          if(xl.lt.0.)goto 2
          if(xl.gt.xlold)goto 2

c         limits of the mirror
          det=0.
          des=0.
          do 3 i=1,3
            bv(i)=av(i)+xl*uv(i)
            dv(i)=bv(i)-cv(i,jmir)
3           vv(i)=dv(i)/rmir
          proj=0.
          do 4 i=1,3
            vn(i)=tv(i,jmir)*cos(axp)-ts(i)*sin(axp)
            gv(i)=qmir(i,jmir)+.5*tv(i,jmir)*dxmir-bv(i)
4           proj=proj+vn(i)*gv(i)
          if(proj.lt.0.)goto 2
          proj=0.
          do 5 i=1,3
            vn(i)=tv(i,jmir)*cos(axm)-ts(i)*sin(axm)
            gv(i)=qmir(i,jmir)-.5*tv(i,jmir)*dxmir-bv(i)
5           proj=proj+vn(i)*gv(i)
          if(proj.gt.0.)goto2
          proj=0.
          do 6 i=1,3
            vn(i)=sv(i,jmir)*cos(azp)+ts(i)*sin(azp)
            gv(i)=qmir(i,jmir)+.5*sv(i,jmir)*dzmir-bv(i)
6           proj=proj+vn(i)*gv(i)
          if(proj.lt.0.)goto2
          proj=0.
          do 7 i=1,3
            vn(i)=sv(i,jmir)*cos(azm)+ts(i)*sin(azm)
            gv(i)=qmir(i,jmir)-.5*sv(i,jmir)*dzmir-bv(i)
7           proj=proj+vn(i)*gv(i)
          if(proj.gt.0.)goto2
          do 9 i=1,3
            bvold(i)=bv(i)
9           vvold(i)=vv(i)
          xlold=xl
          iold=jmir
2       continue
        return
        end


        subroutine rungk(eject,indve)
c*************************************************************************
c                                                                        *
c                       r u n g k                                        *
c                                                                        *
c*************************************************************************
c
c             carry a charged particule through the magnetic field
c                given by the subroutine "field" by integrating
c                the equation of motion by the rung-khuta method.
c                the integration stops on a "end-plane".
c
c             input parameters:
c                pin:input vector(x,y,z,cx,cy,cz,px,py,pz,p)
c         (mm,gev/c)
c                jflag=1:autom.adjust. of the step when field change
c                     =2:step fixed to "hstmin"(mm)
c                prec:threshold value of rel. field change for step adj.
c                hstmax,hstmin:max and min value of step
c                sinc:abs.val. of the relative step change
c                fin:coord. of the end-plane
c                idir=1:the end-plane is orthogonal to x axis
c                    =2: "   "    "    "     "       " y  "
c                    =3: "   "    "    "     "       " z  "
c                epsil:acuracy of the end-plane reachig
c                h:current value of the step length(mm)
c             output parameters:
c                pout:output vector(x,y,z,cx,cy,cz,px,py,pz,p)
c         (mm,gev/c)
c             output flag:
c                ifail=0:o.k.
c                     =1:can't reach the end-plane
c                eject:the subr. field has been called outside its
c                  definition region, or the partic. can't reach the end-plane
c
        include 'sSnake.inc'
        parameter (tmax=.017)
        parameter (turn=2*pi)
        parameter (xme=.510999e-3,re=2.817940e-12)
        parameter (coef=(2./3.)*re*xme)
        logical eject
        logical lpplot
        common/alter/jflag,prec,hstmax,hstmin,sinc
        logical lsynchro,lspin
        common/logic/lsynchro,lspin
        common/particle/z_particle,xm
        common/spin/x_mag_mon,e_over_m,twice_spin,xmref
        character diag*100,diagt*100
        common/cdiag/diag(maxve),diagt
        common/ndiag/iepdead(maxve),iredead(maxve)
        logical lbox,luser,bare
        logical logx,logy,logz
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @              ,iplpi(maxplv),iplpf(maxplv)
     @              ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @              ,luser,nclic,lclic(maxre),clight(4,maxre)
     @              ,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @              ,nrs(maxre)
        common/params/pin(10),fin,ctep,step,idir,epsil,pout(10)
     @                ,h,ifail,pal
        dimension f(3)
        dimension secxs(4),secys(4),seczs(4),cold(3),csync(3)
        dimension spinxs(4),spinys(4),spinzs(4)
        dimension hsave(3),hnext(3),hmean(3),hgrad(3)
        dimension xyz(3),xyzt(3),abc(3),xz(9),xx(9),pxyz(3)
        equivalence (x,xyz(1)),(y,xyz(2)),(z,xyz(3)),(xt,xyzt(1))
     @              ,(yt,xyzt(2)),(zt,xyzt(3)),(a,abc(1)),(b,abc(2))
     @              ,(c,abc(3)),(xz(1),xyz(1)),(xz(4),abc(1))
     @              ,(px,pxyz(1)),(py,pxyz(2)),(pz,pxyz(3))
     @              ,(xz(7),pxyz(1))
        logical start
        parameter (thrd=1./3.)
        data xlight_velocity/.299792458e-3/
        data logx/.false./,logy/.false./,logz/.false./
        xm2=xm*xm
        dp=0.
        call ucopy(pin ,xx,9)

c       convention: zero means infinite momentum --> no curvature:
c         zero mom. rays are plotted at each step
        p7=pin(10)
        if(p7.eq.0.)then
          pinv=0.
          gama=1.
        else
          pinv=xlight_velocity/p7
          gama=sqrt((p7*p7+xm2)/xm2)
        endif

c       computation of the coefficients for the spin precession:
        if (lspin) then
          gsur2=x_mag_mon*xm/(twice_spin*xmref)
          beta1 = 1. + gama * ( gsur2 - 1.)
          alpha1 = (gsur2 - beta1)
        endif
        ctmax=cos(tmax)
        tarn=0.
        pal=0.

c       tilted end-plane: the tilt axis is idir+1 . so the e-p stands
c         between the plane xx(idir)=0. (for tilt angle=0.)
c         and the plane xx(jdir)=0. (for tilt angle=90.deg)
c         where jdir=idir-1
        jdir=mod(idir+1,3)+1
        vtil=-xx(jdir)*step+xx(idir)*ctep
        do=vtil-fin
        sdo=sign(1.,do)
        iflag=1
        call ucopy(pin(4),csync,3)
        call ucopy(pin(4),cold,3)
        call field(eject,xx,f,logx,logy,logz,.true.)
        if(eject)then
          call ucopy(pin,xz,9)
          goto300
        endif
        if(jflag.eq.1) call ucopy(f,hsave,3)
        start=.true.
10      vtil=-xx(jdir)*step+xx(idir)*ctep
        vpil=-xx(jdir+3)*step+xx(idir+3)*ctep

c       kill the track if injected after the e.p:
        if(start.and.vtil.gt.fin)then
          eject=.true.
          call ucopy(pin,pout,10)
          return
        endif

c       plot the first point (to check field continuity):
        if(start)then
          call plorot(xx(1),xx(2),xx(3),xa,ya,za,xb,yb,zb)
          call plorot(f(1),f(2),f(3),xc,yc,zc,bx,by,bz)
          call plorot(px,py,pz,xc,yc,zc,pxa,pya,pza)
          call plotin(indve,xa,ya,za,bx,by,bz,pxa,pya,pza,.false.)
        endif
        start=.false.
        s=h*vpil+vtil
        d2=s-fin

c       don't test the end if y<y end or backw. direct.:
        if(d2.lt.0..or.vpil.lt.0.)goto 30

20      h=h-d2/vpil
        iflag=iflag+1

30      continue
        call ucopy(xx,xz,9)
        h2=0.5*h
        h4=0.25*h
        ph=pinv*h
        ph2=0.5*ph
        secxs(1)=(b*f(3)-c*f(2))*ph2*z_particle
        secys(1)=(c*f(1)-a*f(3))*ph2*z_particle
        seczs(1)=(a*f(2)-b*f(1))*ph2*z_particle

c       spin deb.
        if(lspin)then
          alpha2 = alpha1 * (a*f(1) + b*f(2) + c*f(3))
          spinxs(1) = (py*(alpha2*c + beta1*f(3))
     @              -  pz*(alpha2*b + beta1*f(2))) * ph2
          spinys(1) = (pz*(alpha2*a + beta1*f(1))
     @              -  px*(alpha2*c + beta1*f(3))) * ph2
          spinzs(1) = (px*(alpha2*b + beta1*f(2))
     @              -  py*(alpha2*a + beta1*f(1))) * ph2
          pxt=px+spinxs(1)
          pyt=py+spinys(1)
          pzt=pz+spinzs(1)
        endif

c       spin fin
        xt=x+h2*a+h4*secxs(1)
        yt=y+h2*b+h4*secys(1)
        zt=z+h2*c+h4*seczs(1)
        at=a+secxs(1)
        bt=b+secys(1)
        ct=c+seczs(1)
        if(iflag.le.2)call field(eject,xyzt,f,logx,logy,logz,.true.)
        if(eject)goto300
        secxs(2)=(bt*f(3)-ct*f(2))*ph2*z_particle
        secys(2)=(ct*f(1)-at*f(3))*ph2*z_particle
        seczs(2)=(at*f(2)-bt*f(1))*ph2*z_particle

c       spin deb.
        if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(2) = (pyt*(alpha2*ct + beta1*f(3))
     @              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(2) = (pzt*(alpha2*at + beta1*f(1))
     @              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(2) = (pxt*(alpha2*bt + beta1*f(2))
     @              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          pxt=px+spinxs(2)
          pyt=py+spinys(2)
          pzt=pz+spinzs(2)
        endif

c       spin fin
        at=a+secxs(2)
        bt=b+secys(2)
        ct=c+seczs(2)
        secxs(3)=(bt*f(3)-ct*f(2))*ph2*z_particle
        secys(3)=(ct*f(1)-at*f(3))*ph2*z_particle
        seczs(3)=(at*f(2)-bt*f(1))*ph2*z_particle

c       spin deb.
        if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(3) = (pyt*(alpha2*ct + beta1*f(3))
     @              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(3) = (pzt*(alpha2*at + beta1*f(1))
     @              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(3) = (pxt*(alpha2*bt + beta1*f(2))
     @              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          pxt=px+spinxs(3)+spinxs(3)
          pyt=py+spinys(3)+spinys(3)
          pzt=pz+spinzs(3)+spinzs(3)
        endif

c       spin fin
        xt=x+h*a+h*secxs(3)
        yt=y+h*b+h*secys(3)
        zt=z+h*c+h*seczs(3)
        at=a+secxs(3)+secxs(3)
        bt=b+secys(3)+secys(3)
        ct=c+seczs(3)+seczs(3)
        if(jflag.le.2) call field(eject,xyzt,f,logx,logy,logz,.true.)
        if(eject)goto300
        if(iflag.eq.1) call ucopy(f,hnext,3)
        secxs(4)=(bt*f(3)-ct*f(2))*ph2*z_particle
        secys(4)=(ct*f(1)-at*f(3))*ph2*z_particle
        seczs(4)=(at*f(2)-bt*f(1))*ph2*z_particle

c       spin deb.
        if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(4) = (pyt*(alpha2*ct + beta1*f(3))
     @              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(4) = (pzt*(alpha2*at + beta1*f(1))
     @              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(4) = (pxt*(alpha2*bt + beta1*f(2))
     @              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          px=px+(spinxs(1)+spinxs(4)+2.0*(spinxs(2)+spinxs(3)))*thrd
          py=py+(spinys(1)+spinys(4)+2.0*(spinys(2)+spinys(3)))*thrd
          pz=pz+(spinzs(1)+spinzs(4)+2.0*(spinzs(2)+spinzs(3)))*thrd
        endif

c       spin fin
        z=z+(c+(seczs(1)+seczs(2)+seczs(3))*thrd)*h
        y=y+(b+(secys(1)+secys(2)+secys(3))*thrd)*h
        x=x+(a+(secxs(1)+secxs(2)+secxs(3))*thrd)*h
        a=a+(secxs(1)+secxs(4)+2.0*(secxs(2)+secxs(3)))*thrd
        b=b+(secys(1)+secys(4)+2.0*(secys(2)+secys(3)))*thrd
        c=c+(seczs(1)+seczs(4)+2.0*(seczs(2)+seczs(3)))*thrd
        pal=pal+h

c       plot if there is enough change in direction:
        ctot=0.
        lpplot=.false.
        du2=0.
        do 1 i=1,3
          du=abc(i)-csync(i)
          du2=du2+du*du
1         ctot=ctot+cold(i)*abc(i)

c       sync.rad.aver.ener.loss:
        if(lsynchro.and.p7.ne.0.)then
          p2=p7*p7
          p3=p2*p7
          dp=coef*du2*p3*gama/(h*xm*xm2)
          p7=p7-dp
          pinv=xlight_velocity/p7
          gama=sqrt((p7*p7+xm2)/xm2)
          call ucopy(abc,csync,3)
        endif

c       check for trajectories trapped in the field
        if(ctot.lt.ctmax.or.p7.eq.0.)then
          tarn=tarn+acos(ctot)
          if(tarn.ge.turn)then
            eject=.true.

c           build the diag:
            write(diag(indve),'(3(a,e10.3))')
     @            'total rotation angle exceeds max.:'
     @            ,tarn, ' >',turn,' rad.'
            goto300
          endif
          call ucopy(abc,cold,3)
          lpplot=.true.
          call plorot(x,y,z,xa,ya,za,xb,yb,zb)
          call plorot(f(1),f(2),f(3),xc,yc,zc,bx,by,bz)
          call plorot(px,py,pz,xc,yc,zc,pxa,pya,pza)
          call plotin(indve,xa,ya,za,bx,by,bz,pxa,pya,pza,.false.)
        endif

200     format(8f10.4)
        if(iflag.gt.1)goto90
        if(jflag.eq.1)goto40
        call ucopy(xz,xx,9)
        goto 10
40      do 50 i=1,3
          hmean(i)=(hnext(i)+hsave(i))*0.5
50        hgrad(i)=hnext(i)-hsave(i)
        ht=vmodul(hmean)
        hratio=0
        sh=sign(1.,h)

c       if the field is smaller than 10**-2 tesla (100 gauss)
c         the step can increase or stay constant,else it can
c         also decrease
        if(ht.le.1.e-2) goto 60
        hg=vmodul(hgrad)
        hratio=hg/ht
        sh=sign(1.,h)
        if(abs(h).le.hstmin) goto 60
        if(hratio.gt.prec) goto80
60      call ucopy(xz,xx,9)
        if(hratio.ge.0.1*prec) goto 70
        h=h+sinc*sh
        h=amin1(abs(h),hstmax)*sh
70      call ucopy(hnext,hsave,3)
        goto 10
80      pal=pal-h
        if(p7.ne.0.)p7=p7+dp
        if(lpplot.and.(indve+nplvold).le.maxplv)then
          if(nplp(indve+nplvold).lt.maxplp)then
            nplp(indve+nplvold)=nplp(indve+nplvold)-1
          endif
        endif
        h=h-sinc*sh
        h=amax1(abs(h),hstmin)*sh
        call ucopy(hsave,f,3)
        goto 10
90      vtil=-xz(jdir)*step+xz(idir)*ctep
        d1=vtil-fin
        hold=abs(h)
        h=-d1/vpil
        if(abs(h).lt.epsil)goto100
        if(abs(h).gt.hold)goto120
        call ucopy(xz,xx,9)
        iflag=iflag+1
        goto 30
100     ifail=0
        do 110 mm=1,3
          pout(mm)=xyz(mm)+abc(mm)*h
          pout(mm+3)=abc(mm)
110       pout(mm+6)=pxyz(mm)
        pout(10)=p7
        pal=pal+h
        if(indve.le.maxplv)then
          call plorot(f(1),f(2),f(3),xa,ya,za,bx,by,bz)
          bplot(1,indve)=bx
          bplot(2,indve)=by
          bplot(3,indve)=bz
        endif
        goto 130
120     ifail=2

c       build the diag:
        write(diag(indve),'(a)')
     @        'unsuccessfull landing on the end plane'
        eject=.true.
        goto 300

130     return
300     continue
        do 301 mm=1,9
301       pout(mm)=xz(mm)
        pout(10)=p7
        return
        end

        function vmodul(a)
        dimension a(3)
        vmodul=sqrt(a(1)**2+a(2)**2+a(3)**2)
        return
        end



        subroutine ucopy(a,b,n)
        dimension a(*),b(*)
        if(n.eq.0) return
        do 2 ii=1,n
2         b(ii)=a(ii)
        return
        end


        subroutine circle(xc,yc,zc,r,nuc)
c*************************************************************************
c                                                                        *
c                       c i r c l e                                      *
c                                                                        *
c*************************************************************************
        include 'sSnake.inc'
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @              ,iplpi(maxplv),iplpf(maxplv)
     @              ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @              ,luser,nclic,lclic(maxre),clight(4,maxre)
     @              ,epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3)
     @              ,nrs(maxre)
        if(nulu.eq.nulumax)then
          write(6,*)' warning circle: too many user line'
          return
        endif
        nulu=nulu+1
        dt=2.*pi/real(nuseumax-1)
        t=0.
        do 1 i=1,nuseumax
          t=t+dt
          x=xc+r*cos(t)
          y=yc+r*sin(t)
          z=zc
          call plorot(x,y,z,xa,ya,za,xb,yb,zb)
          plotus(i,nuc,1)=xa
          plotus(i,nuc,2)=ya
1         plotus(i,nuc,3)=za
        nuseu(nuc)=nuseumax
        return
        end



        subroutine analyf(eject,xyz,f,index,indve)
c*************************************************************************
c                                                                        *
c                       a n a l y f                                      *
c                                                                        *
c*************************************************************************
c
c                find the 3 components of the field f(3) at the point xyz(3)
c         analyf is called by field in the case of analytic field definition
c         .true.output value for eject causes the death of the particule.
c
c       input:
c             xyz(x y or z):relative coord. of the point
c             index:field index of the region
c       output:
c             f(x y or z):field component along x y or z
c             eject:value .true. cause the death of the particule
c       for method=3, type=2
c         index=1 -> dipole with grad.
c         index=16/17 -> arditrary field numerical derivative
c         index=2-14 -> raytrace library
c         index=15 -> clam
        include 'sSnake.inc'
        character diag*100,diagt*100
        common/cdiag/diag(maxve),diagt
        common/ndiag/iepdead(maxve),iredead(maxve)
        common/aject/ieject
        common/dipgra/rm,alp,bet
        common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @              ,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @              ,tbound,dst
        common/raytra/fact1(5),qrad1,fact2(5),qrad2,dl
        common/clam1/d0,dx,dy,bdref,ta
        dimension f(3),xyz(3)
        logical eject
        real x,y,z
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
        goto(101,102,103,104,105,106,107,108,109
     @      ,110,111,112,113,114,115,116,117,118,119,120),index

c       index=1:field with gradient (alpha and beta)
c         (x=0.,y=0.,z=0.)=center of symmetry, z=0.:symmetry plane
c         rectangular-->cylindrical coord.:
101     r=sqrt(x*x+y*y)
        rn=r/rm
        te=atan2(y,x)
        zn=z/rm
        dr=rn-1.

c       z component of the field in sym.pl.and its derivatives/rn:
        b0=1.-alp*dr+bet*dr*dr
        db=-alp+2.*bet*dr
        d2b=+2.*bet

c       extrapolation from the symmetry plane:
        cr=zn*db
        f(3)=b0-.5*zn*zn*(d2b+db/rn)

c       cylindrical--->rectangular field components:
        f(1)=cr*cos(te)
        f(2)=cr*sin(te)
        return

c       index=2:unit gradient (1.t/mm ) along x and z (qpole axis=y )
102     f(1)=z
        f(2)=0.
        f(3)=x
        return

c       index=3:quadrupole
103     continue
        call quad(xyz(1),xyz(3),fact1(1),qrad1,f(1),f(3))
        f(2)=0.
        return

c       index=4:dipole entrance fringe
104     ydum=-xyz(2)
        xdum=xyz(1)
        zdum=xyz(3)
        call dfringe(xdum,ydum,zdum,fact1(1),qrad1,f(1),f(2),f(3))
        f(2)=-f(2)
        return

c       index=5:dipole exit fringe
105     ydum=xyz(2)
        xdum=xyz(1)
        zdum=xyz(3)
        call dfringe(xdum,ydum,zdum,fact1(1),qrad1,f(1),f(2),f(3))
        return

c       index=6:uniform dipole
106     f(1)=fact1(1)
        f(2)=fact1(2)
        f(3)=fact1(3)
        return

c       index=7:entrance multipole
107     xdum=xyz(1)
        ydum=-xyz(2)
        zdum=xyz(3)
        call mpoles(1,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        f(2)=-f(2)
        return

c       index=8:uniform field multipole
108     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call mpoles(2,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        return

c       index=9:exit multipole
109     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call mpoles(3,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        return

c       index=10 qq overlapping fringe fields
110     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call qqfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2
     @                ,qrad1,qrad2)
        return

c       index=11 qd overlapping fringe fields
111     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call qdfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2
     @                ,qrad1,qrad2)
        return

c       index=12 dd overlapping fringe fields
112     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call ddfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2
     @                ,qrad1,qrad2)
        return

c       index=13 dq overlapping fringe fields (old version,PV 4/9/96)
113     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call dqfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2
     @                ,qrad1,qrad2)
        return

c       index=14 dq overlapping fringe fields (new version,PV 4/9/96)
114     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call dqfringe2(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2
     @                ,qrad1,qrad2)
        return

c       index=15:central field of a clam:
115     d1=abs(d0+x*dx+y*dy)
        d2=d1*d1
        eject=d1.le.0..or.abs(z).gt.(d1*ta)
        if(eject)then

c         built the diag:
          write(diagt,'(a)')
     @      'try to go in a forbidden region of a clam'
          return
        endif
        ccf=bdref/(d2+z*z)
        f(1)=-ccf*z*dx
        f(2)=-ccf*z*dy
        f(3)=ccf*d1
        return

c       index=16:dipole with gradient numerical deriv:
116     call dipole(x,y,b00)
        call dipole(x+dst,y,bp0)
        call dipole(x-dst,y,bm0)
        call dipole(x,y+dst,b0p)
        call dipole(x,y-dst,b0m)
        div1=1./(2.*dst)
        div2=1./(dst*dst)
        grx1=(bp0-bm0)*div1
        gry1=(b0p-b0m)*div1
        b002=2.*b00
        grx2=bp0+bm0-b002
        gry2=b0p+b0m-b002
        f(1)=z*grx1
        f(2)=z*gry1
        f(3)=b00-.5*z*z*(grx2+gry2)*div2
        return

c       index=17:raytrace style dipole,old version:
117     xdum=xyz(1)/10.
        ydum=xyz(2)/10.
        zdum=xyz(3)/10.
        call dipolert(xdum,zdum,ydum,f(1),f(3),f(2))

c       get the right sign on the field
        do 12 i=1,3
12        f(i)=-f(i)
        return

c       index=18:raytrace style dipole, new version
118     xdum=xyz(1)/10.
        ydum=xyz(2)/10.
        zdum=xyz(3)/10.
        call dipolert2(xdum,zdum,ydum,f(1),f(3),f(2))

c       get the right sign on the field
        do 13 i=1,3
13        f(i)=-f(i)
        return

c       index=19:polarized target
119     ydum=xyz(2)/1000.
        xdum=xyz(1)/1000.
        zdum=xyz(3)/1000.
        call hlmhltz(xdum,ydum,zdum,qrad1,fact1(1),f(1),f(2),f(3))
        return

c       index=20: Solenoid
120     ydum=xyz(2)/1000.
        xdum=xyz(1)/1000.
        zdum=xyz(3)/1000.
        call solenoid(xdum,ydum,zdum,fact1(1),fact1(2),fact1(3)
     @                ,bx,by,bz)
        f(1)=bx
        f(2)=by
        f(3)=bz
        return
        end



        subroutine mpoles(in,nx,ny,nz,ngrad,nrad,nbx,nby,nbz)
c       stolen from raytrace jjl 10/8/86                                                                
c       calculation of multipole(poles) field components                  
c         2 - quadrupole  (grad1)                                           
c         3 - hexapole    (grad2)                                           
c         4 - octapole    (grad3)                                           
c         5 - decapole    (grad4)                                           
c         6 - dodecapole  (grad5)                                           
        implicit real*8(a-h,o-z)
        real nx,ny,nz,ngrad(5),nrad,nbx,nby,nbz
        common  /blck91/  c0, c1, c2, c3, c4, c5
        data c0, c1, c2, c3, c4, c5/0.1039d0,6.27108d0,-1.51247d0
     @                              ,3.59946d0,-2.1323d0,1.7230d0/
        if (nrad.eq.0.)then
          write(6,*)' error in mpoles,  nrad= 0.'
          call exit(0)
        endif
        grad1 = -ngrad(1)/nrad
        grad2 =  ngrad(2)/nrad**2
        grad3 = -ngrad(3)/nrad**3
        grad4 =  ngrad(4)/nrad**4
        grad5 = -ngrad(5)/nrad**5
        rad=nrad
        d = 2. * rad
        frh  = 1.d0
        fro  = 1.d0
        frd  = 1.d0
        frdd = 1.d0
        dh  = frh *d
        do  = fro *d
        dd  = frd *d
        ddd = frdd*d
        x = nx
        y = ny
        z = nz
        x2 = x*x
        x3 = x2*x
        x4 = x3*x
        x5 = x4*x
        x6 = x5*x
        x7 = x6*x
        y2 = y*y
        y3 = y2*y
        y4 = y3*y
        y5 = y4*y
        y6 = y5*y
        y7 = y6*y
        go to ( 2, 1, 2 ) , in
        print 3, in
3       format( '  error in bpoles in= ', i5 ///)
        call exit(0)
1       continue
        b2x = grad1*y
        b2y = grad1*x
        b3x = grad2*2.*x*y
        b3y = grad2*(x2-y2)
        b4x = grad3*(3.*x2*y-y3)
        b4y = grad3*(x3-3.*x*y2)
        b5x = grad4*4.*(x3*y-x*y3)
        b5y = grad4*(x4-6.*x2*y2+y4)
        b6x = grad5*(5.*x4*y-10.*x2*y3+y5)
        b6y = grad5*(x5-10.*x3*y2+5.*x*y4)
        bx = b2x + b3x + b4x + b5x + b6x
        by = b2y + b3y + b4y + b5y + b6y
        bz = 0.
        bt =   sqrt( bx*bx + by*by )
        nbx=bx
        nby=by
        nbz=bz
        return
2       s=z/d
        call bpls( 2, d, s, re, g1, g2, g3, g4, g5, g6 )
        b2x = grad1*( re*y - (g2/12.)*(3.*x2*y + y3)
     @      + (g4/384.)*(5.*x4*y + 6.*x2*y3 + y5 )
     @      - (g6/23040.)*(7.*x6*y + 15.*x4*y3 + 9.*x2*y5 + y7)  )
        b2y = grad1*( re*x - (g2/12.)*(x3 + 3.*x*y2)
     @      +(g4/384.)*(x5 + 6.*x3*y2 + 5.*x*y4 )
     @      - (g6/23040.)*(x7 + 9.*x5*y2 + 15.*x3*y4 + 7.*x*y6) )
        b2z = grad1*( g1*x*y - (g3/12.)*(x3*y + x*y3 )
     @      + (g5/384.)*(x5*y +2.*x3*y3 + x*y5)  )
        ss = z/dh  + dsh
        call bpls( 3, dh, ss, re, g1, g2, g3, g4, g5, g6 )
        b3x = grad2*( re*2.*x*y - (g2/48.)*(12.*x3*y + 4.*x*y3 ) )
        b3y = grad2*(re*(x2-y2) - (g2/48.)*(3.*x4 + 6.*x2*y2 - 5.*y4))
        b3z = grad2*(g1*(x2*y - y3/3.)-(g3/48.)*(3.*x4*y+2.*x2*y3-y5))
        ss = z/do  + dso
        call bpls( 4, do, ss, re, g1, g2, g3, g4, g5, g6 )
        b4x = grad3*( re*(3.*x2*y - y3) - (g4/80.)*(20.*x4*y - 4.*y5))
        b4y = grad3*( re*(x3 - 3.*x*y2) - (g4/80.)*(4.*x5-20.*x*y4 ) )
        b4z = grad3*g1*(x3*y - x*y3 )
        ss = z/dd  + dsd
        call bpls( 5, dd, ss, re, g1, g2, g3, g4, g5, g6 )
        b5x = grad4*re*(4.*x3*y - 4.*x*y3)
        b5y = grad4*re*(x4 - 6.*x2*y2 + y4 )
        b5z = grad4*g1*(x4*y - 2.*x2*y3 + y5/5.)
        ss = z/ddd + dsdd
        call bpls( 6, ddd,ss, re, g1, g2, g3, g4, g5, g6 )
        b6x = grad5*re*(5.*x4*y - 10.*x2*y3 + y5 )
        b6y = grad5*re*(x5 - 10.*x3*y2 + 5.*x*y4 )
        b6z = 0.
        bx = b2x + b3x + b4x + b5x + b6x
        by = b2y + b3y + b4y + b5y + b6y
        bz = b2z + b3z + b4z + b5z + b6z
        bt =   sqrt( bx*bx + by*by + bz*bz )
        nbx=bx
        nby=by
        nbz=bz
        return
        end



        subroutine bpls ( igp, d, s, re, g1, g2, g3, g4, g5, g6 )
        implicit real*8 (a-h,o-z)
        common  /blck91/  c0, c1, c2, c3, c4, c5
        s2 = s*s
        s3 = s2*s
        s4 = s2*s2
        s5 = s4*s
        cs = c0 + c1*s + c2*s2 + c3*s3 + c4*s4 + c5*s5
        cp1 =(c1 + 2.*c2*s + 3.*c3*s2 + 4.*c4*s3 + 5.*c5*s4) / d
        cp2 = (2.*c2 + 6.*c3*s + 12.*c4*s2 + 20.*c5*s3  ) / (d*d)
        cp3 = ( 6.*c3 + 24.*c4*s + 60.*c5*s2 ) / (d**3)
        cp4 = ( 24.*c4 + 120.*c5*s ) / (d**4)
        cp5 = 120.*c5/(d**5)
        if( abs(cs) .gt. 70. )  cs = sign(70.d0, cs )
        e = exp(cs)
        re = 1./(1. + e)
        ere = e*re
        ere1= ere*re
        ere2= ere*ere1
        ere3= ere*ere2
        ere4= ere*ere3
        ere5= ere*ere4
        ere6= ere*ere5
        cp12 = cp1*cp1
        cp13 = cp1*cp12
        cp14 = cp12*cp12
        cp22 = cp2*cp2
        cp15 = cp12*cp13
        cp16 = cp13*cp13
        cp23 = cp2*cp22
        cp32 = cp3*cp3
        if( igp .eq. 6 ) return
        g1 = -cp1*ere1
        if( igp .eq. 5 ) return
        if( igp .eq. 4 ) go to 1
        g2 =-( cp2+cp12   )*ere1    + 2.*cp12 * ere2
        g3 =-(cp3 + 3.*cp1*cp2 + cp13  ) * ere1
     @      +6.*(cp1*cp2 + cp13)*ere2 - 6.*cp13*ere3
        if( igp .eq. 3 ) return
1       g4 = -(cp4 + 4.*cp1*cp3 + 3.*cp22 + 6.*cp12*cp2 + cp14)*ere1
     @      + (8.*cp1*cp3 + 36.*cp12*cp2 + 6.*cp22 + 14.*cp14)*ere2
     @      - 36.*(cp12*cp2 + cp14)*ere3       + 24.*cp14*ere4
        if( igp .ne. 2 ) return
        g5 = (-cp5 - 5.*cp1*cp14 - 10.*cp2*cp3 - 10.*cp12*cp3
     @      - 15.*cp1*cp22 - 10.*cp13*cp2 - cp15)*ere1
     @      + (10.*cp1*cp4 +20.*cp2*cp3 +60.*cp12*cp3 + 90.*cp1*cp22
     @      + 140.*cp13*cp2 +30.*cp15)*ere2 + (-60.*cp12*cp3
     @      - 90.*cp1*cp22 - 360.*cp13*cp2 - 150.*cp15)*ere3
     @      + (240.*cp13*cp2 +240.*cp15)*ere4 + (-120.*cp15)*ere5
        g6 = (-6.*cp1*cp5 - 15.*cp2*cp4 - 15.*cp12*cp4 - 10.*cp32
     @      - 60.*cp1*cp2*cp3 - 20.*cp13*cp3 - 15.*cp23 - 45.*cp12*cp22
     @      - 15.*cp14*cp2 - cp16)*ere1 + (12.*cp1*cp5 + 30.*cp2*cp4
     @      + 90.*cp12*cp4 +20.*cp32 + 360.*cp1*cp2*cp3 +280.*cp13*cp3
     @      + 90.*cp23 + 630.*cp12*cp22 + 450.*cp14*cp2 + 62.*cp16)*ere2
     @      + (-90.*cp12*cp4 - 360.*cp1*cp2*cp3 -720.*cp13*cp3 -90.*cp23
     @      - 1620.*cp12*cp22 -2250.*cp14*cp2 - 540.*cp16)*ere3
     @      + (480.*cp13*cp3 + 1080.*cp12*cp22 + 3600.*cp14*cp2
     @      + 1560.*cp16)*ere4 + (-1800.*cp14*cp2 - 1800.*cp16)*ere5
     @      + 720.*cp16*ere6
        return
        end



        subroutine quad(x,y,fact,a,fx,fy)
c       calculates quadrupole field components fx and fy for given
c         input coordinates x,y     jjl 10/2/86
c
c       fact = max field strength at r=a
c       a = aperture radius 
        real x,y,fx,fy,fact,a
        if(x**2+y**2.gt.a**2)go to 900
        fx=-y*fact/a
        fy=-x*fact/a
        return
900     write(6,*)'point outside quadrupole aperture'
        return
        end


        subroutine dipolert (x,y,z,bbx,bby,bbz)
c*************************************************************************
c                                                                        *
c                       d i p o l e r t                                  *
c                                                                        *
c*************************************************************************
c
c  dipolert- analytic calculation of dipole magnetic fields
c            stolen from raytrace adapted for use in snake
c                                         -jjl  8/21/89
c    data(i) are read in snake in ifield
c    dipolert is called from snake in analyf
c                snake provides x,y,z in the c-axis system
c                dipolert returns bx,by,bz in the same system
        implicit real*8(a-h,o-z)
        real x,y,z,bbx,bby,bbz,data(75),pii
        real*8  ndx, k
        common  /diprt/data
        common  /block10old/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(6), dtc(6)
        data pii/0.0174532925/
        lf1  = dble(data(  1  ))
        lu1  = dble(data(  2  ))
        lf2  = dble(data(  3  ))
        dg   = dble(data(  4  ))
        mtyp = dble(data(  5  ))
        a    = dble(data( 11  ))
        b    = dble(data( 12  ))
        d    = dble(data( 13  ))
        rb   = dble(data( 14  ))
        bf   = dble(data( 15  ))
        phi  = dble(data( 16  ))
        alpha= dble(data( 17  ))
        beta = dble(data( 18  ))
        ndx  = dble(data( 19  ))
        bet1 = dble(data( 20  ))
        gama = dble(data( 21  ))
        delt = dble(data( 22  ))
        z11  = dble(data( 25  ))
        z12  = dble(data( 26  ))
        z21  = dble(data( 27  ))
        z22  = dble(data( 28  ))
        br1  = dble(data( 41  ))
        br2  = dble(data( 42  ))
        xcr1 = dble(data( 43  ))
        xcr2 = dble(data( 44  ))
        xpric=x+(2.*rb*sin(pii*phi/2.)*sin(pii*((phi/2.)-beta)))
        zpric=z+(2.*rb*sin(pii*phi/2.)*cos(pii*((phi/2.)-beta)))
        cosa=cos(pii*(phi-alpha-beta))
        sina=sin(pii*(phi-alpha-beta))
        tc(1)=(-xpric*cosa)+(zpric*sina)
        tc(2)=dble(y)
        tc(3)=(-xpric*sina)-(zpric*cosa)

c       test for region
        if(tc(3).gt.z11)then
          bx=0.
          by=0.
          bz=0.
          go to 10
        endif
        if(tc(3).ge.z12) go to 1
        if(tc(3).lt.z12) go to 2

c       in designates magnet regions for bfun
1       in = 1
        xc= rb*cos( alpha/ 57.29578 )
        zc=-rb*sin( alpha/ 57.29578 )
        c0   = dble(data( 29  ))
        c1   = dble(data( 30  ))
        c2   = dble(data( 31  ))
        c3   = dble(data( 32  ))
        c4   = dble(data( 33  ))
        c5   = dble(data( 34  ))
        dels = dble(data( 45  ))
        rca  = dble(data( 47  ))
        csc = cos( alpha/57.29578 )
        scor = dble(data(49 ))
        s2   = dble(data( 51  )) / rb    + rca/2.d0
        s3   = dble(data( 52  )) / rb**2
        s4   = dble(data( 53  )) / rb**3 + rca**3/8.d0
        s5   = dble(data( 54  )) / rb**4
        s6   = dble(data( 55  )) / rb**5 + rca**5/16.d0
        s7   = dble(data( 56  )) / rb**6
        s8   = dble(data( 57  )) / rb**7 + rca**7/25.6d0
        call ndipold

c       transform to c-axis system
        bxdum=bx
        bzdum=bz
        bx=(-bxdum*cosa)-(bzdum*sina)
        bz=(+bxdum*sina)-(bzdum*cosa)
        go to 10

c       transform to second vfb coord system
2       copab =cos( (phi-alpha-beta)/57.29578)
        sipab =sin( (phi-alpha-beta)/57.29578)
        cospb =cos( (phi/2.-beta)/57.29578 )
        sinpb =sin( (phi/2.-beta)/57.29578 )
        sip2 =sin( (phi/2.)/57.29578 )
        xt = tc(1)
        zt = tc(3)
        vxt = tc(4)
        vzt = tc(6)
        tc(3) = - zt  *copab +  xt  *sipab -2.*rb*sip2*cospb
        tc(1) = - zt  *sipab -  xt  *copab -2.*rb*sip2*sinpb
        tc(6) = - vzt *copab +  vxt *sipab
        tc(4) = - vzt *sipab -  vxt *copab
        if(tc(3).lt.z21)go to 3
        if(tc(3).gt.z22)then
          bx=0.
          by=0.
          bz=0.
          go to 10
        endif
        if(tc(3).le.z22)go to 4

c       uniform field integration region
3       in = 2
        xc=-rb*cos( beta / 57.29578 )
        zc=-rb*sin( beta / 57.29578 )
        call ndipold
        go to 10

c       setup for second fringe field and integration
4       br   = br2
        c0   = dble(data( 35  ))
        c1   = dble(data( 36  ))
        c2   = dble(data( 37  ))
        c3   = dble(data( 38  ))
        c4   = dble(data( 39  ))
        c5   = dble(data( 40  ))
        dels = dble(data( 46  ))
        rca  = dble(data( 48  ))
        scor = dble(data(50 ))
        csc = cos( beta /57.29578 )
        s2   = dble(data( 58  )) / rb    + rca/2.d0
        s3   = dble(data( 59  )) / rb**2
        s4   = dble(data( 60  )) / rb**3 + rca**3/8.d0
        s5   = dble(data( 61  )) / rb**4
        s6   = dble(data( 62  )) / rb**5 + rca**5/16.d0
        s7   = dble(data( 63  )) / rb**6
        s8   = dble(data( 64  )) / rb**7 + rca**7/25.6d0
        in = 3
        xc=-rb*cos( beta / 57.29578 )
        zc=-rb*sin( beta / 57.29578 )
        call ndipold
10      continue
        bbx=bx
        bby=by
        bbz=bz
        return
        end



        subroutine ndipold
c       mtyp = 3 or 4
c       this version of bfun is mainly for nonuniform field magnets
c         the central field region is represented to 3'rd order on-and-
c         off the midplane by analytic expressions. see slac no. 75
c         fringe field regions represented by fermi type fall-off
c         along with radial fall-off components of 'b' in fringe region
c         evaluated by numerical methods
c
c       The relationship between b0, ...... b12 and b(i,j) relative to
c         axes (z,x) is given by
c
c         b0  = b( 0, 0 )
c         b1  = b( 1, 0 )
c         b2  = b( 2, 0 )
c         b3  = b( 1, 1 )
c         b4  = b( 1,-1 )
c         b5  = b( 0, 1 )
c         b6  = b( 0, 2 )
c         b7  = b( 0,-1 )
c         b8  = b( 0,-2 )
c         b9  = b(-1, 0 )
c         b10 = b(-2, 0 )
c         b11 = b(-1, 1 )
c         b12 = b(-1,-1 )
        implicit real*8(a-h,o-z)
        real*8  ndx, k
        common  /block10old/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(6), dtc(6)
        x = tc(1)
        y = tc(2)
        z = tc(3)
        dx = x - xc
        dz = z - zc
        rp =sqrt( dx**2 + dz**2 )
        dr = rp - rb
        go to ( 1, 2, 3, 14 ), in
7       print 8, in, mtyp
        call exit(0)
8       format ('0 error -go to -  in bfun   in=', i3, '   mtyp=',i4 )
2       drr1 = dr/rb
        drr2 = drr1*drr1
        drr3 = drr2*drr1
        drr4 = drr3*drr1
        if( y .ne. 0. )  go to 4

c       mid-plane uniform field region
        bx = 0.
        by = 0.
        if( mtyp .eq. 3)then
           by=bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 )
        endif
        if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
        bz = 0.
        bt = by
        return

c       non mid-plane uniform field region
4       yr1 = y/rb
        yr2 = yr1*yr1
        yr3 = yr2*yr1
        yr4 = yr3*yr1
        rr1 = rb/rp
        rr2 = rr1*rr1
        rr3 = rr2*rr1
        if( mtyp .eq. 3 ) go to 11
        if( mtyp .eq. 4 ) go to 12
        go to 7

c       mtyp = 3
11      brr = bf*( (-ndx + 2.*bet1*drr1 + 3.*gama*drr2 + 4.*delt*drr3)
     @      * yr1 - (ndx*rr2 + 2.*bet1*rr1*(1.-rr1*drr1)
     @      + 3.*gama*( 2. + 2.*rr1*drr1 - rr2*drr2 )
     @      + 4.*delt*( 6.*drr1 + 3.*rr1*drr2 - rr2*drr3 ))*yr3/6. )
        by = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4
     @      - .5*yr2*( -ndx*rr1 + 2.*bet1*( 1. + rr1*drr1)
     @      + 3.*gama*drr1*( 2. + rr1*drr1) + 4.*delt*drr2*(3.
     @      + rr1*drr1) ) + yr4*( -ndx*rr3 + 2.*bet1*( rr3*drr1 - rr2)
     @      + 3.*gama*( 4.*rr1 - 2.*rr2*drr1 + rr3*drr2 )
     @      + 4.*delt*( 6. + 12.*rr1*drr1 - 3.*rr2*drr2
     @      + rr3*drr3 ) )/24. )
        go to 13

c       mtyp = 4
12      dnr1 = 1. + ndx*drr1
        dnr2 = dnr1*dnr1
        dnr3 = dnr2*dnr1
        dnr4 = dnr3*dnr1
        dnr5 = dnr4*dnr1
        brr = bf*ndx*( -yr1/dnr2 + yr3*( 6.*ndx*ndx/dnr4
     @      - 2.*ndx*rr1/dnr3 - rr2/dnr2 ) /6.  )
        by = bf*( 1./dnr1 + .5*yr2*ndx*( -2.*ndx/dnr3 + rr1/dnr2)
     @      + yr4*ndx*( 24.*ndx**3 /dnr5 - 12.*ndx*ndx*rr1/dnr4
     @      - 2.*ndx*rr2/dnr3 - rr3/dnr2 ) /24.  )
13      bx = brr*dx/rp
        bz = brr*dz/rp
        bt  =sqrt(bx*bx + by*by + bz*bz)
        return
1       sine = -1.
        goto 5
3       sine = 1.
5       if( z  .gt. 0. ) dr = x * sine*csc
        call ndppold( b0, z, x, y, dr      )
        if( y  .ne. 0. )  go to 6

c       mid-plane fringing field region
        bx = 0.
        by = b0
        bz = 0.
        bt   = b0
        return

c       non mid-plane fringing field region
6       if( z .gt. 0. )  go to 9
        dr1  = (sqrt( dx**2 + (dz+dg)**2 ) - rb )
        dr2  = (sqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
        dr3  = (sqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
        dr4  = (sqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
        dr5  = (sqrt( (dx+dg)**2 + dz**2 ) - rb )
        dr6  = (sqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
        dr7  = (sqrt( (dx-dg)**2 + dz**2 ) - rb )
        dr8  = (sqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
        dr9  = (sqrt( dx**2 + (dz-dg)**2 ) - rb )
        dr10 = (sqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
        dr11 = (sqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
        dr12 = (sqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
        go to 10
9       dr1  = sine* x*csc
        dr2  = dr1
        dr9  = dr1
        dr10 = dr1
        dr3  = sine* ( x + dg )*csc
        dr5  = dr3
        dr11 = dr3
        dr4  = sine*( x - dg )*csc
        dr7  = dr4
        dr12 = dr4
        dr6  = sine* ( x + 2.*dg )*csc
        dr8  = sine* ( x - 2.*dg )*csc
10      call ndppold ( b1 , z + dg, x , y , dr1 )
        call ndppold ( b2 , z + 2.*dg, x , y , dr2 )
        call ndppold ( b3 , z + dg, x + dg , y , dr3 )
        call ndppold ( b4 , z + dg, x - dg , y , dr4 )
        call ndppold ( b5 , z , x + dg , y, dr5 )
        call ndppold ( b6 , z , x + 2.*dg , y , dr6 )
        call ndppold ( b7 , z , x - dg , y, dr7 )
        call ndppold ( b8 , z , x - 2.*dg , y , dr8 )
        call ndppold ( b9 , z - dg, x , y , dr9 )
        call ndppold ( b10, z - 2.*dg, x, y, dr10 )
        call ndppold ( b11, z - dg, x + dg , y , dr11 )
        call ndppold ( b12, z - dg, x - dg , y , dr12 )
        yg1 = y/dg
        yg2 = yg1**2
        yg3 = yg1**3
        yg4 = yg1**4
        bx = yg1 * ( (b5-b7)*2./3. - (b6-b8)/12. )
     @      + yg3*( (b5-b7)/6. - (b6-b8)/12.
     @      - (b3 + b11 - b4 - b12 - 2.*b5 + 2.*b7 ) / 12. )
        by = b0 - yg2*( ( b1 + b9 + b5 + b7 - 4.*b0 ) *2./3.
     @      - ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. )
     @      + yg4* (-( b1 + b9 + b5 + b7 - 4.*b0 ) / 6.
     @      + ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24.
     @      + ( b3 + b11 + b4 + b12 - 2.*b1 - 2.*b9
     @      - 2.*b5 - 2.*b7 + 4.*b0 ) / 12. )
        bz = yg1*( (b1 - b9 ) *2./3. - ( b2 - b10 ) /12. )
     @      + yg3*( ( b1 - b9 ) / 6. - ( b2 - b10 ) / 12.
     @      - ( b3 + b4 - b11 - b12 - 2.*b1 + 2.*b9 ) / 12.  )
        bt  =sqrt(bx*bx + by*by + bz*bz)
        return
14      bx = 0.
        by = br
        bz = 0.
        bt = br
        return
        end



        subroutine  ndppold ( bfld, z, x, y , dr )
        implicit real*8(a-h,o-z)
        real*8  ndx, k
        common  /blck10old/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(6), dtc(6)
        drr1 = dr/rb
        drr2 = drr1*drr1
        drr3 = drr2*drr1
        drr4 = drr3*drr1

c       mtyp: modified iterative procedure
        xp = x
        xp2 = xp*xp
        xp3 = xp2*xp
        xp4 = xp3 * xp
        zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2
     @      + s7*xp4*xp3 + s8*xp4*xp4 )
        az = (z-zp)/10.d0
        azmax = sqrt(  x*x + z*z  )
        if( az  .gt.  azmax  )  az = azmax
        zsign = z-zp
        rinv4 = 0.
        do 11 i=1,21
          xp  = x + az*(i-11)
          xp2 = xp*xp
          xp3 = xp2*xp
          xp4 = xp3*xp
          zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2
     @        + s7*xp4*xp3 + s8*xp4*xp4 )
          xxp = x-xp
          zzp = z-zp
          dd = xxp*xxp + zzp*zzp
          if( dd  .lt.  1.d-15 )  dd = 1.d-15
          if( dd  .gt.  1.d15  )  dd = 1.d15
          rinv4 = rinv4 + 1.0d0 / (dd*dd )
11      continue
        dp = sqrt( 1.d0/rinv4 )
        dp = sqrt( dp )
        s = 1.9023d0* sign( 1.d0, zsign ) * dp/d + dels

        cs=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))))
        if( abs(cs)  .gt.  70.  )  cs =sign( 70.d0 ,cs  )
        e=exp(cs)
        p0 = 1.0 + e
        db=bf-br
        bfld = 0.
        if( mtyp .eq. 3 ) bfld
     @    = br +( 1. - ndx*drr1 + bet1*drr2+gama*drr3+delt*drr4)*db/p0
        if( mtyp .eq. 4 ) bfld = br + ( 1./(1. +ndx*drr1) )*db/p0
        return
        end






        subroutine qqfringe(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qqfringe -                                                            *
c      calculates combined fringe fields of two adjacent quads, q1 & q2. *
c      the exit of q1 is located at (x, y, z) and the entrance of q2 is  *
c      dy downstream from the exit of q1.                                *
c      facti describes qi - bquad, bhex, boct, bdec, bddec               *
c      qradi radius of qi                                                *
c*************************************************************************
        real fact1(5),fact2(5),b(3)
        bx=0.
        by=0.
        bz=0.

c       exit of q1 field
        xdum=x
        ydum=y
        zdum=z
        call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
        bx=b(1)
        by=b(2)
        bz=b(3)

c       entrance of q2 field
        ydum=dy-ydum
        call mpoles(1,xdum,zdum,ydum,fact2,qrad2,b(1),b(3),b(2))
        bx=bx+b(1)
        by=by-b(2)
        bz=bz+b(3)
        return
        end




        subroutine qdfringe(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qdfringe -                                                            *
c      calculates combined fringe fields of one quad, q1, followed by a  *
c      dipole, d.                                                        *
c      the exit of q1 is located at x=y=z=0 and the entrance of d is     *
c      dy downstream from the exit of q1 and rotated by angle th.        *
c      fact1 describes q1 - bquad, bhex, boct, bdec, bddec               *
c      fact2 describes dipole - bdip, theta
c*************************************************************************
        real fact1(5),fact2(5),b(3),bpri(3)
        data pifac/1.745329e-02/
        th=fact2(2)*pifac
        bx=0.
        by=0.
        bz=0.

c       exit of q1 field
        xdum=x
        ydum=y
        zdum=z
        call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
        bx=b(1)
        by=-b(2)
        bz=b(3)

c       entrance dipole field
c         translate and reflect y axis
        ypri=y-dy
        xpri=x
        zpri=z

c       rotate by th
        xdpri=(xpri*cos(th))+(ypri*sin(th))
        ydpri=(-xpri*sin(th))+(ypri*cos(th))
        zdpri=z
        ydpri=-ydpri
        call dfringe(xdpri,ydpri,zdpri,fact2(1),qrad2,b(1),b(2),b(3))
        b(2)=-b(2)
        bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
        bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
        bpri(3)=b(3)
        bx=bx+bpri(1)
        by=by+bpri(2)
        bz=bz+bpri(3)
        return
        end



        subroutine ddfringe(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  ddfringe -                                                            *
c      calculates combined fringe fields of two dipoles                  *
c      both d1 and d2 pole faces are rotated by thi with respect to the  *
c      optic axis making a total bend from the center of d1 to the center*
c      of d2 of th=th1+th2. (see diagrams in notebook 29 april 1987)     *
c      the exit of d1 is located at x=y=z=0 and the entrance of d2 is    *
c      dl downstream from the exit of d1 and rotated by angle th.        *
c      fact1 describes d1 - bdip, theta1                                 *
c      fact2 describes d2 - bdip, theta2                                 *
c*************************************************************************
        real fact1(5),fact2(5),b(3)
        data pifac/1.745329e-02/
        th2=fact2(2)*pifac
        th1=fact1(2)*pifac
        bx=0.
        by=0.
        bz=0.

c       exit of d1 field
        xdum=x
        ydum=y
        zdum=z
        call dfringe(xdum,ydum,zdum,fact1(1),qrad,b(1),b(2),b(3))
        by=-by

c       entrance dipole field
c         translate x & y axes
        dy=dl*sin(th1)
        dx=dl*cos(th1)
        ypri=y-dy
        xpri=x+dx
        zpri=z

c       rotate by th
        th=th1+th2
        xdpri=(xpri*cos(th))+(ypri*sin(th))
        ydpri=(-xpri*sin(th))+(ypri*cos(th))
        zdpri=z
        call dfringe(xdpri,ydpri,zdpri,fact2(1),qrad,b(1),b(2),b(3))
        bx=bx+b(1)
        by=by+b(2)
        bz=bz+b(3)
        return
        end




        subroutine dqfringe(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c***************************************************************************
c    dqfringe - computes overlap fringe fields for a dipole quadrupole     *
c    pair.  the dipole center exit is located at x=y=z=0. the optic axis   *
c    makes an angle th with the dipole face.  the quad is located a        *
c    distance dl along the optic axis from the dipole.                     *
c    diagrams in notebook 29 april 1987                                    *
c                                                                          *
c      fact1 describes dipole - bdip, th                                   *
c      fact2 describes quad - bquad, bhex, boct, bdec, bddec               *
c***************************************************************************
        real fact1(5),fact2(5),b(3)
        data pifac/1.745329e-02/
        th=fact1(2)*pifac
        bx=0.
        by=0.
        bz=0.

c       exit of d1 field
        xdum=x
        ydum=y
        zdum=z
        call dfringe(xdum,ydum,zdum,fact1(1),qrad1,bx,by,bz)
        by=-by

c       entrance quad field
c         translate x & y axes
        dy=dl*sin(th)
        dx=dl*cos(th)
        ypri=y-dy
        xpri=x+dx
        zpri=z

c       rotate by th
        xdpri=(xpri*cos(th))+(ypri*sin(th))
        ydpri=(-xpri*sin(th))+(ypri*cos(th))
        zdpri=z
        ydpri=-ydpri
        call mpoles(1,xdpri,zdpri,ydpri,fact2,qrad2,b(1),b(3),b(2))
        bx=bx+b(1)
        by=by+b(2)
        bz=bz+b(3)
        return
        end




        subroutine dqfringe2(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c***************************************************************************
c    dqfringe - computes overlap fringe fields for a dipole quadrupole     *
c    pair.  the dipole center exit is located at x=y=z=0. the optic axis   *
c    makes an angle th with the dipole face.  the quad is located a        *
c    distance dl along the optic axis from the dipole.                     *
c    diagrams in notebook 29 april 1987                                    *
c    modified version of dqfringe. altered to include dipolert 8/23/89 jjl *
c                                                                          *
c      fact1 describes dipole - bdip, th                                   *
c      fact2 describes quad - bquad, bhex, boct, bdec, bddec               *
c***************************************************************************
        real fact1(5),fact2(5),b(3),bpri(3)
        data pifac/1.745329e-02/,ten/10./
        th=fact1(2)*pifac
        bx=0.
        by=0.
        bz=0.

c       exit of d1 field
        xdum=x/ten
        ydum=y/ten
        zdum=z/ten
        call dipolert2(xdum,zdum,ydum,bx,bz,by)

c       get the right sign on the field
        bx=-bx
        by=-by
        bz=-bz

c       entrance quad field
c         translate x & y axes
        dx=dl*sin(th)
        dy=dl*cos(th)
        ypri=y-dy
        xpri=x+dx
        zpri=z

c       rotate by th
        xdpri=(xpri*cos(th))+(ypri*sin(th))
        ydpri=(-xpri*sin(th))+(ypri*cos(th))
        zdpri=z
        ydpri=-ydpri
        call mpoles(1,xdpri,zdpri,ydpri,fact2,qrad2,b(1),b(3),b(2))
        b(2)=-b(2)
        bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
        bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
        bpri(3)=b(3)
        bx=bx+bpri(1)
        by=by+bpri(2)
        bz=bz+bpri(3)
        return
        end




        subroutine qdfringe2(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qdfringe -                                                            *
c      calculates combined fringe fields of one quad, q1, followed by a  *
c      dipole, d.                                                        *
c      the exit of q1 is located at x=y=z=0 and the entrance of d is     *
c      dy downstream from the exit of q1 and rotated by angle th.        *
c      fact1 describes q1 - bquad, bhex, boct, bdec, bddec               *
c      fact2(1)=dx
c      fact2(2)=dy
c      fact2(3)=th
c*************************************************************************
        real fact1(5),fact2(5),b(3),bpri(3)
        data pifac/1.745329e-02/,ten/10./
        th=fact2(3)*pifac
        dx=fact2(1)
        dy=fact2(2)
        bx=0.
        by=0.
        bz=0.

c       exit of q1 field
        xdum=x
        ydum=y
        zdum=z
        call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
        bx=b(1)
        by=b(2)
        bz=b(3)

c       entrance dipole field
c         translate and reflect y axis
        ypri=y-dy
        xpri=x-dx
        zpri=z

c       rotate by th
        xdpri=((xpri*cos(th))+(ypri*sin(th)))/ten
        ydpri=((-xpri*sin(th))+(ypri*cos(th)))/ten
        zdpri=z/ten
        call dipolert2(xdpri,zdpri,ydpri,b(1),b(3),b(2))

c       get the right sign on the field
        do 10 i=1,3
10        b(i)=-b(i)
        bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
        bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
        bpri(3)=b(3)
        bx=bx+bpri(1)
        by=by+bpri(2)
        bz=bz+bpri(3)
        return
        end




        subroutine dipolert2(x,y,z,bbx,bby,bbz)
c       dipolert2- analytic calculation of dipole magnetic fields
c         stolen from raytrace adapted for use in snake
c
c       data(i) are read in snake with ifield
c       dipolert2 is called from snake in analyf
c         snake provides x,y,z in the c-axis system in cm
c         dipolert2 returns bx,by,bz in the same system
        implicit real*8(a-h,o-z)
        real x,y,z,bbx,bby,bbz,data(75),pii
        real*8  ndx, k
        common  /diprt/data
        common  /blck10/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(2,6), dtc(6), bb(2,0:12)
        data pii/0.0174532925/
        lf1  = dble(data(  1  ))
        lu1  = dble(data(  2  ))
        lf2  = dble(data(  3  ))
        dg   = dble(data(  4  ))
        mtyp = dble(data(  5  ))
        a    = dble(data(  11  ))
        b    = dble(data(  12  ))
        d    = dble(data(  13  ))
        rb   = dble(data(  14  ))
        bf   = dble(data(  15  ))
        phi  = dble(data(  16  ))
        alpha= dble(data(  17  ))
        beta = dble(data(  18  ))
        ndx  = dble(data(  19  ))
        bet1 = dble(data(  20  ))
        gama = dble(data(  21  ))
        delt = dble(data(  22  ))
        z11  = dble(data(  25  ))
        z12  = dble(data(  26  ))
        z21  = dble(data(  27  ))
        z22  = dble(data(  28  ))
        br1  = dble(data(  41  ))
        br2  = dble(data(  42  ))
        xcr1 = dble(data(  43  ))
        xcr2 = dble(data(  44  ))
        xpric=x+(2.*rb*dsin(pii*phi/2.)*dsin(pii*((phi/2.)-beta)))
        zpric=z+(2.*rb*dsin(pii*phi/2.)*dcos(pii*((phi/2.)-beta)))
        cosa=dcos(pii*(phi-alpha-beta))
        sina=dsin(pii*(phi-alpha-beta))
        tc(1,1)=(-xpric*cosa)+(zpric*sina)
        tc(1,2)=dble(y)
        tc(1,3)=(-xpric*sina)-(zpric*cosa)

c       transform to second vfb coord system
        copab =dcos( (phi-alpha-beta)/57.29578)
        sipab =dsin( (phi-alpha-beta)/57.29578)
        cospb =dcos( (phi/2.-beta)/57.29578 )
        sinpb =dsin( (phi/2.-beta)/57.29578 )
        sip2 =dsin( (phi/2.)/57.29578 )
        xt = tc(1,1)
        zt = tc(1,3)
        vxt = tc(1,4)
        vzt = tc(1,6)
        tc(2,3) = - zt  *copab +  xt  *sipab -2.*rb*sip2*cospb
        tc(2,1) = - zt  *sipab -  xt  *copab -2.*rb*sip2*sinpb
        tc(2,6) = - vzt *copab +  vxt *sipab
        tc(2,4) = - vzt *sipab -  vxt *copab
        tc(2,2)=dble(y)

c       test for region
        if(tc(1,3).gt.z11)then
          bx=0.
          by=0.
          bz=0.
          go to 10
        endif
        if((tc(1,3).ge.z12).and.(tc(2,3).lt.z21)) go to 1  !entrance fringe
        if((tc(1,3).ge.z12).and.(tc(2,3).ge.z21)) go to 25 !overlapping entr&exit
        if(tc(1,3).lt.z12) go to 2                         !uniform field or exit

c       in designates magnet regions for bfun
1       in = 1
        xc= rb*dcos( alpha/ 57.29578 )
        zc=-rb*dsin( alpha/ 57.29578 )
        c0   = dble(data(  29  ))
        c1   = dble(data(  30  ))
        c2   = dble(data(  31  ))
        c3   = dble(data(  32  ))
        c4   = dble(data(  33  ))
        c5   = dble(data(  34  ))
        dels = dble(data(  45  ))
        rca  = dble(data(  47  ))
        csc = dcos( alpha/57.29578 )
        scor = dble(data( 49 ))
        s2   = dble(data(  51  )) / rb    + rca/2.d0
        s3   = dble(data(  52  )) / rb**2
        s4   = dble(data(  53  )) / rb**3 + rca**3/8.d0
        s5   = dble(data(  54  )) / rb**4
        s6   = dble(data(  55  )) / rb**5 + rca**5/16.d0
        s7   = dble(data(  56  )) / rb**6
        s8   = dble(data(  57  )) / rb**7 + rca**7/25.6d0
        call ndip(1)

c       transform to c-axis system
        bxdum=bx
        bzdum=bz
        bx=(-bxdum*cosa)-(bzdum*sina)
        bz=(+bxdum*sina)-(bzdum*cosa)
        go to 10
2       continue
        if(tc(2,3).lt.z21)go to 3      !uniform field
        if(tc(2,3).gt.z22)then         !outside of fringe region
          bx=0.
          by=0.
          bz=0.
          go to 10
        endif
        if(tc(2,3).le.z22)go to 4      !exit fringe

c       uniform field integration region
3       in = 2
        xc= rb*dcos( beta / 57.29578 )
        zc=-rb*dsin( beta / 57.29578 )
        call ndip(1)
        go to 10

c       setup for second fringe field and integration
4       br   = br2
        c0   = dble(data(  35  ))
        c1   = dble(data(  36  ))
        c2   = dble(data(  37  ))
        c3   = dble(data(  38  ))
        c4   = dble(data(  39  ))
        c5   = dble(data(  40  ))
        dels = dble(data(  46  ))
        rca  = dble(data(  48  ))
        scor = dble(data( 50 ))
        csc = dcos( beta /57.29578 )
        s2   = dble(data(  58  )) / rb    + rca/2.d0
        s3   = dble(data(  59  )) / rb**2
        s4   = dble(data(  60  )) / rb**3 + rca**3/8.d0
        s5   = dble(data(  61  )) / rb**4
        s6   = dble(data(  62  )) / rb**5 + rca**5/16.d0
        s7   = dble(data(  63  )) / rb**6
        s8   = dble(data(  64  )) / rb**7 + rca**7/25.6d0
        in = 3
        xc=-rb*dcos( beta / 57.29578 )
        zc=-rb*dsin( beta / 57.29578 )
        call ndip(2)
10      continue
        bbx=bx
        bby=by
        bbz=bz
        return

c       overlapping entrance and exit fringe fields
25      in = 1                         !setup for entrance field
        xc= rb*dcos( alpha/ 57.29578 )
        zc=-rb*dsin( alpha/ 57.29578 )
        c0   = dble(data(  29  ))
        c1   = dble(data(  30  ))
        c2   = dble(data(  31  ))
        c3   = dble(data(  32  ))
        c4   = dble(data(  33  ))
        c5   = dble(data(  34  ))
        dels = dble(data(  45  ))
        rca  = dble(data(  47  ))
        csc = dcos( alpha/57.29578 )
        scor = dble(data(49 ))
        s2   = dble(data(  51  )) / rb    + rca/2.d0
        s3   = dble(data(  52  )) / rb**2
        s4   = dble(data(  53  )) / rb**3 + rca**3/8.d0
        s5   = dble(data(  54  )) / rb**4
        s6   = dble(data(  55  )) / rb**5 + rca**5/16.d0
        s7   = dble(data(  56  )) / rb**6
        s8   = dble(data(  57  )) / rb**7 + rca**7/25.6d0
        call nndip(1,bb)

c       setup for exit fringe field
        br   = br2
        c0   = dble(data(  35  ))
        c1   = dble(data(  36  ))
        c2   = dble(data(  37  ))
        c3   = dble(data(  38  ))
        c4   = dble(data(  39  ))
        c5   = dble(data(  40  ))
        dels = dble(data(  46  ))
        rca  = dble(data(  48  ))
        scor = dble(data(50 ))
        csc = dcos( beta /57.29578 )
        s2   = dble(data(  58  )) / rb    + rca/2.d0
        s3   = dble(data(  59  )) / rb**2
        s4   = dble(data(  60  )) / rb**3 + rca**3/8.d0
        s5   = dble(data(  61  )) / rb**4
        s6   = dble(data(  62  )) / rb**5 + rca**5/16.d0
        s7   = dble(data(  63  )) / rb**6
        s8   = dble(data(  64  )) / rb**7 + rca**7/25.6d0
        in = 3
        xc=-rb*dcos( beta / 57.29578 )
        zc=-rb*dsin( beta / 57.29578 )
        call nndip(2,bb)
        if(y.eq.0.)then
          bx=0.
          bz=0.
          by=bb(1,0)+bb(2,0)-bf
        else
          do 26 i=0,12
            bb(1,i)=bb(1,i)+bb(2,i)-bf
26        continue
          yg1 = y/dg
          yg2 = yg1**2
          yg3 = yg1**3
          yg4 = yg1**4
          bx = yg1 * ( (bb(1,5)-bb(1,7))*2./3. - (bb(1,6)-bb(1,8))
     @        / 12. ) + yg3*( (bb(1,5)-bb(1,7))/6. - (bb(1,6)-bb(1,8))
     @        / 12. - (bb(1,3) + bb(1,11) - bb(1,4) - bb(1,12)
     @        - 2.*bb(1,5) + 2.*bb(1,7) ) / 12. )
          by = bb(1,0) - yg2*( ( bb(1,1) + bb(1,9) + bb(1,5) + bb(1,7)
     @        - 4.*bb(1,0) ) *2./3. - ( bb(1,2) + bb(1,10) + bb(1,6)
     @        + bb(1,8)-4.*bb(1,0))/24.) + yg4*(-(bb(1,1)+ bb(1,9)
     @        + bb(1,5)+ bb(1,7)- 4.*bb(1,0) )/6.+ ( bb(1,2)+ bb(1,10)
     @        + bb(1,6) + bb(1,8) - 4.*bb(1,0))/24. + ( bb(1,3)+bb(1,11)
     @        + bb(1,4) +bb(1,12)-2.*bb(1,1)-2.*bb(1,9)-2.*bb(1,5) - 2.
     @        * bb(1,7) + 4.*bb(1,0) ) / 12. )
          bz = yg1*((bb(1,1)-bb(1,9))*2./3. - (bb(1,2)-bb(1,10) ) /12. )
     @        + yg3*( (bb(1,1)- bb(1,9))/6. - (bb(1,2) - bb(1,10) )
     @        / 12. - ( bb(1,3) + bb(1,4) - bb(1,11) - bb(1,12)
     @        - 2.*bb(1,1) + 2.*bb(1,9) ) / 12.  )
          bt  =dsqrt(bx*bx + by*by + bz*bz)
        endif
        go to 10
        end




        subroutine ndip(l)
c       mtyp = 3 or 4
c       this version of bfun is mainly for nonuniform field magnets
c         the central field region is represented to 3'rd order on-and-
c         off the midplane by analytic expressions. see slac no. 75
c         fringe field regions represented by fermi type fall-off
c         along with radial fall-off
c       components of 'b' in fringe region evaluated by numerical methods
c       the relationship between b0, ......... b12 and b(i,j) relative to
c         axes (z,x) is given by
c         b0  = b( 0, 0 )
c         b1  = b( 1, 0 )
c         b2  = b( 2, 0 )
c         b3  = b( 1, 1 )
c         b4  = b( 1,-1 )
c         b5  = b( 0, 1 )
c         b6  = b( 0, 2 )
c         b7  = b( 0,-1 )
c         b8  = b( 0,-2 )
c         b9  = b(-1, 0 )
c         b10 = b(-2, 0 )
c         b11 = b(-1, 1 )
c         b12 = b(-1,-1 )
        implicit real*8(a-h,o-z)
        real*8  ndx, k
        common  /blck10/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(2,6), dtc(6)
        x = tc(l,1)
        y = tc(l,2)
        z = tc(l,3)
        dx = x - xc
        dz = z - zc
        rp =dsqrt( dx**2 + dz**2 )
        dr = rp - rb
        go to ( 1, 2, 3, 14 ), in
7       print 8, in, mtyp
        call exit(0)
8       format ('0 error -go to -  in bfun   in=', i3, '   mtyp=',i4 )
2       drr1 = dr/rb
        drr2 = drr1*drr1
        drr3 = drr2*drr1
        drr4 = drr3*drr1
        if( y .ne. 0. )  go to 4

c       mid-plane uniform field region
        bx = 0.
        by = 0.
        if( mtyp .eq. 3) by
     @      = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4)
        if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
        bz = 0.
        bt = by
        return

c       non mid-plane uniform field region
4       yr1 = y/rb
        yr2 = yr1*yr1
        yr3 = yr2*yr1
        yr4 = yr3*yr1
        rr1 = rb/rp
        rr2 = rr1*rr1
        rr3 = rr2*rr1
        if( mtyp .eq. 3 ) go to 11
        if( mtyp .eq. 4 ) go to 12
        go to 7

c       mtyp = 3
11      brr = bf*( (-ndx + 2.*bet1*drr1 + 3.*gama*drr2 + 4.*delt*drr3)
     @      * yr1 - (ndx*rr2 + 2.*bet1*rr1*(1.-rr1*drr1)
     @      + 3.*gama*( 2. + 2.*rr1*drr1 - rr2*drr2 )
     @      + 4.*delt*( 6.*drr1 + 3.*rr1*drr2 - rr2*drr3 ))*yr3/6. )
        by = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4
     @      - .5*yr2*( -ndx*rr1 + 2.*bet1*( 1. + rr1*drr1)
     @      + 3.*gama*drr1*( 2. + rr1*drr1) + 4.*delt*drr2
     @      * (3. + rr1*drr1) )+yr4*( -ndx*rr3 + 2.*bet1
     @      * ( rr3*drr1 - rr2)+3.*gama*( 4.*rr1 - 2.*rr2*drr1
     @      + rr3*drr2 )+4.*delt*( 6. + 12.*rr1*drr1 - 3.*rr2*drr2
     @      + rr3*drr3 ) )/24. )
        go to 13

c       mtyp = 4
12      dnr1 = 1. + ndx*drr1
        dnr2 = dnr1*dnr1
        dnr3 = dnr2*dnr1
        dnr4 = dnr3*dnr1
        dnr5 = dnr4*dnr1
        brr = bf*ndx*( -yr1/dnr2 + yr3*( 6.*ndx*ndx/dnr4
     @      - 2.*ndx*rr1/dnr3 - rr2/dnr2 ) /6.  )
        by = bf*( 1./dnr1 + .5*yr2*ndx*( -2.*ndx/dnr3 + rr1/dnr2)
     @      + yr4*ndx*( 24.*ndx**3 /dnr5 - 12.*ndx*ndx*rr1/dnr4
     @      - 2.*ndx*rr2/dnr3 - rr3/dnr2 ) /24.  )

13      bx = brr*dx/rp
        bz = brr*dz/rp
        bt  =dsqrt(bx*bx + by*by + bz*bz)
        return
1       sine = -1.
        go to 5
3       sine = 1.
5       if( z  .gt. 0. ) dr = x * sine*csc
        call ndpp( b0, z, x, y, dr      )
        if( y  .ne. 0. )  go to 6

c       mid-plane fringing field region
        bx = 0.
        by = b0
        bz = 0.
        bt   = b0
        return

c       non mid-plane fringing field region
6       if( z .gt. 0. )  go to 9
        dr1  = (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
        dr2  = (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
        dr3  = (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
        dr4  = (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
        dr5  = (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
        dr6  = (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
        dr7  = (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
        dr8  = (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
        dr9  = (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
        dr10 = (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
        dr11 = (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
        dr12 = (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
        go to 10
9       dr1  = sine* x*csc
        dr2  = dr1
        dr9  = dr1
        dr10 = dr1
        dr3  = sine* ( x + dg )*csc
        dr5  = dr3
        dr11 = dr3
        dr4  = sine*( x - dg )*csc
        dr7  = dr4
        dr12 = dr4
        dr6  = sine* ( x + 2.*dg )*csc
        dr8  = sine* ( x - 2.*dg )*csc
10      call ndpp ( b1 , z + dg, x , y , dr1 )
        call ndpp ( b2 , z + 2.*dg, x , y , dr2 )
        call ndpp ( b3 , z + dg, x + dg , y , dr3 )
        call ndpp ( b4 , z + dg, x - dg , y , dr4 )
        call ndpp ( b5 , z , x + dg , y, dr5 )
        call ndpp ( b6 , z , x + 2.*dg , y , dr6 )
        call ndpp ( b7 , z , x - dg , y, dr7 )
        call ndpp ( b8 , z , x - 2.*dg , y , dr8 )
        call ndpp ( b9 , z - dg, x , y , dr9 )
        call ndpp ( b10, z - 2.*dg, x, y, dr10 )
        call ndpp ( b11, z - dg, x + dg , y , dr11 )
        call ndpp ( b12, z - dg, x - dg , y , dr12 )
        yg1 = y/dg
        yg2 = yg1**2
        yg3 = yg1**3
        yg4 = yg1**4

c       fudge to fix dipole fringe field -4/8/02  jjl
        fudge=1.30986
        bx = fudge*yg1 * ( (b5-b7)*2./3. - (b6-b8)/12. )
     @      + yg3*( (b5-b7)/6. - (b6-b8)/12.
     @      - (b3 + b11 - b4 - b12 - 2.*b5 + 2.*b7 ) / 12. )
        by = b0 - yg2*( ( b1 + b9 + b5 + b7 - 4.*b0 ) *2./3.
     @      - ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. )
     @      + yg4* (-( b1 + b9 + b5 + b7 - 4.*b0 ) / 6.
     @      + ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24.
     @      + ( b3 + b11 + b4 + b12 - 2.*b1 - 2.*b9
     @      - 2.*b5 - 2.*b7 + 4.*b0 ) / 12. )
        bz = fudge*yg1*( (b1 - b9 ) *2./3. - ( b2 - b10 ) /12. )
     @      + yg3*( ( b1 - b9 ) / 6. - ( b2 - b10 ) / 12.
     @      - ( b3 + b4 - b11 - b12 - 2.*b1 + 2.*b9 ) / 12.  )
        bt  =dsqrt(bx*bx + by*by + bz*bz)
        return
14      bx = 0.
        by = br
        bz = 0.
        bt = br
        return
        end





        subroutine nndip(l,bb)
c       mtyp = 3 or 4
c       this version of bfun is mainly for nonuniform field magnets
c         the central field region is represented to 3'rd order on-and-
c         off the midplane by analytic expressions. see slac no. 75
c         fringe field regions represented by fermi type fall-off
c         along with radial fall-off
c       components of 'b' in fringe region evaluated by numerical methods
c       the relationship between b0, ......... b12 and b(i,j) relative to
c         axes (z,x) is given by
c         bb(l,0)  = b( 0, 0 )
c         bb(l,1)  = b( 1, 0 )
c         bb(l,2)  = b( 2, 0 )
c         bb(l,3)  = b( 1, 1 )
c         bb(l,4)  = b( 1,-1 )
c         bb(l,5)  = b( 0, 1 )
c         bb(l,6)  = b( 0, 2 )
c         bb(l,7)  = b( 0,-1 )
c         bb(l,8)  = b( 0,-2 )
c         bb(l,9)  = b(-1, 0 )
c         bb(l,10) = b(-2, 0 )
c         bb(l,11) = b(-1, 1 )
c         bb(l,12) = b(-1,-1 )
c       modified for use with snake in order to deal with everlapping
c         entrance and exit fringing fields
        implicit real*8(a-h,o-z)
        real*8  ndx, k
        common  /blck10/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(2,6), dtc(6),bb(2,0:12)
        x = tc(l,1)
        y = tc(l,2)
        z = tc(l,3)
        dx = x - xc
        dz = z - zc
        rp =dsqrt( dx**2 + dz**2 )
        dr = rp - rb
        go to ( 1, 7, 3, 7 ), in
7       print 8, in, mtyp
        call exit(0)
8       format ('0 error -go to -  in nndip   in=', i3, '   mtyp=',i4)
1       sine = -1.
        go to 5
3       sine = 1.
5       if( z  .gt. 0. ) dr = x * sine*csc
        call ndpp( bb(l,0), z, x, y, dr      )
        if( y  .ne. 0. )  go to 6

c       mid-plane fringing field region
        return

c       non mid-plane fringing field region
6       if( z .gt. 0. )  go to 9
        dr1  =       (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
        dr2  =       (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
        dr3  =       (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
        dr4  =       (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
        dr5  =       (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
        dr6  =       (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
        dr7  =       (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
        dr8  =       (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
        dr9  =       (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
        dr10 =       (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
        dr11 =       (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
        dr12 =       (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
        go to 10
9       dr1  = sine* x*csc
        dr2  = dr1
        dr9  = dr1
        dr10 = dr1
        dr3  = sine* ( x + dg )*csc
        dr5  = dr3
        dr11 = dr3
        dr4  = sine*( x - dg )*csc
        dr7  = dr4
        dr12 = dr4
        dr6  = sine* ( x + 2.*dg )*csc
        dr8  = sine* ( x - 2.*dg )*csc
10      call ndpp ( bb(l,1) , z + dg, x , y , dr1 )
        call ndpp ( bb(l,2) , z + 2.*dg, x , y , dr2 )
        call ndpp ( bb(l,3) , z + dg, x + dg , y , dr3 )
        call ndpp ( bb(l,4) , z + dg, x - dg , y , dr4 )
        call ndpp ( bb(l,5) , z , x + dg , y, dr5 )
        call ndpp ( bb(l,6) , z , x + 2.*dg , y , dr6 )
        call ndpp ( bb(l,7) , z , x - dg , y, dr7 )
        call ndpp ( bb(l,8) , z , x - 2.*dg , y , dr8 )
        call ndpp ( bb(l,9) , z - dg, x , y , dr9 )
        call ndpp ( bb(l,10), z - 2.*dg, x, y, dr10 )
        call ndpp ( bb(l,11), z - dg, x + dg , y , dr11 )
        call ndpp ( bb(l,12), z - dg, x - dg , y , dr12 )
        return
        end




        subroutine  ndpp ( bfld, z, x, y , dr )
        implicit real*8(a-h,o-z)
        real*8  ndx, k
        common  /blck10/  bx, by, bz, k, tc, dtc
        common  /blck20/  ndx,bet1,gama,delt,csc
        common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
        common  /blck22/  d, dg, s, bf, bt
        common  /blck23/  c0, c1, c2, c3, c4, c5
        common  /blck24/  rb, xc, zc
        common  /blck25/  in, mtyp
        dimension tc(2,6), dtc(6)
        drr1 = dr/rb
        drr2 = drr1*drr1
        drr3 = drr2*drr1
        drr4 = drr3*drr1

c       mtyp: modified iterative procedure
        xp = x
        xp2 = xp*xp
        xp3 = xp2*xp
        xp4 = xp3 * xp
        zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2
     @      + s7*xp4*xp3 + s8*xp4*xp4 )
        az = (z-zp)/10.d0
        azmax = dsqrt(  x*x + z*z  )
        if( az  .gt.  azmax  )  az = azmax
        zsign = z-zp
        rinv4 = 0.
        do 11 i=1,21
          xp   = x + az*(i-11)
          xp2 = xp*xp
          xp3 = xp2*xp
          xp4 = xp3*xp
          zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2
     @        + s7*xp4*xp3 + s8*xp4*xp4 )
          xxp = x-xp
          zzp = z-zp
          dd = xxp*xxp + zzp*zzp
          if( dd  .lt.  1.d-15 )  dd = 1.d-15
          if( dd  .gt.  1.d15  )  dd = 1.d15
          rinv4 = rinv4 + 1.0d0 / (dd*dd )
11      continue
        dp = dsqrt( 1.d0/rinv4 )
        dp = dsqrt( dp )
        s = 1.9023d0* dsign( 1.d0, zsign ) * dp/d + dels
        cs=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))))
        if( dabs(cs)  .gt.  70.  )  cs =dsign( 70.d0 ,cs  )
        e=dexp(cs)
        p0 = 1.0 + e
        db=bf-br
        bfld = 0.
        if( mtyp .eq. 3 ) bfld
     @    = br +( 1. - ndx*drr1 + bet1*drr2+gama*drr3+delt*drr4)*db/p0
        if( mtyp .eq. 4 ) bfld = br + ( 1./(1. +ndx*drr1) )*db/p0
        return
        end




        subroutine hlmhltz(x,y,z,r,amp,bx,by,bz)
c************************************************************************
c hlmhltz - calculates field for a pair of helmholtz coils
c           x,y,z = field point coords, origin exactly between coils (meters)
c               r = coil radius and separation
c             amp = current in coils
c              bi = x,y,z components of field (tesla)
c 7/8/87
c*************************************************************************
        real x,y,z,r,amp,bx,by,bz,rho,br1,br2,by1,by2,theta,y1,y2
        y1=(r/2)+y
        y2=y-(r/2)
        rho=sqrt(x**2+z**2)
        if (rho.eq.0.)then
          theta=0.
          go to 1
        endif
        theta=atan2(z,x)
1       call loop(r,amp,rho,y1,br1,by1)
        call loop(r,amp,rho,y2,br2,by2)
        by=by1+by2
        br=br1+br2
        bx=br*cos(theta)
        bz=br*sin(theta)
        return
        end




        subroutine solenoid(x,y,z,r,rl,amp,bx,by,bz)
c************************************************************************
c solenoid - calculates approximate field for a solenoid (sum over many current loops)
c           x,y,z = field point coords, origin center of coil (meters)
c               r = coil radius
c              rl = length
c             amp = current in coils
c              bi = x,y,z components of field (tesla)
c 1/7/05
c*************************************************************************
        parameter (loops=1000) !even number of loops
        real x,y,z,r,rl,amp,bx,by,bz,rho,br1,by1,theta,yy(loops)
        dl=rl/(loops-1) !space between loops
        loops2 = loops/2
        bx = 0.
        by = 0.
        bz = 0.
        br = 0.
        do 10 i=1,loops
        j=i-1
        yy(i)=(dl*(j-loops2))+y
        rho=sqrt(x**2+z**2)
        if (rho.eq.0.)then
        theta=0.
        go to 1
        endif
        theta=atan2(z,x)
1       call loop(r,amp,rho,yy(i),br1,by1)
        by=by+by1
        bx=bx+br1*cos(theta)
        bz=Bz+br1*sin(theta)
10      continue
        return
        end




        subroutine loop(a,amp,rho,y,br,by)
c************************************************************************
c  loop - calculates the magnetic field for a current loop
c         a = radius of loop
c       amp = current in loop
c       rho = distance from loop axis (meters)
c         y = distance from loop along axis (meters)
c        br = radial component of field (tesla)
c        by = axial component of field
c
c  reference: garrett, journal of applied physics, 34 sept 1963, p. 2567
c  7/8/87
c******************************************************************
        real a,amp,rho,y,br,by,r1sq,r2sq,ksq,kpsq
        real*8 dksq,dvk,dve

c       evaluate parameters
        r1sq=(a+rho)**2+y**2
        r2sq=(a-rho)**2+y**2
        ksq=4.*a*rho/r1sq
        kpsq=1-ksq

c       evaluate elliptical integrals
        dksq=ksq
        call fb01ad(dksq,dvk,dve)

c       calculate fields
        if(r2sq.eq.0.)then
          write(6,*)' bull''s-eye!!! you just hit a coil.'
          by=1.69e+38
          go to 15
        endif
        by=(amp*2.0e-07/sqrt(r1sq))*((dvk-dve)+(a-rho)*2*a*dve/r2sq)
15      if(rho.eq.0.)then
          br=0.
          return
        endif
        if(kpsq.eq.0.)kpsq=0.3e-36
        br=-(amp*1.e-07*y)/(rho*sqrt(r1sq)*kpsq)*((2-ksq)*(dvk-dve)
     @    -ksq*dvk)
        return
        end





        subroutine fb01ad(c,  vk,ve)
        implicit real*8(a-h,o-z)
        real * 8 xlg
        data xlg/'7fffffffffffffff'x/
        d=1d0-c
        if(d .gt. 0d0)e=-log(d)
        if(c .ge. 1d0)go to 2
        ve=e*((((((((((
     @    3.18591956555015718d-5*d  +.989833284622538479d-3)*d
     @    +.643214658643830177d-2)*d +.16804023346363385d-1)*d
     @    +.261450147003138789d-1)*d +.334789436657616262d-1)*d
     @    +.427178905473830956d-1)*d +.585936612555314917d-1)*d
     @    +.937499997212031407d-1)*d +.249999999999901772d0)*d)
     @    +(((((((((
     @    .149466217571813268d-3*d  +.246850333046072273d-2)*d
     @    +.863844217360407443d-2)*d+.107706350398664555d-1)*d
     @    +.782040406095955417d-2)*d +.759509342255943228d-2)*d
     @    +.115695957452954022d-1)*d +.218318116761304816d-1)*d
     @    +.568051945675591566d-1)*d +.443147180560889526d0)*d
     @    +1d0

c**** routine modified to calculate vk and ve always
        vk=e*((((((((((
     @    .297002809665556121d-4*d   +.921554634963249846d-3)*d
     @    +.597390429915542916d-2)*d  +.155309416319772039d-1)*d
     @    +.239319133231107901d-1)*d  +.301248490128989303d-1)*d
     @    +.373777397586236041d-1)*d  +.48828041906862398d-1)*d
     @    +.703124997390383521d-1)*d  +.124999999999908081d0)*d
     @    +.5d0)+(((((((((
     @    .139308785700664673d-3*d   +.229663489839695869d-2)*d
     @    +.800300398064998537d-2)*d  +.984892932217689377d-2)*d
     @    +.684790928262450512d-2)*d  +.617962744605331761d-2)*d
     @    +.878980187455506468d-2)*d  +.149380135326871652d-1)*d
     @    +.308851462713051899d-1)*d  +.965735902808562554d-1)*d
     @    +1.38629436111989062d0
        return
2       ve=1d0
        vk=xlg
        return
        end




        subroutine dfringe(x,y,z,fact,a,fx,fy,fz)
c   calculates  entrance fringe field of dipole magnet
c   in the same manner as raytrace, for use in snake. jjl 2/17/87
c   n.b. can only be used for bending in the x-y plane.
c
c          fact = central field 
c             a = gap size (radius)
c         x,y,z = coordinates relative to magnet entrance
c                  (set entrance=0 in the relative ref frame)
c            fz = returned field value
c         c0(i) = fringe field coefficients
c     
c            fz = fact/(1+exp(s))
c             s = c0(0) + c0(1)*(y/a) + c0(2)*((y/a)**2) ....etc.
c****************************************************************
c
c  set up grid of b's for expansion (see raytrace manual p. 11-12)
        real x,y,z,fx,fy,fz,fact,a,delta
        real c0(0:5),s,b(-2:2,-2:2)
        data c0/0.2383,1.7395,-0.4768,0.5288,-0.1299,0.0222/
        data delta/80./
        do 20 j=-2,2
        do 20 k=-2,2
          if(abs(j)+abs(k).ge.3) go to 20
          s=c0(0)
          do 10 i=1,5

c         when curved efb's are introduced y will depend on x
10        s=s+(c0(i)*(((y+(j*delta))/a)**i))
20      b(j,k)=fact/(1+exp(s))

c       calculate fields
        fz=b(0,0)-((z**2/delta**2)
     @    * (((2./3.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     @    - ((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0)))))
     @    + ((z**4/delta**4)
     @    + (((-1./6.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     @    + ((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0))))
     @    + ((1./12.)*(b(1,1)+b(-1,1)+b(1,-1)+b(-1,-1)))
     @    - ((1./6.0)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(2.*b(0,0)))))))
        fy=((z/delta)
     @    * (((2./3.)*(b(1,0)-b(-1,0)))
     @    - ((1./12.)*(b(2,0)-b(-2,0)))))
     @    + (((z**3)/(delta**3))
     @    * (((1./6.)*(b(1,0)-b(-1,0)))
     @    - ((1./12.)*(b(2,0)-b(-2,0)))
     @    - ((1./12.)*(b(1,1)+b(1,-1)-b(-1,1)-b(-1,-1)
     @    - (2.*b(1,0))+(2.*b(-1,0))))))
        fx=((z/delta)
     @    * (((2./3.)*(b(0,1)-b(0,-1)))
     @    - ((1./12.)*(b(0,2)-b(0,-2)))))
     @    + (((z**3)/(delta**3))
     @    * (((1./6.)*(b(0,1)-b(0,-1)))
     @    - ((1./12.)*(b(0,2)-b(0,-2)))
     @    - ((1./12.)*(b(1,1)+b(-1,1)-b(1,-1)-b(-1,-1)
     @    - (2.*b(0,1))+(2.*b(0,-1))))))
        return
        end




        subroutine dipole(x,y,b0)
c*************************************************************************
c                                                                        *
c                       d i p o l e                                      *
c                                                                        *
c*************************************************************************
c       for meth=3 , type=2 , index=16+17 (arbitrary field numeric.deriv.)
c       field with gradient (alpha and beta)
c         (x=0.,y=0.,z=0.)=center of symmetry, z=0.:symmetry plane
c       this routine compute the field only in the symm. plane
c       the fringing field of the faces is computed
c       out of plane field needs numerical derivation.

c       rectangular-->cylindrical coord.:
        common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @              ,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @              ,tbound,dst
101     r=sqrt(x*x+y*y)
        rn=r/rms
        te=atan2(y,x)
        dr=rn-1.

c       z component of the field in sym.pl.:
        b0=1.-alps*dr+bets*dr*dr

c       input or exit fringing field:
        if(te.ge.tbound)then
          i=2
        else
          i=1
        endif
        ar0=abs(r0(i))
        sign=r0(i)/ar0
        er=e0(i)/b0
        rr=sqrt((x-xc(i))**2+(y-yc(i))**2)
        drf=(rr-ar0)*sign/er
        dr2=drf*drf
        sr=s0(i)+s1(i)*drf+s2(i)*dr2+s3(i)*dr2*drf
        ch=1./(1.+exp(sr))
        b0=b0*ch
1       return
        end




        subroutine geo(xyz,eject)
c*************************************************************************
c                                                                        *
c                       g e o                                            *
c                                                                        *
c*************************************************************************
c       eject=.true. if the point xyz is not in a free space.
c       coordinates are relative to the frame of the region indre.
        include 'sSnake.inc'
        parameter (dpi=2.*pi)
        parameter (ngalette=8)
        parameter (alphs2=pi/ngalette)
        parameter (alpha=2.*alphs2)
        parameter (r1=200.)
        parameter (r2=1200.)
        parameter (r3=(r2-r1)/2.)
        parameter (dr=50.)
        parameter (dz=10.)
        parameter (rmoy=(r1+r2)/2.)
        parameter (margr=50.)
        parameter (margz=50.)
        parameter (drt=dr+margr)
        parameter (dzt=dz+margz)
        parameter (rmaxt=r3+drt)
        parameter (rmaxt2=rmaxt**2)
        logical eject,first
        dimension xyz(3)
        character diag*100,diagt*100
        common/cdiag/diag(maxve),diagt
        common/ndiag/iepdead(maxve),iredead(maxve)
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ver(maxve,12),xyz0(maxre,3)
     @                ,frbox(2,maxre,3)
     @                ,zang(maxre),xang(maxre),yang(maxre)
     @                ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                ,indrer(maxre),tept(maxepr,maxre)
     @                ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                ,colli(maxepr,maxre,6)
        common/geoh/zcol,xcol,rcol,yaa,yaaa,zlim,zmin,ay2max,zmax
        common/geod/xtab(11),ytab(11),ztab(11),ttheta,zv,zh
     @              ,xa,ya,za,xb,yb,zb,xc,yc,zc,ya1,za1,yc1,zc1,yc2,zc2
        data first/.true./
        if(rname(indre).eq.'tor')goto2000
        eject=.false.
        return
2000    if(first)then
          ta=tan(alphs2)
          cas2=cos(alphs2)
          sas2=sin(alphs2)
          first=.false.
        endif
        eject=.false.
        ax=abs(xyz(1))
        y=xyz(2)
        az=abs(xyz(3))
        if(ax.eq.0.)goto10
        eject=az/ax.gt.ta
        if(eject)then

c         build the diag:
          write(diagt,'(2a,2(e10.3,a))')
     @          ' try to cross the'
     @          , ' plane of a coil ; angle ='
     @          ,az/ax,' max =',ta,' (rad)'
        return
        endif
10      continue
        z1=az*cas2-ax*sas2
        x1=ax*cas2+az*sas2
        if((x1+2.67*y+4200.).lt.0.)return
        if((x1-.49*y-620.).lt.0.)return
        if(abs(z1).gt.dzt)return
        eject=.true.
        az1=abs(z1)

c       built the diag:
        write(diagt,'(2a,2(e10.3,a))')
     @        'try to go inside the cryostat'
     @        , ' ; dist./coil plane ='
     @        ,az1,' min =',dzt
        return
20      x1=ax*cas2+az*sas2-rmoy
        ray2=y*y+x1*x1
        eject=ray2.lt.rmaxt2
        if(eject)write(diagt,'(2a,2(e10.3,a))')
     @                'try to go inside the cyl.'
     @                , ';dist./coil plane ='
     @                ,az1,' min =',dzt
        return
        end





c*********************************************************************
c                                                                    *
c         a p h a s a , r a n f , r a n g e t  and  r a n s e t      *
c                                                                    *
c*********************************************************************
        function aphasa(a1,a2)
        common /rand/ix
        aphasa=a1+ranf*(a2-a1)
        return
        end

        function ranf(ibidon)
        common /rand/ix
        equivalence (r,ir)
        iz=ix*899
        if(iz.le.0)then
          iz=iz+2147483647+1
        endif
        r=iz
        ranf=r/2147483647.  
        ix=iz
        return
        end

        subroutine ranget(irand)
        common /rand/ix
        irand=ix
        return
        end

        subroutine ranset(irand)
        common /rand/ix
        ix=irand
        return
        end

