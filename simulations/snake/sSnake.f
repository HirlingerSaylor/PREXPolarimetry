c*************************************************************************
c                                                                        *
c                       s n a k e                                        *
c                                                                        *
c*************************************************************************
c
c        general ray tracing program for multi regions problems
c                       p.vernin 21/07/86
c                version 2 16/07/88
c                version 3 31/05/91
c
c       units:mm,gev/c,t
c       input files:
c             -geometry and field are defined region after region in
c                the "description file"
c             -the description file may itself refer to "field files"
c       output file:
c             -vectors for mudifi (cern program for multi dimensional fit)
c                may be stored on a file
c       parameters:
c             maxve:maximum number of input vector(s) that can be defined
c             maxre:   "      "    "  region(s)       "
c             maxep:   "      "    "  end-plane(s)    "
c             maxepr:  "      "    "  end-plane(s) by region
c             maxplv:  "      "    "  trajectories   that can be plotted
c             maxplp:  "      "    "  point per traj.    "
c             maxq:    "      "    "  tic marks          "
c             maxdat: "      "    "  additional data to describe a region
c       common /absfra/:
c             nuve:number of input vector(s) defined
c             nuep:number of end-plane(s) defined
c             nure:number of region(s)  defined
c             vea(iv,iep,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
c                iep<=nuep<=maxep:end-plane index
c                ic=1:x coord. in absolute frame
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum in absol.frame
c                   5:y            "
c                   6:z            "
c                   7:p module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c            10 to 15:the same as 1 to 6 but in relative frame coordinates
c		   16:x component of the spin
c		   17:y          "
c		   18:z          "
c	     19 to 21: the same as 16 to 18 but in relative frame coord.
c       common /relfra/:
c             indre=current region   index
c             indep=   "    end-plane "
c             ilive=   "    number of alive vector(s)
c             ver(iv,ic):description of the particules in the current region
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord.relative to the current region(in relative frame)
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum(in relat.frame)
c                   5:y            "
c                   6:z            "
c                   7:p:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c		   10:x component of the spin (in relat.frame)
c		   11:y          "
c		   12:z          "
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             frbox(maxi or mini,region index,x y or z):relat.coord.
c                of the free-box
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c             adata(data index < maxdat,region index):additional data to
c                describe a region
c       common /fidef/:
c             method(region index):methode to use to carry particles
c                in this region
c             itype(region index):type in this method
c             indfi(region index):index in this type and in this method
c             fact(region index):multipl.factor of the field or ref.mommentum
c             fifi(region index):field file name
c       common /plot/:
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(submit index,traj.index,x y z bx by bz...
c                    ...px py pz):absol.coord.of submit,
c				  absol. field components and
c				  absol. spin components at this point
c             nq<=maxq:number of tic marks to be plotted
c             h1(submit index,tic mark index):absol.horiz.coord of the submit
c                of the tic mark
c             v1(submit index,tic mark index):absol.vertic.coord of the submit
c                of the tic mark
c             ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c             ih:       "       "           "           "      horizontal "
c             iv:       "       "           "           "      vertical   "
c             kp= 1:observer is at ip>0.
c                -1:   "     "  "  IP<0.
c             plobox(submit index,horiz. or vertic.,region index):abs.coord
c                of submit of free-box(special submit list for 3d projection)
				program snake

				character*20 fide

				include 'sSnake.inc'
				character diag*100,diagt*100
				common/cdiag/diag(maxve),diagt
				common/ndiag/iepdead(maxve),iredead(maxve)
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
				common/regin/veri(maxve,maxre,12)
				character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
				common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @								,title(maxre),cerfiln(maxepr,maxre)
				common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @								,frbox(2,maxre,3)
     @								,zang(maxre),xang(maxre),yang(maxre)
     @								,indve,indre,indep,ilive,adata(maxdat,maxre)
     @								,indrer(maxre),tept(maxepr,maxre)
     @								,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @								,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @								,colli(maxepr,maxre,6)
				common /comlight/cl(3)
				character*20 fidi,fibl,fidit,fide,fimu,fimut,firay,outData,outPts
				character*80 line
				character*1 rep,rip
				character rep2*2
				character*24 daet 
				logical cyl,lfil,lmap,ladd
				character*80 fifi
				common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @					,cyl(maxre),lfil
				common /cfidef/ fifi(maxre)
				common/alter/jflag,prec,hstmax,hstmin,sinc
				common ntraj
				dimension xmax(8),xmin(8),plottn(2,maxplv,6)


				data fide/'d.dat'/
				data outData/'data.txt'/
				data outPts/'pts.txt'/


				data fibl/'                    '/
				data rip/'r'/
				data fimu/'mud.mud'/
				data firay/'temp.dat'/
				data cl/1.,1.,1./
				lfil=.false.
				name(1)='x'
				name(2)='y'
				name(3)='z'
				name(4)='bx'
				name(5)='by'
				name(6)='bz'
				name(7)='sx'
				name(8)='sy'
				name(9)='sz'
				name(10)='l'
				name(11)='b'
				name(12)='r'
				name(13)='t'
				name(14)='rt'
				name(15)='s'
				unit(1)='mm'
				unit(2)='mm'
				unit(3)='mm'
				unit(4)='T'
				unit(5)='T'
				unit(6)='T'
				unit(7)=' '
				unit(8)=' '
				unit(9)=' '
				unit(10)='mm'
				unit(11)='T'
				unit(12)='mm'
				unit(13)='deg'
				unit(14)='mm'
				unit(15)=' '
				fidi=fide
				nplv=0
				nplvold=0





c				write(6,'(/,''$ directive-file name(def='',a,'')'')')fide
				open(2,file=fide,err=428,status='old')
				goto429

428			continue
				write(6,'('' can''''t open directive file'')')
				stop

c				read directives (descin reads 2):
429			call descin
				close(2)

c				descin leaves indre as max region index (starts at 1)
				nure=indre

c       fill the absolute coord. buff.,set region index to 0:
				call spring

c				set end-plane index to 0:
				indep=0
				indved=0
				nulu=0
				luser=.false.

c**************************new region:*************************

				do 200 indrel=1,nure
					indre=indrel

c					title from descin, can't find ilive
					write(6,*)title(indre)
					if(ilive.le.0)then
						write(6,'(''  alive:'',i4,'' vector, region skipped'')')
     @						ilive
						goto 200
					endif

c       fill the relat. coodr. buff. with the abs. one(last e-p index)
c        accordind to new geometrical specif.
					call rot

c				initialize the field routine according to new field specif.
					lmap=.true.
					call ifield(lmap,ladd)
					if(ladd)call add_field

c       loop on the end-planes of this region:
					do 1 iep=1,nep(indre)
						if(epsw(indre).eq.'loop')then
							kep=2			
							yep=yepmin(indre)+yepstp(indre)*(iep-1)
							tep=0.
						else
							kep=kept(iep,indre)
							yep=yept(iep,indre)
							tep=tept(iep,indre)
						endif

c				send the particle to the next end-plane...
c        ...by working on the relat.coord.buffer
						call send(kep,yep,tep,iep)

c				increment the end-plane index and fill the abs.coo.buff.:
						call unrot
						write(6,'(''  alive:'',i4,'' vector(s)'')')ilive
1					continue
200			continue

c*************************end of the ray-tracing****************


c       save the number of end-planes:
				nuep=indep
c				unclip the plot:
				call unclip
600			continue
				goto 530



c530			write(6,'(''$Number of trajectories (def='',i4,'')'')')ntraj
530			open(19,file=outData)
				open(20,file=outPts)
c				read(5,'(i3)')indve  
				if(nuve.le.0)then
					write(6,*)' snake:nothing to write'
					goto600
				endif
				do l=1,ntraj
					do k=1,nplp(l)
						write(19,*)l,k,vea(l,k,1),vea(l,k,2),vea(l,k,3)
					enddo
				enddo
				do l=1,ntraj
					do k=1,nplp(l)
						write(20,*)l,k,plott(k,l,1),plott(k,l,2),plott(k,l,3)
					enddo
				enddo
				write(6,*)
				stop


				end
