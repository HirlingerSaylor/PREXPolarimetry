

c
c       indfi=5:3 dim table
c
c
805	continue
c index 6: same initialization as index 5:
806     read(3,*)ndx,xmin,xmax
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
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in x'
        read(3,*)ndy,ymin,ymax
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
c 		if cylindrical , convert degree to radian :
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
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
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
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in z'
c   read the format:
       		 read(3,'(a)')form
        if(.not.lmap)then
		close(3)
		return
	endif
c   read the field
       		 if(form(1:3).eq.'bin')then
       		     read(3)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3)
     @                   ,i=1,ndx),j=1,ndy),k=1,ndz)
       		 else
        	    read(3,form)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3)
     @                   ,i=1,ndx),j=1,ndy),k=1,ndz)
       		 endif
	goto100


