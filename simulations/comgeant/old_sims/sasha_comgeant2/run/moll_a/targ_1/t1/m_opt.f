      SUBROUTINE M_OPT(E)
C
C===    Magnet optimization
C
C      E - beam energy
C
      IMPLICIT NONE
      REAL E
C
      VECTOR PAR_MAG(8,4),PAR_COL(5,7)
      VECTOR PZCO(8),PXCO(8)
C
      REAL p0,pcm,thetmx,ctmx,b,a1,a2,bmn(3),bmx(3),bst(3)
     +    ,bcur(3),ctcur,ctcur0,ctcur1,zco(8),xco(8),xsl(8),pl,pt
     +    ,hm(4),zm(4),zcol(7),xycol(2,2,7),zmd(4),rr,dr,dsl
      INTEGER i,j,i1,i2,i3,nstep,ithet,ip,miss,im
C
C     ------------------------------------------------------------------
C
      nstep=2
      bmn(1)=0.
      bmn(2)=60.0
      bmn(3)=0.
      bmx(1)=0.0001
      bmx(2)=60.0001
      bmx(3)=0.0001
      DO i=1,3
         bst(i)=(bmx(i)-bmn(i))/REAL(nstep)
      END DO
      CALL HBOOK2(101,'Thetamax vs Q1(X),Q2=Q3(Y)'
     +               ,nstep,0.,bmx(1),nstep,0.,bmx(2),0.)
C
C===    Dipole settings
C
      p0=E/2.
      pcm=SQRT(2.*0.0005*E)/2.
      a1=(PAR_COL(1,7)-PAR_MAG(1,4))/PAR_COL(5,7)
      a1=SQRT(1.+a1**2)
      a2=(PAR_COL(1,6)-PAR_MAG(1,4))/PAR_COL(4,6)
      a2=SQRT(1.+a2**2)
      b=-2.*p0/(a1+a2)/PAR_MAG(4,4)
      xmx=(a1-a2)/(a1+a2)
      thetmx=ACOS(xmx)
      WRITE(6,*) 'Dipole kGs=',b,xmx,thetmx*180./3.1415-90.
C
      DO i=1,4
         zm(i)=PAR_MAG(1,i)
         zmd(i)=PAR_MAG(3,i)
      END DO
      DO i=1,7
         zcol(i)=PAR_COL(1,i)
         xycol(1,1,i)=PAR_COL(2,i)
         xycol(2,1,i)=PAR_COL(3,i)
         xycol(1,2,i)=PAR_COL(4,i)
         xycol(2,2,i)=PAR_COL(5,i)
      END DO
C
C---   Loop on the fields in Q1,Q2 and Q3.
C---    1) assume Q2=Q3
C
      DO i=1,3
         bcur(i)=0.
      END DO
C      
      DO i1=1,nstep
         bcur(1)=bmn(1)+(i1-1)*bst(1)
C
         DO i2=1,nstep
            bcur(2)=bmn(2)+(i2-1)*bst(2)
C
C            DO i3=1,nstep
               i3=i2
               bcur(3)=bmn(3)+(i3-1)*bst(3)
C
               DO i=1,3
                  hm(i)=bcur(i)*PAR_MAG(4,i)
               END DO
C
C---            Find the max theta to pass
C
               ithet=0
               xcur1=2.
               xcur0=1.
 400           ithet=ithet+1
C               WRITE(6,*)  i1,i2,i3,ithet
               IF(ithet.EQ.1) then
                  xcur=0.
               ELSE
                  xcur=(xcur0+xcur1)/2.
               ENDIF
C
C---           Try two electrons
C
               miss=0
               DO ip=1,2
                  DO i=1,8
                     xco(i)=0.
                     xsl(i)=0.
                  END DO
                  pl=p0*(1.+xcur)
                  pt=pcm*SQRT(1.-xcur**2)
                  zco(1)=0.
                  xco(1)=0.
                  xsl(1)=pt/pl
C
                  DO im=1,3
                     zco(im+1)=zm(im)
                     xco(im+1)=xco(im)+xsl(im)*(zco(im+1)-zco(im))
                     rr=xco(im+1)
                     dsl=-hm(im)*rr/pl
                     xsl(im+1)=xsl(im)+dsl
C
C---                   Correction
C
                     dr=zmd(im)/4.*(xsl(im+1)-xsl(im))
                     WRITE(6,*) im,zmd(im),dsl,dr,rr,dr/rr
                     rr=xco(im+1)+dr
C                     xsl(im+1)=xsl(im)-hm(im)*rr/pl
C
                     IF(xco(im+1).LT.xycol(1,1,im).OR.
     +                  xco(im+1).GT.xycol(2,1,im)) THEN
                        miss=im+10*ip
                        GO TO 500
                     ENDIF
                  END DO
                  DO im=4,7
                     zco(im+1)=zcol(im)
                     xsl(im+1)=xsl(im)
                     xco(im+1)=xco(im)+xsl(im)*(zco(im+1)-zco(im))
                     IF(xco(im+1).LT.xycol(1,1,im).OR.
     +                  xco(im+1).GT.xycol(2,1,im)) THEN
                        miss=im+10*ip
                        GO TO 500
                     ENDIF
                  END DO
 500              CONTINUE
                  WRITE(6,2000) miss,i1,i2,i3,(bcur(i),i=1,3),ip,xcur
     +                         ,zco,xco,xsl
 2000             FORMAT(4I4,3F7.3,I4,F7.3
     +                  /5X,8F10.1/5X,8F10.2/5X,8F10.6)
                  DO im=1,8
                     PZCO(im)=zco(im)
                     PXCO(im)=xco(im)
                  END DO
               END DO
C
               IF(miss.NE.0) THEN
                  IF(xcur.LT.xcur0) xcur0=xcur 
               ELSE
                  IF(xcur.GT.xcur1) xcur1=xcur 
               ENDIF
C
C---            Continue the theta loop?
C
               IF(xcur0.GT.xcur1.AND.
     +            ithet.LT.50000.AND.
     +            ABS(xcur1-xcur0).GT.1.E-4) GO TO 400
C
               CALL HF2(101,bcur(1),bcur(2),xcur)
C
C            ENDDO
         ENDDO
      ENDDO
C
 999  CONTINUE
      END







