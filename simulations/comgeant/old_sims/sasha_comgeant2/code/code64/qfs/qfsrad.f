
      subroutine qfsrad(ebeam,p,theta,target,result,success)
****************************************************************
*     ebeam: beam energy in GeV
*     p: scattered particle momentum in GeV
*     th: scattered angle in sr
*     target: 0 proton, 1 deuterium
*     result: DXs in nb/Gev/sr
*     success: success flag
****************************************************************
      real ebeam,p,result,scale,PI,theta
      real E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      real sigqfza,sigda,sigxa,sigr1a,sigr2a,sig2na,sig,sigradl
      integer target,success
      COMMON/QFSPAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      
      SCALE=1.D-26
      PI=3.1415926
      success=0
      th=theta/PI*180.
      E=ebeam*1000.
      EPS=20.
      EPSD=20.
      PF=130.

      Z=1.
      if (target.eq.0) then
         A=1.
      else if (target.eq.1) then
         A=2.
      else
      endif
         
      W=(ebeam-p)*1000.
      IF(W.LE.1.) success=-1
      SIGQFZA=SIGQFS(E,TH,W,Z,A,EPS,PF)*SCALE
      write(6,*) 'sigqfs done'
      SIGDA=SIGDEL(E,TH,W,A,EPSD,PF)*SCALE
      write(6,*) 'sigdel done'
      SIGXA=SIGX(E,TH,W,A)*SCALE
      write(6,*) 'sigx done'
      SIGR1A=SIGR1(E,TH,W,A,PF)*SCALE
      write(6,*) 'sigr1 done'
      SIGR2A=SIGR2(E,TH,W,A,PF)*SCALE
      write(6,*) 'sigr2 done'
      SIG2NA=SIG2N(E,TH,W,Z,A,PF)*SCALE
      write(6,*) 'sig2n done'
      SIG=SIGQFZA+SIGDA+SIGXA+SIGR1A+SIGR2A+SIG2NA

      write(6,*) SIG
      CALL RADIATE(E,TH,W,SIG,SIGRADL)

      result=SIGRADL*1.D+36
     
      END
      
      



