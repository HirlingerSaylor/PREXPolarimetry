
      SUBROUTINE ROM(A,B,EPS,ANS,K)
      IMPLICIT REAL (A-H,O-Z)
C  ROMBERG METHOD OF INTEGRATION
      DIMENSION W(50,50)
      H=B-A
      K=0
      CALL VALY(A,FA)
      CALL VALY(B,FB)
      W(1,1)=(FA+FB)*H/2.
    4 K=K+1
      IF(K.GE.49)GO TO 5
      H=H/2.
      SIG=0.
      M=2**(K-1)
      write(6,*) 'VALY loop ',M,EPS,A,B,FA,FB
      DO J=1,M
         J1=2*J-1
         X=A+FLOAT(J1)*H
         CALL VALY(X,F)
         IF(K.EQ.1) write(6,*) ' ROM loop j,x,f ',J,X,F
C         stop
         SIG=SIG+F
      ENDDO
      write(6,*) 'ROM SIG=',SIG
      stop
      W(K+1,1)=W(K,1)/2.+H*SIG
      DO L=1,K
         IU=K+1-L
         IV=L+1
         W(IU,IV)=(4.**(IV-1)*W(IU+1,IV-1)-W(IU,IV-1))/(4.**(IV-1)-1.)
      ENDDO
      E=(W(IU,IV)-W(IU,IV-1))/W(IU,IV)
      write(6,*) 'ROM loop, K,E,EPS=',K,E,EPS
C
      IF(ABS(E)-EPS.GT.0.) GO TO 4
C
      ANS=W(1,IV)
      write(6,*) 'ROM end, K=',K
      RETURN
    5 PRINT 100
      write(6,*) 'ROM error'
  100 FORMAT(' K OVERFLOW')
      CALL EXIT
      END
