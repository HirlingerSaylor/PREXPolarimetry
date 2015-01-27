      REAL FUNCTION BOSTED(XX)
C
C      REAL XX
C      XX=BOSTED(20.)
C      WRITE(6,*) 'XX=',XX
C      END
C ! Calculate phtoproduction singles rates from Wiser fit for Moller

      IMPLICIT NONE
      REAL XX,YY
      REAL E0,RADLEN,P_PI,THETA_DEG,TH,DOMEGA,DEP,SIGMA,TARG
      REAL RATE,FLUX,SUM,SUMP(100),SUMT(100),PT,SUMM
      INTEGER IPP,ITH
      REAL    tfac
      INTEGER i

      tfac=1.
      tfac=XX
      RADLEN=12.5 ! effective rad. length in percent for target
C                 ! (half of full target)
      FLUX=6.E11 ! beam in 1 pulse
      TARG=10. * 0.6 ! 1.5 m LH2 (in barns)
      E0=48.
C
      RADLEN=RADLEN*tfac
      TARG=TARG*tfac
C      open(unit=6,file='wise.out')
C      open(unit=22,file='wise.plt1')


C ! get yield at angle th and integrate up to 12 mr
      DO i=1,100
         SUMP(i)=0.
         SUMT(i)=0.
      ENDDO
      SUM=0.
      SUMM=0.
      ITH=1
      DO ITH=1,12
         IPP=1
        DO IPP=10,40
          P_PI = 1.0 * float(ipp)
C          P_PI = XX
          TH = 0.0005 * FLOAT(2*ITH-1)
C          TH=YY
          DOMEGA = 2. * 3.14 * TH * 0.001 
          THETA_DEG = TH * 57.295
          DEP = 1.0
          PT= P_PI * TH

          CALL WISER_ALL_SIG (E0,P_PI,THETA_DEG,RADLEN,5,SIGMA)
          RATE = SIGMA * DEP * DOMEGA *
     >         1.E-9 * FLUX * TARG 
          SUMP(IPP)=SUMP(IPP)+RATE
          SUMT(ITH)=SUMT(ITH)+RATE
          SUM = SUM + RATE
          IF(PT.GT.0.08.AND.PT.LT.0.12) SUMM = SUMM + RATE
        ENDDO
      ENDDO
C      WRITE(6,'(1X,''P,PI/SPILL='',F5.1,F10.0)') 
C     >   (FLOAT(IPP),SUMP(IPP),IPP=1,40)
C      WRITE(6,'(1X,''TH,PI/SPILL='',F8.4,F10.0)') 
C     >   (0.0005*FLOAT(2*ITH-1),SUMT(ITH),ITH=1,17)
C      WRITE(6,'(1X,''TOTAL PI/SPILL='',F10.0)') SUM
C      WRITE(6,'(1X,''TOTAL PI/SPILL (0.08<pt<0.12)='',F10.0)') SUMM
      BOSTED=SUM
C      WRITE(6,*) 'Fun=',BOSTED
C
 99   CONTINUE
      END

                                                    

      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   1717.   
      DOUBLE PRECISION FUNCT,LOWER,UPPER,EPSLON                         1718.   
         INTEGER NLVL                                                   1719.   
      INTEGER   LEVEL,MINLVL,MAXLVL,RETRN(50),I                         1720.   
      DATA      MINLVL/3/,MAXLVL/24/
      DOUBLE PRECISION VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),  1721.   
     1   FMRX(50), ESTRX(50), EPSX(50)                                  1722.   
      DOUBLE PRECISION  R, FL, FML, FM, FMR, FR, EST, ESTL
     +      , ESTR, ESTINT,L,                                           1723.   
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         1724.   
         LEVEL = 0                                                      1725.   
         NLVL = 0                                                       1726.   
         ABAREA = 0.0                                                   1727.   
         L = LOWER                                                      1728.   
         R = UPPER                                                      1729.   
         FL = FUNCT(L)                                                  1730.   
         FM = FUNCT(0.5*(L+R))                                          1731.   
         FR = FUNCT(R)                                                  1732.   
         EST = 0.0                                                      1733.   
         EPS = EPSLON                                                   1734.   
  100 LEVEL = LEVEL+1                                                   1735.   
      M = 0.5*(L+R)                                                     1736.   
      COEF = R-L                                                        1737.   
      IF(COEF.NE.0) GO TO 150                                           1738.   
         ROMBRG = EST                                                   1739.   
         GO TO 300                                                      1740.   
  150 FML = FUNCT(0.5*(L+M))                                            1741.   
      FMR = FUNCT(0.5*(M+R))                                            1742.   
      ESTL = (FL+4.0*FML+FM)*COEF                                       1743.   
      ESTR = (FM+4.0*FMR+FR)*COEF                                       1744.   
      ESTINT = ESTL+ESTR                                                1745.   
      AREA=DABS(ESTL)+DABS(ESTR)                                        1746.   
      ABAREA=AREA+ABAREA-DABS(EST)                                      1747.   
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     1748.   
         NLVL = NLVL+1                                                  1749.   
         ROMBRG = ESTINT                                                1750.   
         GO TO 300                                                      1751.   
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         1752.   
     1         (LEVEL.LT.MINLVL))  GO TO 400                            1753.   
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             1754.   
  300    LEVEL = LEVEL-1                                                1755.   
         I = RETRN(LEVEL)                                               1756.   
         VALINT(LEVEL, I) = ROMBRG                                      1757.   
         GO TO (500, 600), I                                            1758.   
  400    RETRN(LEVEL) = 1                                               1759.   
         MX(LEVEL) = M                                                  1760.   
         RX(LEVEL) = R                                                  1761.   
         FMX(LEVEL) = FM                                                1762.   
         FMRX(LEVEL) = FMR                                              1763.   
         FRX(LEVEL) = FR                                                1764.   
         ESTRX(LEVEL) = ESTR                                            1765.   
         EPSX(LEVEL) = EPS                                              1766.   
         EPS = EPS/1.4                                                  1767.   
         R = M                                                          1768.   
         FR = FM                                                        1769.   
         FM = FML                                                       1770.   
         EST = ESTL                                                     1771.   
         GO TO 100                                                      1772.   
  500    RETRN(LEVEL) = 2                                               1773.   
         L = MX(LEVEL)                                                  1774.   
         R = RX(LEVEL)                                                  1775.   
         FL = FMX(LEVEL)                                                1776.   
         FM = FMRX(LEVEL)                                               1777.   
         FR = FRX(LEVEL)                                                1778.   
         EST = ESTRX(LEVEL)                                             1779.   
         EPS = EPSX(LEVEL)                                              1780.   
         GO TO 100                                                      1781.   
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          1782.   
      IF(LEVEL.GT.1) GO TO 300                                          1783.   
      QUADMO = ROMBRG /12.0D0                                           1784.   
      RETURN                                                            1785.   
      END                                                               1786.   

      Subroutine WISER_ALL_SIG(E0,P,THETA_DEG,RAD_LEN,TYPE,SIGMA)

C!------------------------------------------------------------------------------
C! Calculate pi,K,p  cross section for electron beam on a proton target
C! IntegrateQs over function WISER_FIT using integration routine QUADMO
C! E0         is electron beam energy, OR max of Brem spectra
C! P,E       is scattered particle  momentum,energy
C! THETA_DEG  is kaon angle in degrees
C! RAD_LEN (%)is the radiation length of target, including internal
C!                (typically 5%)
C!               = .5 *(target radiation length in %) +5.
C!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***
C! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
C! SIGMA      is output cross section in nanobars/GeV-str
C!------------------------------------------------------------------------------

      IMPLICIT NONE       
      REAL E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE ,M_X,U_MAN
      REAL E,P_COM,COST,P_T,M_X,U_MAN
      INTEGER TYPE_COM,PARTICLE
C!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      REAL A5(6),A6(6),MASS2(3),MASS(3),MP,MP2,RADDEG
      DATA A5/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      DATA A6/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      DATA MASS2/.019488, .2437, .8804/
      DATA MASS/.1396, .4973, .9383/ 
      DATA MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      REAL  M_L,E_GAMMA,SIG_E
      DOUBLE PRECISION E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08
     +                ,EPSILON
      DATA            EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE
                                            
      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
C!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
C!     >     E_GAMMA_MIN,P_T
C!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION WISER_ALL_FIT(E_GAMMA)

C!---------------------------------------------------------
C! Calculates  pi, k, p  cross section for gamma + p -> k
C!  It is already divided by E_GAMMA, the bremstrulung spectra
C! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
C! Can be called from WISER_SIG using integration routine QUADMO
C! E,P are KAON energy and momentum
C! P_t is KAON transverse momentum
C! P_CM is KAON center of mass momentum
C! P_CM_L is KAON center of mass longitudinal momentum
C! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
C! E_GAMMA is photon energy.
C!             Steve Rock 2/21/96
C!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      REAL  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
C!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      REAL A1(6),A2(6),A3(6),A4(6),A6,A7,MASS2(3),MASS(3),MP2,MP,RADDEG
      DATA A1/566.,  486.,   368., 18.2,  1.33E5,  1.63E3 / 
      DATA A2/829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 / 
      DATA A3/1.79,  1.77,   1.91, 0.98,  1.41,    1.79 / 
      DATA A4/2.10,  2.18,   1.15, 1.83,   .72,    2.24 /
      DATA A6/1.90/,A7/-.0117/ !proton only
      DATA MASS2/.019488, .2437, .8804/
      DATA MASS/.1396, .4973, .9383/ 
      DATA MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      REAL X_R,  SIGMA, SIG_E,S,B_CM, GAM_CM,  P_CM
      REAL P_CM_MAX, P_CM_L,PART
      DOUBLE PRECISION E_GAMMA
                                            

C!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

C!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END
