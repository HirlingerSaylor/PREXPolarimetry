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
      END
