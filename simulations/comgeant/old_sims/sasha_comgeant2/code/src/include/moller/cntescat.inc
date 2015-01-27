* $Header:$
* $Log:$
*
#ifndef CERNLIB_CNTESCAT_INC
* #define CERNLIB_CNTESCAT_INC
C
C---           e- scattering ntuple
C
      COMMON/CNTESCAT/ NEWOR,IETYPE,THETCM,PHICM,WGCROS,ANPOWER
C
      INTEGER          NEWOR    ! number of words
     +                ,IETYPE   ! the type of the process:
C                                 =1 - Moller 
      REAL             THETCM   ! theta in CM of the first particle in the process 
     +                ,PHICM    ! phi   in CM of the first particle in the process 
     +                ,WGCROS   ! event weight (=1. if the generator takes it into 
C                                                account)
     +                ,ANPOWER  ! the analyzing power (for Moller scattering)  
C
#if defined CERNLIB_COMMENT
C
#endif
C
#endif
