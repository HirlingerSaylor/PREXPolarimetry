* $Header:$
* $Log:$
*
#ifndef CERNLIB_OMCUNIT_INC
* #define CERNLIB_OMCUNIT_INC
C
C---           e- scattering flags
C
      COMMON/CESCATT/ THETMOL(2),PHIMOL(2),IMSHELL
     +               ,IDMOLTET
     +               ,THETMOT(2),PHIMOT(2)
     +               ,IMOTTWG,MOTTARG,IDMOTT,ALPHMOTT
     +               ,IRADMOT,RADMOTL(2),IDRADMOT
     +               ,FERMIMOM
     +               ,IPHOTPR,THETPHPR(2),PMOMPHPR(2),IPHPRGEA(6)
     +               ,IDPHOTPR
     +               ,NTWRHIT
     +               ,IPHOTVMP,IPHOTVMT,PHOTVMTT
     +               ,IPARONE,PAVONE(3),DPAVONE(3)
C
      INTEGER         IMSHELL
     +               ,IDMOLTET
     +               ,IMOTTWG,MOTTARG,IDMOTT,IRADMOT,IDRADMOT
     +               ,IPHOTPR,IPHPRGEA,IDPHOTPR
     +               ,NTWRHIT
     +               ,IPHOTVMP,IPHOTVMT
     +               ,IPARONE
C
      REAL            THETMOL,PHIMOL,THETMOT,PHIMOT,ALPHMOTT,RADMOTL
     +               ,FERMIMOM
     +               ,THETPHPR,PMOMPHPR
     +               ,PHOTVMTT
     +               ,PAVONE,DPAVONE
C
#if defined CERNLIB_COMMENT
C
C=F=   IMSHELL         = 0 - scatterinmg on an external shell
C                        1 - L - shell
C                        2 - K - shell
C---   IDMOLTET          ID for a histogram of Theta CM dependence
C=F=   THETMOL(1-2)    lower/upper limits for CM theta of Moller scattering, degree
C=F=   PHIMOL(1-2)     lower/upper limits for CM phi of Moller scattering  , degree
C
C=F=   THETMOT(1-2)    lower/upper limits for LAB theta of Mott scattering , degree
C=F=   PHIMOT(1-2)     lower/upper limits for LAB phi of Mott scattering   , degree
C=F=   MOTTARG         - the target particle number for Mott
C                            A>0 a spinless nucleus with a certain radius is considered
C
C=F=   IMOTTWG          = 0 - simulate events without a weight.
C                           The Mott cross-section is approximated with an exponential 
C                           A substitution theta=-log(x)/alpha, or x=exp(-alpha*theta) is made.
C                           The Mott cross-section on the variable x is stored in the histograms.
C                           Events are simulated using HRNDM1   
C                       = 1 - "light" weights: as before, but take a even distrubution of x
C                           and assign a weight to the event (cross-section on x) 
C                       = 1 - "heavy" weights: even distribution in theta and the cross-section
C                              used for weight. These weights differ by many orders of magnitude.
C      IDMOTT       - ID for a histogram of x(theta) LAB dependence in the 
C    ALPHMOTT       - the alpha slope of the exponential
C
C=F= IRADMOT        >0 - Radiative Mott. It uses the Levchuk's code instead of the regular Mott code.
C                        works only at IMOTTWG=0
C
C=F= FERMIMOM        (default=0.045GeV) Fermi max momentum. Fermi motion is taken into account 
C                                      for the Mott scattering and vector meson photoproduction
C
C=F= RADMOTL(1:2)   -    Limits on the energy fraction retained by the electron                        
C    IDRADMOT       - ID for a histogram of Rad. loss. dependence (for Rad. Mott) 
C
C=F=   IPHOTPR      =n - photoproduction in the target (Bosted and Wiser's program)
C                          n=1,2,3,4,5,6 - pi+,pi-,K+,K-,p,pbar
C=F=   THETPHPR     - theta (LAB) limits for "pions"
C=F=   PMOMPHPR     - momentum/Ebeam (LAB) limits for "pions"
C=F=   IPHPRGEA     - corresponding GEANT numbers of pions,kaons.. 
C=F=   IDPHOTPR     - ID of a histogram for photoproduction
C
C----              Vector meson photoproduction
C=F=   IPHOTVMP      - the vector meson particle type 
C=F=   IPHOTVMT      - the target particle type
C                         (say Atarg=1: A(target material)=12 - Fermi motion is taken into account)
C=F=   PHOTVMTT     - <t> - slope for the t- dependence: dsigma/dt = sigma0*exp(t/PHOTVMTT)
C
C=F=   NTWRHIT      =0 - write all events to NTUPLES
C                   >0 - write only events with hits
C
C----      One particle for spectrometer tests
C=F==  IPARONE       - particle type 
C=F==  PAVONE(1:3)   - central P, theta, phi (GeV, rad, rad)
C=F==  DPAVONE(1:3)  - +/- range on P,thet,phi, say range in thet is PAVONE(1)-DPAVPONE(1):PAVONE(1)+DPAVPONE(1)
C
#endif
C
#endif
