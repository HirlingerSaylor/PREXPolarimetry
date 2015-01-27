* $Header:$
* $Log:$
*
#ifndef CERNLIB_CMAGPAR_INC
* #define CERNLIB_CMAGPAR_INC
C
C---           Magnets (general)
C
      INTEGER    MXMAGN,MXMAGNP
      PARAMETER (MXMAGN=12,MXMAGNP=12)
      COMMON/CMAGPAR/ NMAGN,MAGNPAR(MXMAGNP,MXMAGN)
C
      INTEGER         NMAGN,MAGNPAR
      REAL            AMAGNPAR(MXMAGNP,MXMAGN)
      EQUIVALENCE (MAGNPAR(1,1),AMAGNPAR(1,1))
C
#if defined CERNLIB_COMMENT
C
C===   NMAGN  - number of defined magnets , m=1,NMAGN
C=F=   MAGNPAR( 1,m)   - type =0 - map (not used yet), =1 -quad, =2 - dipole (Bx=B), =-2 dipole (By=B)
C                             =3 - solenoid, 4 - solenoid(map), 
C                             =5 - general map (from snake) - MAD magnets
C             ( 2,m)   - mother volume (say, 'HALL'), ='    ' - magnet not used
C     AMAGNPAR( 3-5,m) - X,Y,Z of the center in DRS of the mother volume
C      MAGNPAR( 6,m)   - irot  - rotation matrix number, =0 - no rotation
C     AMAGNPAR( 7,m)   - X - 1/2 size (Q - radius) of the field   
C             ( 8,m)   - Y - 1/2 size (Q - radius)   
C             ( 9,m)   - Z - 1/2 size 
C             (10,m)   - B (field at the center) 
C      MAGNPAR(11-12,m)   - auxil. flags/parameters
#endif
C
#endif
