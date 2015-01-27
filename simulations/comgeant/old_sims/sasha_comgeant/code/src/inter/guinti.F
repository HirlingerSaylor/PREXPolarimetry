      SUBROUTINE GUINTI
      CHARACTER*1 BS1
      PARAMETER (BS1='\\')
      PARAMETER (MGUIDL=199)
      CHARACTER*80 GUID
      COMMON /KCGUID/ GUID(MGUIDL)
      EXTERNAL OMXDRAW
      EXTERNAL OMXCNTR
      EXTERNAL OUXCNTR
      EXTERNAL OMXPLOT

      CALL KUCMD(' ','OMGEANT','C')
      GUID(1)='Commands for OMGEANT'
      CALL KUGUID('OMGEANT',GUID,1,'S')

      CALL KUCMD('OMGEANT',' ','SW')

      CALL KUCMD(' ','DRAW','C')
      GUID(1)='Drawing commands'
      CALL KUGUID('DRAW',GUID,1,'S')

      CALL KUCMD('DRAW',' ','SW')

      CALL KUCMD(' ','DDIGI','C')
      CALL KUPAR('DDIGI','IUSET','Name of the detector SET','CO','S')
      CALL KUPVAL('DDIGI','IUSET',0,0.,'*','D')
      CALL KUPAR('DDIGI','IUDET','Name of the detector','CO','S')
      CALL KUPVAL('DDIGI','IUDET',0,0.,'*','D')
      CALL KUPAR('DDIGI','IUNI','Number of the unit','IO','S')
      CALL KUPVAL('DDIGI','IUNI',0,0.,' ','D')
      GUID(1)='Plots the digitisations in wire chambers (DDIGI [set] [de
     +t] [unit])'
      CALL KUGUID('DDIGI',GUID,1,'S')
      CALL KUACT('DDIGI',OMXDRAW)

      CALL KUCMD(' ','SPARATT','C')
      CALL KUPAR('SPARATT','IPAR','Particle number','IO','S')
      CALL KUPVAL('SPARATT','IPAR',1,0.,' ','L')
      CALL KUPVAL('SPARATT','IPAR',8,0.,' ','H')
      CALL KUPAR('SPARATT','MODE',
     +'Mode of the line (=0 - skip these particles)','IO','S')
      CALL KUPAR('SPARATT','COLOR','Color of the line','IO','S')
      GUID(1)='Set the line mode and color for drawing trajectories of d
     +ifferent particles'
      GUID(2)='          ITRTYP  PARTICLE  | mode | color |       HIGZ '
     +//'        |'
      GUID(3)='          -----------------------------------------------
     +--------'
      GUID(4)='              1   gammas    |  3   |  4    | dotted     '
     +//' blue   |'
      GUID(5)='              2   electrons |  1   |  2    | solid      '
     +//' red    |'
      GUID(6)='              3   neutral   |  4   |  1    | dot-dashed '
     +//' black  |'
      GUID(7)='              4   hadrons   |  1   |  2    | solid      '
     +//' red    |'
      GUID(8)='              5   muons     |  2   |  3    | dashed     '
     +//' green  |'
      GUID(9)='              6   geantino  |  1   |  1    | dot-dashed '
     +//' black  |'
      GUID(10)='              7   cerenkov  |  3   |  6    | dotted    '
     +//' magenta |'
      CALL KUGUID('SPARATT',GUID,10,'S')
      CALL KUACT('SPARATT',OMXCNTR)

      CALL KUCMD(' ','ITX','C')
      CALL KUPAR('ITX','X','X - coordinate','R','S')
      CALL KUPVAL('ITX','X',0,10.,' ','D')
      CALL KUPAR('ITX','Y','Y - coordinate','R','S')
      CALL KUPVAL('ITX','Y',0,10.,' ','D')
      CALL KUPAR('ITX','TEXT','Text to be drawn','C','S')
      CALL KUPVAL('ITX','TEXT',0,0.,' ','D')
      GUID(1)='Draw a text (ITX X Y text)'
      CALL KUGUID('ITX',GUID,1,'S')
      CALL KUACT('ITX',OMXDRAW)

      CALL KUCMD('..',' ','SW')

      CALL KUCMD(' ','CONTROL','C')
      GUID(1)='Control flags and commands'
      CALL KUGUID('CONTROL',GUID,1,'S')

      CALL KUCMD('CONTROL',' ','SW')

      CALL KUCMD(' ','TKINE1','C')
      CALL KUPAR('TKINE1','tkin1','Type of simulated particle','I','S')
      CALL KUPVAL('TKINE1','tkin1',8,0.,' ','D')
      CALL KUPAR('TKINE1','tkin2',
     +'Number of simulated particles of this type','I','S')
      CALL KUPVAL('TKINE1','tkin2',1,0.,' ','D')
      CALL KUPAR('TKINE1','tkin3','P minimum /GeV','R','S')
      CALL KUPVAL('TKINE1','tkin3',0,1.,' ','D')
      CALL KUPAR('TKINE1','tkin4','P maximum /GeV','R','S')
      CALL KUPVAL('TKINE1','tkin4',0,100.,' ','D')
      CALL KUPAR('TKINE1','tkin5','Pt mean /GeV','R','S')
      CALL KUPVAL('TKINE1','tkin5',0,0.,' ','D')
      GUID(1)=' Parameters of 1-st class of particles to simulate'
      CALL KUGUID('TKINE1',GUID,1,'S')
      CALL KUACT('TKINE1',OMXCNTR)

      CALL KUCMD(' ','TKINE2','C')
      CALL KUPAR('TKINE2','tkin1','Type of simulated particle','I','S')
      CALL KUPVAL('TKINE2','tkin1',8,0.,' ','D')
      CALL KUPAR('TKINE2','tkin2',
     +'Number of simulated particles of this type','I','S')
      CALL KUPVAL('TKINE2','tkin2',1,0.,' ','D')
      CALL KUPAR('TKINE2','tkin3','P minimum /GeV','R','S')
      CALL KUPVAL('TKINE2','tkin3',0,1.,' ','D')
      CALL KUPAR('TKINE2','tkin4','P maximum /GeV','R','S')
      CALL KUPVAL('TKINE2','tkin4',0,100.,' ','D')
      CALL KUPAR('TKINE2','tkin5','Pt mean /GeV','R','S')
      CALL KUPVAL('TKINE2','tkin5',0,0.,' ','D')
      GUID(1)=' Parameters of 2-nd class of particles to simulate'
      CALL KUGUID('TKINE2',GUID,1,'S')
      CALL KUACT('TKINE2',OMXCNTR)

      CALL KUCMD(' ','TKINE3','C')
      CALL KUPAR('TKINE3','tkin1','Type of simulated particle','I','S')
      CALL KUPVAL('TKINE3','tkin1',8,0.,' ','D')
      CALL KUPAR('TKINE3','tkin2',
     +'Number of simulated particles of this type','I','S')
      CALL KUPVAL('TKINE3','tkin2',1,0.,' ','D')
      CALL KUPAR('TKINE3','tkin3','P minimum /GeV','R','S')
      CALL KUPVAL('TKINE3','tkin3',0,1.,' ','D')
      CALL KUPAR('TKINE3','tkin4','P maximum /GeV','R','S')
      CALL KUPVAL('TKINE3','tkin4',0,100.,' ','D')
      CALL KUPAR('TKINE3','tkin5','Pt mean /GeV','R','S')
      CALL KUPVAL('TKINE3','tkin5',0,0.,' ','D')
      GUID(1)=' Parameters of 3-rd class of particles to simulate'
      CALL KUGUID('TKINE3',GUID,1,'S')
      CALL KUACT('TKINE3',OMXCNTR)

      CALL KUCMD(' ','TKINE4','C')
      CALL KUPAR('TKINE4','tkin1','Type of simulated particle','I','S')
      CALL KUPVAL('TKINE4','tkin1',8,0.,' ','D')
      CALL KUPAR('TKINE4','tkin2',
     +'Number of simulated particles of this type','I','S')
      CALL KUPVAL('TKINE4','tkin2',1,0.,' ','D')
      CALL KUPAR('TKINE4','tkin3','P minimum /GeV','R','S')
      CALL KUPVAL('TKINE4','tkin3',0,1.,' ','D')
      CALL KUPAR('TKINE4','tkin4','P maximum /GeV','R','S')
      CALL KUPVAL('TKINE4','tkin4',0,100.,' ','D')
      CALL KUPAR('TKINE4','tkin5','Pt mean /GeV','R','S')
      CALL KUPVAL('TKINE4','tkin5',0,0.,' ','D')
      GUID(1)=' Parameters of 4-th class of particles to simulate'
      CALL KUGUID('TKINE4',GUID,1,'S')
      CALL KUACT('TKINE4',OMXCNTR)

      CALL KUCMD(' ','TKINE5','C')
      CALL KUPAR('TKINE5','tkin1','Type of simulated particle','I','S')
      CALL KUPVAL('TKINE5','tkin1',8,0.,' ','D')
      CALL KUPAR('TKINE5','tkin2',
     +'Number of simulated particles of this type','I','S')
      CALL KUPVAL('TKINE5','tkin2',1,0.,' ','D')
      CALL KUPAR('TKINE5','tkin3','P minimum /GeV','R','S')
      CALL KUPVAL('TKINE5','tkin3',0,1.,' ','D')
      CALL KUPAR('TKINE5','tkin4','P maximum /GeV','R','S')
      CALL KUPVAL('TKINE5','tkin4',0,100.,' ','D')
      CALL KUPAR('TKINE5','tkin5','Pt mean /GeV','R','S')
      CALL KUPVAL('TKINE5','tkin5',0,0.,' ','D')
      GUID(1)=' Parameters of 5-th class of particles to simulate'
      CALL KUGUID('TKINE5',GUID,1,'S')
      CALL KUACT('TKINE5',OMXCNTR)

      CALL KUCMD(' ','OPDET','C')
      CALL KUPAR('OPDET','IUSET','Name of the detector SET','CO','S')
      CALL KUPVAL('OPDET','IUSET',0,0.,'*','D')
      CALL KUPAR('OPDET','IUDET','Name of the detector','CO','S')
      CALL KUPVAL('OPDET','IUDET',0,0.,'*','D')
      GUID(1)='Prints user information on detectors (OPDET [set] [det])'
      CALL KUGUID('OPDET',GUID,1,'S')
      CALL KUACT('OPDET',OMXCNTR)

      CALL KUCMD(' ','OPRIPA','C')
      CALL KUPAR('OPRIPA','IPAR','the particle number, =0 - all','IO',
     +'S')
      CALL KUPVAL('OPRIPA','IPAR',0,0.,' ','D')
      GUID(1)='Prints the particle data (OMPRIPA [ipar])'
      CALL KUGUID('OPRIPA',GUID,1,'S')
      CALL KUACT('OPRIPA',OMXCNTR)

      CALL KUCMD(' ','OPRIEV','C')
      CALL KUPAR('OPRIEV','IFLAG','Flag - what to print','IO','S')
      CALL KUPVAL('OPRIEV','IFLAG',0,0.,' ','D')
      CALL KUPAR('OPRIEV','ECUT','Cut on the energy','RO','S')
      CALL KUPVAL('OPRIEV','ECUT',0,0.,' ','D')
      GUID(1)='Prints the event data (OMPRIEV [iflag] [cut])'
      GUID(2)='           iflag=0 - the header only'
      GUID(3)='                =1 - list of tracks with Ekin>ECUT'
      CALL KUGUID('OPRIEV',GUID,3,'S')
      CALL KUACT('OPRIEV',OMXCNTR)

      CALL KUCMD(' ','OPHITS','C')
      CALL KUPAR('OPHITS','IUSET','Name of the detector SET','CO','S')
      CALL KUPVAL('OPHITS','IUSET',0,0.,'*','D')
      CALL KUPAR('OPHITS','IUDET','Name of the detector','CO','S')
      CALL KUPVAL('OPHITS','IUDET',0,0.,'*','D')
      CALL KUPAR('OPHITS','ITRA','Track number, =0 - all','IO','S')
      CALL KUPVAL('OPHITS','ITRA',0,0.,' ','D')
      CALL KUPAR('OPHITS','OPT','Option: what to print','CO','S')
      CALL KUPVAL('OPHITS','OPT',0,0.,' ','D')
      GUID(1)='Prints hits (OPHITS [set] [det] [itra] [opt])'
      GUID(2)='       OPT='' '' - call GEANT routine GPHITS (all tracks)
     +'
      GUID(3)='           ''D'' - prints GEANT banks from OMGEANT, inclu
     +de the reference'
      GUID(4)='                 to digitisations'
      CALL KUGUID('OPHITS',GUID,4,'S')
      CALL KUACT('OPHITS',OMXCNTR)

      CALL KUCMD(' ','OPDIGI','C')
      CALL KUPAR('OPDIGI','IUSET','Name of the detector SET','CO','S')
      CALL KUPVAL('OPDIGI','IUSET',0,0.,'*','D')
      CALL KUPAR('OPDIGI','IUDET','Name of the detector','CO','S')
      CALL KUPVAL('OPDIGI','IUDET',0,0.,'*','D')
      CALL KUPAR('OPDIGI','ITRA','Track number, =0 - all','IO','S')
      CALL KUPVAL('OPDIGI','ITRA',0,0.,' ','D')
      CALL KUPAR('OPDIGI','OPT','Option: what to print','CO','S')
      CALL KUPVAL('OPDIGI','OPT',0,0.,' ','D')
      GUID(1)='Prints digitizings (OPDIGI [set] [det] [itra] [opt] )'
      GUID(2)='       OPT='' '' - call GEANT routine GPDIGI (all tracks)
     +'
      GUID(3)='           ''G'' - prints GEANT banks from OMGEANT'
      GUID(4)='           ''O'' - prints OMGEANT banks (all detectors)'
      CALL KUGUID('OPDIGI',GUID,4,'S')
      CALL KUACT('OPDIGI',OMXCNTR)

      CALL KUCMD(' ','OPOHIT','C')
      CALL KUPAR('OPOHIT','ITRA','Track number, =0 - all','IO','S')
      CALL KUPVAL('OPOHIT','ITRA',0,0.,' ','D')
      GUID(1)='Prints local hit bank (OPOHIT [itra])'
      CALL KUGUID('OPOHIT',GUID,1,'S')
      CALL KUACT('OPOHIT',OMXCNTR)

      CALL KUCMD(' ','OPTRIG','C')
      CALL KUPAR('OPTRIG','IFLA','Flag (not used at the moment)','IO',
     +'S')
      CALL KUPVAL('OPTRIG','IFLA',0,0.,' ','D')
      GUID(1)='Prints the trigger detectors hits (OPTRIG [ifla])'
      CALL KUGUID('OPTRIG',GUID,1,'S')
      CALL KUACT('OPTRIG',OMXCNTR)

      CALL KUCMD(' ','OPLUND','C')
      CALL KUPAR('OPLUND','IFLA','flag for LULIST','IO','S')
      CALL KUPVAL('OPLUND','IFLA',1,0.,' ','D')
      GUID(1)='Prints the LUND structure using LULIST(ifla) (OPLUND [ifl
     +a])'
      CALL KUGUID('OPLUND',GUID,1,'S')
      CALL KUACT('OPLUND',OMXCNTR)

      CALL KUCMD(' ','PHYGET','C')
      GUID(1)='Read a PHYNIX event (/PHTRAC/) from LUN=21'
      CALL KUGUID('PHYGET',GUID,1,'S')
      CALL KUACT('PHYGET',OMXCNTR)

      CALL KUCMD(' ','NTOPROC','C')
      CALL KUPAR('NTOPROC','BEAM','!=0 - fill the beam ntuple','IO','S')
      CALL KUPVAL('NTOPROC','BEAM',0,0.,' ','D')
      CALL KUPAR('NTOPROC','KINE','!=0 - fill the kinematics ntuple',
     +'IO','S')
      CALL KUPVAL('NTOPROC','KINE',0,0.,' ','D')
      CALL KUPAR('NTOPROC','HITS','!=0 - fill the hits ntuple','IO','S')
      CALL KUPVAL('NTOPROC','HITS',0,0.,' ','D')
      GUID(1)='  KINE=-1  - write all the vertices and tracks'
      GUID(2)='  KINE=-2  - write only the vertex #1 and the tracks stem
     +ming from it'
      GUID(3)='  KINE=n>0 - write only the track #n and all descendent t
     +racks and verices'
      CALL KUGUID('NTOPROC',GUID,3,'S')
      CALL KUACT('NTOPROC',OMXCNTR)

      CALL KUCMD(' ','NTOOPEN','C')
      GUID(1)='  Opens the output NTUPLE: see NTOPROC'
      CALL KUGUID('NTOOPEN',GUID,1,'S')
      CALL KUACT('NTOOPEN',OMXCNTR)

      CALL KUCMD(' ','NTOCLOSE','C')
      GUID(1)='  Close the output NTUPLE'
      CALL KUGUID('NTOCLOSE',GUID,1,'S')
      CALL KUACT('NTOCLOSE',OMXCNTR)

      CALL KUCMD(' ','NTIPROC','C')
      CALL KUPAR('NTIPROC','BEAM','!=0 - read the beam ntuple','IO','S')
      CALL KUPVAL('NTIPROC','BEAM',0,0.,' ','D')
      CALL KUPAR('NTIPROC','KINE','!=0 - read the kinematics ntuple',
     +'IO','S')
      CALL KUPVAL('NTIPROC','KINE',0,0.,' ','D')
      CALL KUPAR('NTIPROC','HITS','!=0 - read the hits ntuple','IO','S')
      CALL KUPVAL('NTIPROC','HITS',0,0.,' ','D')
      CALL KUPAR('NTIPROC','LUND','!=0 - read the LUND ntuple','IO','S')
      CALL KUPVAL('NTIPROC','LUND',0,0.,' ','D')
      GUID(1)='  LUND=1   - instead of simulating a LUND event - take it
     + from the NTUPLE'
      CALL KUGUID('NTIPROC',GUID,1,'S')
      CALL KUACT('NTIPROC',OMXCNTR)

      CALL KUCMD(' ','NTIOPEN','C')
      GUID(1)='  Opens the input NTUPLE: see NTIPROC'
      CALL KUGUID('NTIOPEN',GUID,1,'S')
      CALL KUACT('NTIOPEN',OMXCNTR)

      CALL KUCMD(' ','NTICLOSE','C')
      GUID(1)='  Close the input NTUPLE'
      CALL KUGUID('NTICLOSE',GUID,1,'S')
      CALL KUACT('NTICLOSE',OMXCNTR)

      CALL KUCMD(' ','GETSEED','C')
      CALL KUPAR('GETSEED','EVENT',
     +' serial event (entry) number in the file','I','S')
      CALL KUPVAL('GETSEED','EVENT',0,0.,' ','D')
      CALL KUPVAL('GETSEED','EVENT',1,0.,' ','L')
      CALL KUPVAL('GETSEED','EVENT',999999999,0.,' ','H')
      CALL KUPAR('GETSEED','IPRI',' <>0 - print the seed values','I',
     +'S')
      CALL KUPVAL('GETSEED','IPRI',0,0.,' ','D')
      CALL KUPAR('GETSEED','NAME',' filename (up to 80 characters)','C',
     +'S')
      CALL KUPVAL('GETSEED','NAME',0,0.,'none','D')
      GUID(1)='  Reads a set of random generator seeds to re-generate a'
     +//' particular event.'
      GUID(2)='  One has to define the filename (direct access unformatt
     +ed file, written by COMGEANT,'
      GUID(3)='  and the serial event (entry) number.'
      CALL KUGUID('GETSEED',GUID,3,'S')
      CALL KUACT('GETSEED',OUXCNTR)

      CALL KUCMD('..',' ','SW')

      CALL KUCMD(' ','OPLOT','C')
      GUID(1)='Plot/print different values'
      CALL KUGUID('OPLOT',GUID,1,'S')

      CALL KUCMD('OPLOT',' ','SW')

      CALL KUCMD(' ','PMAGF','C')
      CALL KUPAR('PMAGF','ID','ID of the histogram. Should be booked',
     +'I','S')
      CALL KUPAR('PMAGF','MAG','Number of the magnet','I','S')
      CALL KUPVAL('PMAGF','MAG',1,0.,' ','D')
      CALL KUPVAL('PMAGF','MAG',1,0.,' ','L')
      CALL KUPVAL('PMAGF','MAG',5,0.,' ','H')
      CALL KUPAR('PMAGF','XH',
     +'Projection X of the 1/2-dim histogram (X/Y/Z=1/2/3)','I','S')
      CALL KUPVAL('PMAGF','XH',1,0.,' ','D')
      CALL KUPVAL('PMAGF','XH',1,0.,' ','L')
      CALL KUPVAL('PMAGF','XH',3,0.,' ','H')
      CALL KUPAR('PMAGF','YH',
     +'Projection Y of the   2-dim histogram (X/Y/Z=1/2/3)','I','S')
      CALL KUPVAL('PMAGF','YH',2,0.,' ','D')
      CALL KUPVAL('PMAGF','YH',1,0.,' ','L')
      CALL KUPVAL('PMAGF','YH',3,0.,' ','H')
      CALL KUPAR('PMAGF','X',
     +'X-coordinate (in the magnet reference frame)','R','S')
      CALL KUPVAL('PMAGF','X',0,0.,' ','D')
      CALL KUPAR('PMAGF','Y','Y-coordinate','R','S')
      CALL KUPVAL('PMAGF','Y',0,0.,' ','D')
      CALL KUPAR('PMAGF','Z','Z-coordinate','R','S')
      CALL KUPVAL('PMAGF','Z',0,0.,' ','D')
      CALL KUPAR('PMAGF','AX',' The plotted value: B=BX*AX+BY*AY+BZ*AZ',
     +'R','S')
      CALL KUPVAL('PMAGF','AX',0,0.,' ','D')
      CALL KUPAR('PMAGF','AY',' The plotted value: B=BX*AX+BY*AY+BZ*AZ',
     +'R','S')
      CALL KUPVAL('PMAGF','AY',0,0.,' ','D')
      CALL KUPAR('PMAGF','AZ',' The plotted value: B=BX*AX+BY*AY+BZ*AZ',
     +'R','S')
      CALL KUPVAL('PMAGF','AZ',0,1.,' ','D')
      GUID(1)='  Plots 1 or 2-dimensional plot of the magnetic field in'
     +//' magnets 1-5'
      GUID(2)='  One may choose a field projection to plot.'
      GUID(3)='  The field component can be protted versus X,Y or Z'
      GUID(4)='    (or any pair for 2-dim plots).'
      CALL KUGUID('PMAGF',GUID,4,'S')
      CALL KUACT('PMAGF',OMXPLOT)

      CALL KUCMD(' ','PRIMAG','C')
      CALL KUPAR('PRIMAG','X','X coordinate of the point','R','S')
      CALL KUPVAL('PRIMAG','X',0,0.,' ','D')
      CALL KUPAR('PRIMAG','Y','Y coordinate of the point','R','S')
      CALL KUPVAL('PRIMAG','Y',0,0.,' ','D')
      CALL KUPAR('PRIMAG','Z','Z coordinate of the point','R','S')
      CALL KUPVAL('PRIMAG','Z',0,0.,' ','D')
      GUID(1)='  Prints the magnetic field in a given point'
      CALL KUGUID('PRIMAG',GUID,1,'S')
      CALL KUACT('PRIMAG',OUXCNTR)

      CALL KUCMD(' ','PRIMED','C')
      CALL KUPAR('PRIMED','X','X coordinate of the point','R','S')
      CALL KUPVAL('PRIMED','X',0,0.,' ','D')
      CALL KUPAR('PRIMED','Y','Y coordinate of the point','R','S')
      CALL KUPVAL('PRIMED','Y',0,0.,' ','D')
      CALL KUPAR('PRIMED','Z','Z coordinate of the point','R','S')
      CALL KUPVAL('PRIMED','Z',0,0.,' ','D')
      GUID(1)='  Prints for a given point: medium, volume names'
      CALL KUGUID('PRIMED',GUID,1,'S')
      CALL KUACT('PRIMED',OUXCNTR)

      CALL KUCMD(' ','PSCATRAJ','C')
      CALL KUPAR('PSCATRAJ','ID',
     +'ID of the histogram. Should be booked','I','S')
      CALL KUPAR('PSCATRAJ','PARA',
     +'Type of the value histogrammed: 0 - RL, 1 - Int.len.','I','S')
      CALL KUPVAL('PSCATRAJ','PARA',0,0.,' ','D')
      CALL KUPVAL('PSCATRAJ','PARA',0,0.,' ','L')
      CALL KUPVAL('PSCATRAJ','PARA',1,0.,' ','H')
      CALL KUPAR('PSCATRAJ','TYPE','Particle type','I','S')
      CALL KUPVAL('PSCATRAJ','TYPE',14,0.,' ','D')
      CALL KUPVAL('PSCATRAJ','TYPE',1,0.,' ','L')
      CALL KUPVAL('PSCATRAJ','TYPE',20,0.,' ','H')
      CALL KUPAR('PSCATRAJ','X','X-coordinate of the track start','R',
     +'S')
      CALL KUPVAL('PSCATRAJ','X',0,0.,' ','D')
      CALL KUPAR('PSCATRAJ','Y','Y-coordinate','R','S')
      CALL KUPVAL('PSCATRAJ','Y',0,0.,' ','D')
      CALL KUPAR('PSCATRAJ','Z','Z-coordinate','R','S')
      CALL KUPVAL('PSCATRAJ','Z',0,0.,' ','D')
      CALL KUPAR('PSCATRAJ','PX','PX - of the momentum','R','S')
      CALL KUPVAL('PSCATRAJ','PX',0,300.,' ','D')
      CALL KUPAR('PSCATRAJ','PY','PY - of the momentum','R','S')
      CALL KUPVAL('PSCATRAJ','PY',0,0.,' ','D')
      CALL KUPAR('PSCATRAJ','PZ','PZ - of the momentum','R','S')
      CALL KUPVAL('PSCATRAJ','PZ',0,0.,' ','D')
      GUID(1)='  Fills a histogram with the amount of material in RL/..'
     +//' along'
      GUID(2)='  the length of a track for the next event simulated. If'
     +//' this'
      GUID(3)='  flag is set then only one track is simulated for the ne
     +xt event.'
      CALL KUGUID('PSCATRAJ',GUID,3,'S')
      CALL KUACT('PSCATRAJ',OMXPLOT)

      CALL KUCMD(' ','PDEVIATION','C')
      CALL KUPAR('PDEVIATION','IDT','IDTRID identifier of the plane',
     +'I','S')
      CALL KUPAR('PDEVIATION','ID1',
     +'ID of the histogram for slope deviations. Should be booked','I',
     +'S')
      CALL KUPAR('PDEVIATION','ID2',
     +'ID of the histogram for position deviations.','I','S')
      CALL KUPAR('PDEVIATION','ID3',
     +'ID of the 2-dim histogram for slope-position deviations.','I',
     +'S')
      GUID(1)='  Fills histograms with track deviation from a straight l
     +ine'
      GUID(2)='  measured at the exit from a specified detector'
      CALL KUGUID('PDEVIATION',GUID,2,'S')
      CALL KUACT('PDEVIATION',OMXPLOT)

      CALL KUCMD('..',' ','SW')

      CALL KUCMD(' ','MODIFICATION','C')
      GUID(1)='Modifications of some preset values'
      CALL KUGUID('MODIFICATION',GUID,1,'S')

      CALL KUCMD('MODIFICATION',' ','SW')

      CALL KUCMD(' ','MAGFIELD','C')
      CALL KUPAR('MAGFIELD','MAF','the magnet number','I','S')
      CALL KUPAR('MAGFIELD','FIELD','Field kGs','R','S')
      CALL KUPVAL('MAGFIELD','FIELD',0,-9999.,' ','D')
      CALL KUPAR('MAGFIELD','FLAG1',
     +'For dipole: flag1=1 - longitudinal component included','RO','S')
      CALL KUPVAL('MAGFIELD','FLAG1',0,-9999.,' ','D')
      CALL KUPAR('MAGFIELD','FLAG2','flag2','RO','S')
      CALL KUPVAL('MAGFIELD','FLAG2',0,-9999.,' ','D')
      GUID(1)='  Changes the values of the field (AMAGPAR(5-7,mag))'
      CALL KUGUID('MAGFIELD',GUID,1,'S')
      CALL KUACT('MAGFIELD',OUXCNTR)

      CALL KUCMD(' ','MAGSIZE','C')
      CALL KUPAR('MAGSIZE','MAF','the magnet number','I','S')
      CALL KUPAR('MAGSIZE','X','Full size in X','R','S')
      CALL KUPVAL('MAGSIZE','X',0,-999.,' ','D')
      CALL KUPAR('MAGSIZE','Y','Full size in Y','R','S')
      CALL KUPVAL('MAGSIZE','Y',0,-999.,' ','D')
      CALL KUPAR('MAGSIZE','Z','Full size in Z','R','S')
      CALL KUPVAL('MAGSIZE','Z',0,-999.,' ','D')
      GUID(1)='  Changes the dimensions of the magnetic field (PMOLMAG(2
     +-4,mag)'
      GUID(2)='  Negative values are ignored.'
      CALL KUGUID('MAGSIZE',GUID,2,'S')
      CALL KUACT('MAGSIZE',OUXCNTR)

      CALL KUCMD(' ','MOLLANGL','C')
      CALL KUPAR('MOLLANGL','THETMIN','theta minimum','R','S')
      CALL KUPVAL('MOLLANGL','THETMIN',0,0.,' ','D')
      CALL KUPAR('MOLLANGL','THETMAX','theta maximum','R','S')
      CALL KUPVAL('MOLLANGL','THETMAX',0,0.,' ','D')
      CALL KUPAR('MOLLANGL','PHIMIN','phi minimum','R','S')
      CALL KUPVAL('MOLLANGL','PHIMIN',0,0.,' ','D')
      CALL KUPAR('MOLLANGL','PHIMAX','phi maximum','R','S')
      CALL KUPVAL('MOLLANGL','PHIMAX',0,0.,' ','D')
      GUID(1)='  Changes the range of Moller scattering angles in CM'
      CALL KUGUID('MOLLANGL',GUID,1,'S')
      CALL KUACT('MOLLANGL',OUXCNTR)

      CALL KUCMD(' ','MOTTANGL','C')
      CALL KUPAR('MOTTANGL','THETMIN','theta minimum','R','S')
      CALL KUPVAL('MOTTANGL','THETMIN',0,0.,' ','D')
      CALL KUPAR('MOTTANGL','THETMAX','theta maximum','R','S')
      CALL KUPVAL('MOTTANGL','THETMAX',0,0.,' ','D')
      CALL KUPAR('MOTTANGL','PHIMIN','phi minimum','R','S')
      CALL KUPVAL('MOTTANGL','PHIMIN',0,0.,' ','D')
      CALL KUPAR('MOTTANGL','PHIMAX','phi maximum','R','S')
      CALL KUPVAL('MOTTANGL','PHIMAX',0,0.,' ','D')
      GUID(1)='  Changes the range of Mott scattering angles in Lab'
      CALL KUGUID('MOTTANGL',GUID,1,'S')
      CALL KUACT('MOTTANGL',OUXCNTR)

      CALL KUCMD(' ','BEAMMOM','C')
      CALL KUPAR('BEAMMOM','B1','minimum momentum GeV/c','R','S')
      CALL KUPAR('BEAMMOM','B1','minimum momentum GeV/c','R','S')
      CALL KUPAR('BEAMMOM','SB','sigma of the distribution','R','S')
      GUID(1)='  Redefines the beam momentum (it the beam is defined as'
     +//' a gaussiqan'
      CALL KUGUID('BEAMMOM',GUID,1,'S')
      CALL KUACT('BEAMMOM',OUXCNTR)

      CALL KUCMD(' ','CHERFAC','C')
      CALL KUPAR('CHERFAC','FAC',
     +'A factor for Cherenkov photons multiplicity','R','S')
      CALL KUPVAL('CHERFAC','FAC',0,369.81E9,' ','D')
      GUID(1)='  Redefines the multiplicity for Cherenkov photons'
      CALL KUGUID('CHERFAC',GUID,1,'S')
      CALL KUACT('CHERFAC',OUXCNTR)

      CALL KUCMD(' ','KINONE','C')
      CALL KUPAR('KINONE','IPAR',
     +'A factor for Cherenkov photons multiplicity','I','S')
      CALL KUPVAL('KINONE','IPAR',1,0.,' ','L')
      CALL KUPVAL('KINONE','IPAR',200,0.,' ','H')
      CALL KUPAR('KINONE','MOM','Central momentum, GeV','R','S')
      CALL KUPVAL('KINONE','MOM',0,0.,' ','L')
      CALL KUPVAL('KINONE','MOM',0,200.,' ','H')
      CALL KUPAR('KINONE','THET','Central theta angle (to Z axis), rad',
     +'R','S')
      CALL KUPVAL('KINONE','THET',0,0.,' ','L')
      CALL KUPVAL('KINONE','THET',0,3.1416,' ','H')
      CALL KUPAR('KINONE','PHI',
     +'Central phi   angle (proj onto XY, angle to X), rad','R','S')
      CALL KUPVAL('KINONE','PHI',0,0.,' ','L')
      CALL KUPVAL('KINONE','PHI',0,6.29,' ','H')
      CALL KUPAR('KINONE','DMOM','Momentum width: p=p0+/-dp , GeV','R',
     +'S')
      CALL KUPVAL('KINONE','DMOM',0,0.,' ','L')
      CALL KUPVAL('KINONE','DMOM',0,200.,' ','H')
      CALL KUPAR('KINONE','DTHET','Theta    width, rad','R','S')
      CALL KUPVAL('KINONE','DTHET',0,0.,' ','L')
      CALL KUPVAL('KINONE','DTHET',0,3.1416,' ','H')
      CALL KUPAR('KINONE','DPHI','Phi      width, rad','R','S')
      CALL KUPVAL('KINONE','DPHI',0,0.,' ','L')
      CALL KUPVAL('KINONE','DPHI',0,6.29,' ','H')
      GUID(1)='  Parameters of the single particle to generate (OMKINONE
     +)'
      CALL KUGUID('KINONE',GUID,1,'S')
      CALL KUACT('KINONE',OUXCNTR)

      CALL KUCMD('..',' ','SW')

      CALL KUCMD('..',' ','SW')

      CALL KUCMD('/',' ','SW')

      END
