#!/usr/bin/python

import os
import numpy as np
from subprocess import call

# Delete data from previous run
call('cd data; rm run.*; rm figures/*',shell=True)

# Shortcuts for run, data, & include directories
dRun=os.getcwd()
iDat=os.path.dirname(os.path.dirname(dRun))+'/include'

# Retrieve previous run parameters
iPar=open('omg.ffr').readlines()[197].split()
iPar.pop(0)
vPar=np.zeros(len(iPar))
for x in xrange(len(iPar)):
	vPar[x]=float(iPar[x])
nPar=np.array([	'Number of Events     ',
				'Theta CM center      ',
				'Theta CM range       ',
				'Phi CM center        ',
				'Phi CM range         ',
				'Quad 1 field (kG)    ',
				'Quad 2 field (kG)    ',
				'Quad 3 field (kG)    ',
				'Quad 4 field (kG)    ',
				'Dipole field (kG)    ',
				'Target Field (T)     ',
				'Beam Energy center   ',
				'Beam Energy range    ',
				'Beam x-Position range',
				'Beam y-Position range',
				'Beam x-Slope range   ',
				'Beam y-Slope range   ',])
print('\nCode\tName\t\t\tValue\n')
for x in xrange(17):
	print(str(x)+'\t'+str(nPar[x])+'\t'+str(vPar[x]))

# Permit user to change sim parameters
change=raw_input('\nWould you like to change the run parameters? (y/n) ')
if change=='y':
	while 1:
		new_param=raw_input('\nEnter the parameter code, or "q" to run sim: ')
		if new_param=='q': break
		value=raw_input('\nEnter the new '+str(nPar[int(new_param)])+' : ')
		vPar[int(new_param)]=value
		print('\nCode\tName\t\t\tValue\n')
		for x in xrange(17):
			print(str(x)+'\t'+str(nPar[x])+'\t'+str(vPar[x]))

# Initialize sim run parameters
eve=vPar[0]
theMin=vPar[1]-vPar[2]/2
theMax=vPar[1]+vPar[2]/2
phiMin=vPar[3]-vPar[4]/2
phiMax=vPar[3]+vPar[4]/2
mags=np.array([vPar[5],vPar[6],vPar[7],vPar[8],vPar[10]/4,vPar[9]])
beVMin=vPar[11]-vPar[12]/2
beVMax=vPar[11]+vPar[12]/2
bxpMin=0-vPar[13]/2
bxpMax=vPar[13]/2
bypMin=0-vPar[14]/2
bypMax=vPar[14]/2
bxsMin=0-vPar[15]/2
bxsMax=vPar[15]/2
bysMin=0-vPar[16]/2
bysMax=vPar[16]/2
typMag=np.array([1,1,1,1,9,2])
xMag=np.zeros(6)
yMag=np.zeros(6)
zMag=np.array([75.19,140.46,209.08,274.59,0.0,423.20])
xwMag=np.array([5.08, 5.08, 5.08, 5.08, 20.0, 8.001])
ywMag=([5.08, 5.08, 5.08, 5.08, 20.0, 30.001])
zwMag=([18.29, 22.3, 18.37, 18.37, 80.0, 82.5])
 
# Prepare input file for sim
call('cp '+str(iDat)+'/omg_0.ffr omg.ffr',shell=True)
pFile=open('omg.ffr','a')
pFile.write('BEAMYZLIM '+str(bxpMin)+' '+str(bxpMax)+' '+str(bypMin)+' '
			+str(bypMax)+' '+str(bxsMin)+' '+str(bxsMax)+' '+str(bysMin)
			+' '+str(bxsMax)+'\n'
			+'TRIG '+str(eve)+' events\n'
			+'BEAMOMLIM '+str(beVMin)+' '+str(beVMax)+' momentum limits\n'
			+'MAGNET05 10='+str(vPar[10]/4)+'\n'
			+'MOLLIMTHETA '+str(theMin)+' '+str(theMax)+' theta cm range\n'
			+'MOLLIMPHI '+str(phiMin)+' '+str(phiMax)+' phi cm range\n')
for x in xrange(6):
	if x==4: continue
	pFile.write('MAGNET0'+str(x+1)+' '+str(typMag[x])+" 'HALL' "+str(xMag[x])
			+' '+str(yMag[x])+' '+str(zMag[x])+' 0 '+str(xwMag[x])+' '
			+str(ywMag[x])+' '+str(zwMag[x])+' '+str(mags[x])+' 0 0 Q one\n')
pFile.write('C '+str(vPar[0])+' '+str(vPar[1])+' '+str(vPar[2])+' '
			+str(vPar[3])+' '+str(vPar[4])+' '+str(vPar[5])+' '
			+str(vPar[6])+' '+str(vPar[7])+' '+str(vPar[8])+' '
			+str(vPar[9])+' '+str(vPar[10])+' '+str(vPar[11])+' '
			+str(vPar[12])+' '+str(vPar[13])+' '+str(vPar[14])+' '
			+str(vPar[15])+' '+str(vPar[16])+'\n')
pFile.close()

# Prepare run directory
call('ln -s omg.ffr fort.15;ln -s '+str(iDat)+'/particles.ffr fort.16;ln -s '
			+str(iDat)+'/materials.ffr fort.17;ln -s '+str(iDat)
			+'/geom_global.ffr fort.21;ln -s '
		  +str(iDat)+'/geom_local.ffr fort.22;ln -s '
		  +str(iDat)+'/high_field.map magnet_5.map',shell=True)

# Run Comgeant, store data and log files, and clean up a bit
print('\nNow running COMGEANT!\n')
call('ln -s '+str(iDat)+'/comg.exe comg.exe',shell=True)
call('./comg.exe > csimout.log',shell=True)
call('rm comg.exe fort.* mag*; mv csimout.nt data/run.nt;'
			+'mv omgeant.his data/run.his',shell=True)
call('mv csimout.log data/run.log;rm detectors.dat omg_random.dat',shell=True)

# Run paw analysis
analysis=open('data/analysis.txt','w')
analysis.write('3                                                                              \nexe analysis')
analysis.close()
call('cp analysis.kumac data/analysis.kumac',shell=True)
call('ln -s '+str(iDat)+'/siglg.f data/siglg.f',shell=True)
call('ln -s '+str(iDat)+'/apnum.f data/apnum.f',shell=True)
call('ln -s '+str(iDat)+'/hitraxy.f data/hitraxy.f',shell=True)
call('cd data; cat analysis.txt | paw',shell=True)
call('cd data; mv *.dat python_analysis/; '
			+'rm *.f analysis.* comis* last.kumac',shell=True)
