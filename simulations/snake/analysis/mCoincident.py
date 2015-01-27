#!/usr/bin/python
# Locate Moller scattering events where both electrons are detected

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat


# Analyzing Power
def anp(th):
	return ( (np.sin(th*np.pi/180.)**2. *(7+np.cos(th*np.pi/180.)**2.))
						/((3+np.cos(th*np.pi/180.)**2.)**2.) )

# Forward scattered electron trajectories
fD = np.loadtxt('dFse.txt')

# Reverse scattered electron trajectories
rD = np.loadtxt('dRse.txt')

# Scattering angle data
aD = np.loadtxt('dAngs.txt')


# Sort trajectory data
#  Note that maxep if the endplane number of the final detector aperture,
#  this could change if the directive file is modified, and could give very
#  confusing results if it isn't modified here accordingly.
maxep = 25
fse = np.zeros((max(fD[:,0]),max(fD[:,1]),3))
rse = np.zeros((max(rD[:,0]),max(rD[:,1]),3))

if len(fse) != len(rse):
	print('error: different numbers of forward and reverse scattered electrons')
if len(fse[0]) != len(rse[0]):
	print('error: different numbers of trajectory points')
if len(fse) != len(aD):
	print('error: different number of trajectories and scattering angles')

for e in fD:
	fse[e[0]-1,e[1]-1] = e[2:5]
for e in rD:
	rse[e[0]-1,e[1]-1] = e[2:5]

# Find coincident events
coin = open('dCoincident.txt','w')
for i in xrange(len(fse)):
	if all(fse[i,maxep]!=0.) and all(rse[i,maxep]!=0.):
		coin.write( str(aD[i,0])+' '+str(aD[i,1])+' '+str(anp(aD[i,0]))+' '
			+str(fse[i,maxep,0])+' '+str(rse[i,maxep,0])+' '+str(fse[i,maxep,1])+' '
			+str(rse[i,maxep,1])+' '+str(fse[i,maxep,2])+' '+str(rse[i,maxep,2])+'\n')

