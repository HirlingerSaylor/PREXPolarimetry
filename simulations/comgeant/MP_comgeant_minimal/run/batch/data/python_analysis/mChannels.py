#!/usr/bin/python
# Divided coincident data based on energy in specific detector channels
# Provides method to use channels as apertures
# Could look funny if dipole not correctly adjusted such that high anp events
# impact the vertical center of a channel pair, keep in mind that the vertical
# center of the detector is a channel boundary

import numpy as np
import matplotlib.pyplot as plt

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dCoin = np.loadtxt('dCoincidentData.txt')

# Load energy data
#  format (energy in channel 1,2,3,4,5,6,7,8)
#  units are MeV
dEne = np.loadtxt('dEnergyByChannels.txt')

# Files for data restriced only by left arm channel requirements
#  right arm results SHOULD be the same
dEneL1 = open('dDataChannel1.txt','w')
dEneL2 = open('dDataChannel2.txt','w')
dEneL3 = open('dDataChannel3.txt','w')
dEneL4 = open('dDataChannel4.txt','w')

# Files for data restricted by channel pairs
dEne_1_5 = open('dDataChannel_1_5.txt','w')
dEne_2_6 = open('dDataChannel_2_6.txt','w')
dEne_3_7 = open('dDataChannel_3_7.txt','w')
dEne_4_8 = open('dDataChannel_4_8.txt','w')

# Sort coincident events based on channels with largest energy
#  could (should?) be changed to some majority percent of total energy
for i in xrange(len(dCoin)):

	# Locate channels with highest energies
	q = np.argmax(dEne[i,0:4])
	p = np.argmax(dEne[i,4:8])

	# Build data string
	dList=''
	for j in xrange(len(dCoin[0])):
		dList += str(dCoin[i,j])
		if j < len(dCoin[0]) -1:
			dList += ' '
		else:
			dList += '\n'

	# Disperse data string
	if q == 0:
		dEneL1.write(dList)
	elif q==1:
		dEneL2.write(dList)
	elif q==2:
		dEneL3.write(dList)
	elif q==3:
		dEneL4.write(dList)

	if q==0 and p==0:
		dEne_1_5.write(dList)
	elif q==1 and p==1:
		dEne_2_6.write(dList)
	elif q==2 and p==2:
		dEne_3_7.write(dList)
	elif q==3 and p==3:
		dEne_4_8.write(dList)

dEneL1.close()
dEneL2.close()
dEneL3.close()
dEneL4.close()

dEne_1_5.close()
dEne_2_6.close()
dEne_3_7.close()
dEne_4_8.close()
