#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dCoin = np.loadtxt('dCoincidentData.txt')

# Load channel 2/6 specific coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dChan = np.loadtxt('dDataChannel_2_6.txt')

# Display average analyzing powers
print('All Coincident events <Anp>='+str(np.mean(dCoin[:,2])))
print('Channel 2-6 Coincident events <Anp>='+str(np.mean(dChan[:,2])))

# Generate histograms
b1 = np.linspace(0.7,0.8,30)

ax1 = plt.subplot(111)
ax1.hist(dCoin[:,2],bins=b1,range=[0.7,0.8],histtype='step',lw=4,
					label='All Coincident Events')
ax1.hist(dChan[:,2],bins=b1,range=[0.7,0.8],histtype='step',lw=4,
					ls='dashed',label='Channel 2-6 Coincident Events')
ax1.set_xlabel('Analyzing Power', fontsize=24)
ax1.set_title('Analyzing Power Histograms',fontsize=26)
ax1.tick_params(axis='both', which='major',labelsize=20)
ax1.tick_params(axis='both', which='minor',labelsize=20)
ax1.legend()

plt.show()
