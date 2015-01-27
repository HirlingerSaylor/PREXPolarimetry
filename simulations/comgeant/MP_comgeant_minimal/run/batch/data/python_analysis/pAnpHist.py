#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dCoin = np.loadtxt('dCoincidentData.txt')

# Generate figure
ax = plt.subplot(111)

ax.hist(dCoin[:,2], 50, range=[0.65,0.8])
ax.set_xlabel('Analyzing Power', fontsize=24)
ax.set_title('GEANT Histogram of Analyzing Power: <Anp>='
						+str(np.mean(dCoin[:,2])), fontsize=26)
ax.tick_params(axis='both', which='major',labelsize=20)
ax.tick_params(axis='both', which='minor',labelsize=20)

plt.show()
