#!/usr/bin/python
# Polar plot of initial and coincident events, to visualize acceptance

import numpy as np
import matplotlib.pyplot as plt

# Load initial data
#  format (thetcm,phicm)
#  angles are in degrees
dInit = np.loadtxt('dAngs.txt')

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in degrees, lengths in mm, and anp is analyzing power
dCoin = np.loadtxt('dCoincident.txt')

# Generate figures
ax1 = plt.subplot(121,aspect=1)
ax1.plot(dInit[:,0] *np.cos(dInit[:,1] *np.pi /180.), 
					dInit[:,0] *np.sin(dInit[:,1] *np.pi /180.), 'bo', markersize=3)
ax1.set_title(str(len(dInit))+" Initial Events", fontsize=26)
ax1.set_xlim((-130,130))
ax1.set_ylim((-130,130))
ax1.tick_params(axis='both', which='major',labelsize=22)
ax1.tick_params(axis='both', which='minor',labelsize=22)


ax2 = plt.subplot(122,aspect=1)
ax2.plot(dCoin[:,0] *np.cos(dCoin[:,1] *np.pi /180.), 
					dCoin[:,0] *np.sin(dCoin[:,1] *np.pi /180.), 'bo', markersize=3)
ax2.set_title(str(len(dCoin))+" Coincident Events",fontsize=26)
ax2.set_xlim((-130,130))
ax2.set_ylim((-130,130))
ax2.tick_params(axis='both', which='major',labelsize=22)
ax2.tick_params(axis='both', which='minor',labelsize=22)

plt.suptitle('SNAKE Angular Acceptance: Polar Plot',fontsize=28)
plt.show()
