#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dCoin = np.loadtxt('dCoincidentData.txt')

# Display number of coincident events, as well as angle ranges
print len(dCoin)
print('theta:'+str(min(dCoin[:,0])*180./np.pi)+','
			+str(max(dCoin[:,0])*180./np.pi))
print('phi:'+str(min(dCoin[:,1])*180./np.pi)+','
			+str(max(dCoin[:,1])*180./np.pi))

# Generate figures
ax1 = plt.subplot(121)
ax1.hist(dCoin[:,0] *180. /np.pi, bins=20)
ax1.set_xlim((70,110))
ax1.set_title(r"Scattering Angle $\theta$", fontsize=26)
ax1.set_xlabel(r"$\theta ^\circ$", fontsize=24)
ax1.tick_params(axis='both', which='major',labelsize=22)
ax1.tick_params(axis='both', which='minor',labelsize=22)

ax2 = plt.subplot(122)
ax2.hist(dCoin[:,1] *180. /np.pi, bins=20)
ax2.set_xlim((-10,10))
ax2.set_title(r"Scattering Angle $\phi$", fontsize=26)
ax2.set_xlabel(r"$\phi ^\circ$", fontsize=24)
ax2.tick_params(axis='both', which='major',labelsize=22)
ax2.tick_params(axis='both', which='minor',labelsize=22)

plt.suptitle('COMGEANT Angular Acceptance Histograms: '
							+str(len(dCoin))+' Coincident Events',fontsize=26)
plt.show()
