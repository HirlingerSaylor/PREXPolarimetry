#!/usr/bin/python

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

# Sort energy data into left arm (channels 1-4), right arm (channels 5-8),
# and total (channels 1-8).
# "left" and "right" are somewhat ambiguous for a symmetric detector, and any
# distinction seems inconsequential, so lets leave it that way (for now)
ene_l = np.zeros(len(dEne))
ene_r = np.zeros(len(dEne))
ene_t = np.zeros(len(dEne))
for i in xrange(len(dEne)):
	ene_l[i] = np.sum(dEne[i,0:4])
	ene_r[i] = np.sum(dEne[i,4:8])
	ene_t[i] = np.sum(dEne[i,0:8])

# Create figures
ax1 = plt.subplot(131)
ax1.plot(dCoin[:,0] *180. /np.pi, ene_l, 'bo', markersize=5)
ax1.set_xlabel(r'Scattering Angle $\theta ^\circ$', fontsize=24)
ax1.set_ylabel('Energy (MeV)', fontsize=24)
ax1.set_title('Left Arm Energy', fontsize=26)
ax1.set_xlim((70,110))
ax1.tick_params(axis='both', which='major',labelsize=22)
ax1.tick_params(axis='both', which='minor',labelsize=22)


ax2 = plt.subplot(132)
ax2.plot(dCoin[:,0] *180. /np.pi, ene_r, 'bo', markersize=5)
ax2.set_xlabel(r'Scattering Angle $\theta ^\circ$', fontsize=24)
ax2.set_ylabel('Energy (MeV)', fontsize=24)
ax2.set_title('Right Arm Energy', fontsize=26)
ax2.set_xlim((70,110))
ax2.tick_params(axis='both', which='major',labelsize=22)
ax2.tick_params(axis='both', which='minor',labelsize=22)

ax3 = plt.subplot(133)
ax3.plot(dCoin[:,0] *180. /np.pi, ene_t, 'bo', markersize=5)
ax3.set_xlabel(r'Scattering Angle $\theta ^\circ$', fontsize=24)
ax3.set_ylabel('Energy (MeV)', fontsize=24)
ax3.set_title('Total Energy', fontsize=26)
ax3.set_xlim((70,110))
ax3.tick_params(axis='both', which='major',labelsize=22)
ax3.tick_params(axis='both', which='minor',labelsize=22)

plt.suptitle(r'Energy vs. Scattering Angle $\theta$', fontsize=28)
plt.show()
