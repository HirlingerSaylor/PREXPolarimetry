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

# Generate histograms
b1 = np.linspace(60,120,30)
b2 = np.linspace(-10,10,30)

ax1 = plt.subplot(121)
ax1.hist(dCoin[:,0] *180. /np.pi,bins=b1,range=[60,120],histtype='step',lw=4,
				label='All Coincident Events')
ax1.hist(dChan[:,0] *180. /np.pi,bins=b1,range=[60,120],histtype='step',lw=4,
				label='Channel 2-6 Coincident Events')
ax1.set_xlim((60,120))
ax1.set_title(r"Scattering Angle $\theta$", fontsize=26)
ax1.set_xlabel(r"$\theta ^\circ$", fontsize=24)
ax1.tick_params(axis='both', which='major',labelsize=22)
ax1.tick_params(axis='both', which='minor',labelsize=22)
ax1.legend()

ax2 = plt.subplot(122)
ax2.hist(dCoin[:,1] *180. /np.pi,bins=b2,range=[-10,10],histtype='step',lw=4,
				label='All Coincident Events')
ax2.hist(dChan[:,1] *180. /np.pi,bins=b2,range=[-10,10],histtype='step',lw=4,
				label='Channel 2-6 Coincident Events')
ax2.set_xlim((-20,20))
ax2.set_title(r"Scattering Angle $\phi$", fontsize=26)
ax2.set_xlabel(r"$\phi ^\circ$", fontsize=24)
ax2.tick_params(axis='both', which='major',labelsize=22)
ax2.tick_params(axis='both', which='minor',labelsize=22)
ax2.legend()

plt.show()

