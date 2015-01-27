#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dCoin = np.loadtxt('dCoincidentData.txt')

# Define positions at detector
#  incorporates 7 degree rotation
xlDet = dCoin[:,3]
xrDet = dCoin[:,4]

ylDet = ((dCoin[:,5] +44.9) *np.cos(7 *np.pi /180.) 
				+(dCoin[:,7] -692.6) *np.sin(7 *np.pi /180.))
yrDet = ((dCoin[:,6] +44.9) *np.cos(7 *np.pi /180.) 
				+(dCoin[:,8] -692.6) *np.sin(7 *np.pi /180.))

# Generate figure
ax1 = plt.subplot(111, aspect=1)

im = ax1.scatter(xlDet,ylDet,s=3,c=dCoin[:,2],marker='x',vmin=0.7,vmax=0.78)
ax1.scatter(xrDet,yrDet,s=3,c=dCoin[:,2],marker='x',vmin=0.7,vmax=0.78)

# Detector aperture
#  see include/geom_global.ffr for gorey details
ax1.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax1.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)

ax1.set_xlim((-10,10))
ax1.set_ylim((-20,20))
ax1.set_xlabel('X(cm)', fontsize=24)
ax1.set_ylabel('Y(cm)', fontsize=24)
ax1.set_title('GEANT Position of Electrons on Detector', fontsize=26)
thebar = plt.colorbar(im)
thebar.set_label('Analyzing Power', size=24)
ax1.tick_params(axis='both', which='major',labelsize=20)
ax1.tick_params(axis='both', which='minor',labelsize=20)


plt.show()
