#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load channel 2/6 specific coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
dChan = np.loadtxt('dDataChannel_2_6.txt')

# Generate positions at detector plane
#  includes 7 degree rotation
xlDet = -dChan[:,3]
xrDet = -dChan[:,4]
ylDet = ((dChan[:,5] +44.9) *np.cos(7 *np.pi /180.) 
				+(dChan[:,7] -692.6) *np.sin(7 *np.pi /180.))
yrDet = ((dChan[:,6] +44.9) *np.cos(7 *np.pi /180.) 
				+(dChan[:,8] -692.6) *np.sin(7 *np.pi /180.))

# Generate figure
ax1 = plt.subplot(111, aspect=1)

im = ax1.scatter(xlDet,ylDet,s=3,c=dChan[:,2],marker='x',vmin=0.7,vmax=0.78)
im = ax1.scatter(xrDet,yrDet,s=3,c=dChan[:,2],marker='x',vmin=0.7,vmax=0.78)
ax1.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax1.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax1.set_xlim((-10,10))
ax1.set_ylim((-20,20))
ax1.set_xlabel('X(cm)', fontsize=24)
ax1.set_ylabel('Y(cm)', fontsize=24)
ax1.set_title('Position of Electrons on Detector', fontsize=26)
ax1.tick_params(axis='both', which='major',labelsize=20)
ax1.tick_params(axis='both', which='minor',labelsize=20)
thebar = plt.colorbar(im)
thebar.set_label('Analyzing Power', size=24)
plt.show()
