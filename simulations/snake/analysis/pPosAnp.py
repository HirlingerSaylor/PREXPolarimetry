#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat

# Load coincident data
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in degrees, lengths in mm, and anp is analyzing power
dCoin = np.loadtxt('dCoincident.txt')

# Define positions at detector
#  incorporates 7 degree rotation
xlDet = dCoin[:,7]
xrDet = dCoin[:,8]

ylDet = ((dCoin[:,3] +451.18) *np.cos(7 *np.pi /180.) 
				+(dCoin[:,5] -6872.87) *np.sin(7 *np.pi /180.))
yrDet = ((dCoin[:,4] +451.18) *np.cos(7 *np.pi /180.) 
				+(dCoin[:,6] -6872.87) *np.sin(7 *np.pi /180.))

# Define detector apertures
la = [ [-23.,155.],[-63.,155.],[-63.,-155.],[-23.,-155.] ]
ra = [ [63.,155],[23.,155.],[23.,-155.],[63.,-155.] ]

# Generate figure
ax1 = plt.subplot(111, aspect=1)

im = ax1.scatter(xlDet,ylDet,s=3,c=dCoin[:,2],marker='x',vmin=0.7,vmax=0.78)
ax1.scatter(xrDet,yrDet,s=3,c=dCoin[:,2],marker='x',vmin=0.7,vmax=0.78)

# Detector aperture
ax1.add_patch(pat.Polygon(la,color='k',lw=2,fill=False))
ax1.add_patch(pat.Polygon(ra,color='k',lw=2,fill=False))

ax1.set_xlim((-100,100))
ax1.set_ylim((-200,200))
ax1.set_xlabel('X(mm)', fontsize=24)
ax1.set_ylabel('Y(mm)', fontsize=24)
ax1.set_title('SNAKE Position of Electrons on Detector', fontsize=26)
thebar = plt.colorbar(im)
thebar.set_label('Analyzing Power', size=24)
ax1.tick_params(axis='both', which='major',labelsize=20)
ax1.tick_params(axis='both', which='minor',labelsize=20)

plt.show()
