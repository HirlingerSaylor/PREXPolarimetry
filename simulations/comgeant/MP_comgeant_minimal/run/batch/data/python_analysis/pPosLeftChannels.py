#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load coincident data specific to left arm channels
#  format (thetcm,phicm,anp,x_l,x_r,y_l,y_r,z_l,z_r)
#  angles are in radians, lengths in cm, and anp is analyzing power
c1 = np.loadtxt('dDataChannel1.txt')
c2 = np.loadtxt('dDataChannel2.txt')
c3 = np.loadtxt('dDataChannel3.txt')
c4 = np.loadtxt('dDataChannel4.txt')

# Generate positions at detector plane
#  includes 7 degree rotation
x1f = -c1[:,3]
x1r = -c1[:,4]
y1f = ((c1[:,5] +44.9) *np.cos(7 *np.pi /180.) 
			+(c1[:,7] -692.6) *np.sin(7 *np.pi /180.))
y1r = ((c1[:,6] +44.9) *np.cos(7 *np.pi /180.) 
			+(c1[:,8] -692.6) *np.sin(7 *np.pi /180.))
anp1 = c1[:,2]

x2f = -c2[:,3]
x2r = -c2[:,4]
y2f = ((c2[:,5] +44.9) *np.cos(7 *np.pi /180.) 
			+(c2[:,7] -692.6) *np.sin(7 *np.pi /180.))
y2r = ((c2[:,6] +44.9) *np.cos(7 *np.pi /180.) 
			+(c2[:,8] -692.6) *np.sin(7 *np.pi /180.))
anp2 = c2[:,2]

x3f = -c3[:,3]
x3r = -c3[:,4]
y3f = ((c3[:,5] +44.9) *np.cos(7 *np.pi /180.) 
			+(c3[:,7] -692.6) *np.sin(7 *np.pi /180.))
y3r = ((c3[:,6] +44.9) *np.cos(7 *np.pi /180.) 
			+(c3[:,8] -692.6) *np.sin(7 *np.pi /180.))
anp3 = c3[:,2]

x4f = -c4[:,3]
x4r = -c4[:,4]
y4f = ((c4[:,5] +44.9) *np.cos(7 *np.pi /180.) 
			+(c4[:,7] -692.6) *np.sin(7 *np.pi /180.))
y4r = ((c4[:,6] +44.9) *np.cos(7 *np.pi /180.) 
			+(c4[:,8] -692.6) *np.sin(7 *np.pi /180.))
anp4 = c4[:,2]

fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)

im = ax1.scatter(x1f,y1f,s=3,c=anp1,marker="x",vmin = 0.7, vmax=0.78)
im = ax1.scatter(x1r,y1r,s=3,c=anp1,marker="x",vmin = 0.7, vmax=0.78)
ax1.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax1.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax1.set_xlim((-10,10))
ax1.set_ylim((-20,20))
ax1.set_ylabel('Y(cm)', fontsize=24)
ax1.set_title('Channel 1', fontsize=26)
ax1.tick_params(axis='both', which='major',labelsize=20)
ax1.tick_params(axis='both', which='minor',labelsize=20)


im = ax2.scatter(x2f,y2f,s=3,c=anp2,marker="x",vmin = 0.7, vmax=0.78)
im = ax2.scatter(x2r,y2r,s=3,c=anp2,marker="x",vmin = 0.7, vmax=0.78)
ax2.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax2.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax2.set_xlim((-10,10))
ax2.set_ylim((-20,20))
ax2.set_title('Channel 2', fontsize=26)
ax2.tick_params(axis='both', which='major',labelsize=20)
ax2.tick_params(axis='both', which='minor',labelsize=20)

im = ax3.scatter(x3f,y3f,s=3,c=anp3,marker="x",vmin = 0.7, vmax=0.78)
im = ax3.scatter(x3r,y3r,s=3,c=anp3,marker="x",vmin = 0.7, vmax=0.78)
ax3.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax3.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax3.set_xlim((-10,10))
ax3.set_ylim((-20,20))
ax3.set_xlabel('X(cm)', fontsize=24)
ax3.set_ylabel('Y(cm)', fontsize=24)
ax3.set_title('Channel 3', fontsize=26)
ax3.tick_params(axis='both', which='major',labelsize=20)
ax3.tick_params(axis='both', which='minor',labelsize=20)

im = ax4.scatter(x4f,y4f,s=3,c=anp4,marker="x",vmin = 0.7, vmax=0.78)
im = ax4.scatter(x4r,y4r,s=3,c=anp4,marker="x",vmin = 0.7, vmax=0.78)
ax4.plot([2.3,2.3,6.3,6.3,2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax4.plot([-2.3,-2.3,-6.3,-6.3,-2.3],[-15.5,15.5,15.5,-15.5,-15.5], 'k-', lw=3)
ax4.set_xlim((-10,10))
ax4.set_ylim((-20,20))
ax4.set_xlabel('X(cm)', fontsize=24)
ax4.set_title('Channel 4', fontsize=26)
ax4.tick_params(axis='both', which='major',labelsize=20)
ax4.tick_params(axis='both', which='minor',labelsize=20)

plt.suptitle('Left Arm Channel-wise Events',fontsize=26)
cax = fig.add_axes([0.92,0.1,0.03,0.8])
thebar = plt.colorbar(im,cax=cax)
thebar.set_label('Analyzing Power', size=24)
plt.show()
