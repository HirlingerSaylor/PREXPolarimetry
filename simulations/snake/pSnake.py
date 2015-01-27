#!/usr/bin/python
# Plot polarimeter layout straight from directive file, and overlay
# trajectories from pts file.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat

# Rotation function for rotation of np vector about Z axis
def rotX(vec,cen,ang):
	r = np.array([[np.cos(ang),-np.sin(ang),0.],
								[np.sin(ang),np.cos(ang),0.],
								[0.,0.,1.]])
	return np.dot(r,vec-cen)+cen	
	
# Convert SNAKE box parameters to corner points
def boxPts(a):
	newA = np.zeros((8,3))
	newA[0] = np.array([a[0,0]+a[2,0],a[0,1]+a[2,1],a[0,2]+a[2,2]])
	newA[1] = np.array([a[0,0]+a[3,0],a[0,1]+a[2,1],a[0,2]+a[2,2]])
	newA[2] = np.array([a[0,0]+a[2,0],a[0,1]+a[3,1],a[0,2]+a[2,2]])
	newA[3] = np.array([a[0,0]+a[3,0],a[0,1]+a[3,1],a[0,2]+a[2,2]])
	newA[4] = np.array([a[0,0]+a[2,0],a[0,1]+a[2,1],a[0,2]+a[3,2]])
	newA[5] = np.array([a[0,0]+a[3,0],a[0,1]+a[2,1],a[0,2]+a[3,2]])
	newA[6] = np.array([a[0,0]+a[2,0],a[0,1]+a[3,1],a[0,2]+a[3,2]])
	newA[7] = np.array([a[0,0]+a[3,0],a[0,1]+a[3,1],a[0,2]+a[3,2]])
	if a[1,0] != 0:
		for i in xrange(len(newA)):
			newA[i] = rotX(newA[i],a[0],a[1,0] *np.pi /180.)
	elif a[1,1] != 0 or a[1,2] != 0:
		print('Rotations about X and Y not yet implemented.')
	return newA

# box functions very specific to Hall-A Moller Polarimeter
def boxZ(b):
	newB = np.array([ [b[0,1],b[0,0]],
										[b[1,1],b[1,0]],
										[b[3,1],b[3,0]],
										[b[2,1],b[2,0]] ])
	return newB

def boxX(b):
	newB = np.array([ [b[0,1],b[0,2]],
										[b[4,1],b[4,2]],
										[b[5,1],b[5,2]],
										[b[1,1],b[1,2]],
										[b[3,1],b[3,2]],
										[b[7,1],b[7,2]],
										[b[5,1],b[5,2]],
										[b[1,1],b[1,2]] ])
	return newB

# Box points for open endplanes
def oep(ep,reg):
	crns = np.array([ [reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10] ])
	if reg[1,0] != 0:
		for i in xrange(len(crns)):
			crns[i] = rotX(crns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((2,4,2))
	pts[0] = np.array([ [crns[0,1],crns[0,0]],
											[crns[1,1],crns[1,0]],
											[crns[3,1],crns[3,0]],
											[crns[2,1],crns[2,0]] ])
	pts[1] = np.array([ [crns[0,1],crns[0,2]],
											[crns[1,1],crns[1,2]],
											[crns[3,1],crns[3,2]],
											[crns[2,1],crns[2,2]] ])
	return pts

# Box points for blocking aperture
def bep(ep,reg):
	crns = np.array([ [reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +ep[5]],
										[reg[0,0] +ep[3],reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +ep[3],reg[0,1] +ep[1],reg[0,2] +ep[5]] ])
	if reg[1,0] != 0:
		for i in xrange(len(crns)):
			crns[i] = rotX(crns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((2,4,2))
	pts[0] = np.array([ [crns[0,1],crns[0,0]],
											[crns[1,1],crns[1,0]],
											[crns[3,1],crns[3,0]],
											[crns[2,1],crns[2,0]] ])
	pts[1] = np.array([ [crns[0,1],crns[0,2]],
											[crns[1,1],crns[1,2]],
											[crns[3,1],crns[3,2]],
											[crns[2,1],crns[2,2]] ])
	return pts

# Box points for circular aperture
def cep(ep,reg):
	iCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] -ep[3]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] -ep[3]],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] -ep[2],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] -ep[2],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	oCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	if reg[1,0] != 0:
		for i in xrange(8):
			iCrns[i] = rotX(iCrns[i],reg[0],reg[1,0] *np.pi /180.)
			if i < 4:
				oCrns[i] = rotX(oCrns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((2,2,4,2))
	pts[0,0] = np.array([ [oCrns[0,1],oCrns[0,0]],
												[oCrns[1,1],oCrns[1,0]],
												[iCrns[6,1],iCrns[6,0]],
												[iCrns[4,1],iCrns[4,0]] ])
	pts[0,1] = np.array([ [oCrns[2,1],oCrns[2,2]],
												[oCrns[0,1],oCrns[0,2]],
												[iCrns[0,1],iCrns[0,2]],
												[iCrns[2,1],iCrns[2,2]] ])
	pts[1,0] = np.array([ [oCrns[3,1],oCrns[3,0]],
												[oCrns[2,1],oCrns[2,0]],
												[iCrns[5,1],iCrns[5,0]],
												[iCrns[7,1],iCrns[7,0]] ])
	pts[1,1] = np.array([ [oCrns[3,1],oCrns[3,2]],
												[oCrns[1,1],oCrns[1,2]],
												[iCrns[1,1],iCrns[1,2]],
												[iCrns[3,1],iCrns[3,2]] ])
	return pts

# Box points for trapezoidal apertures
def tep(ep,reg):
	iCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[2]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +ep[2]],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] +ep[6],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[4],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[6],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +ep[4],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	oCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	if reg[1,0] != 0:
		for i in xrange(8):
			iCrns[i] = rotX(iCrns[i],reg[0],reg[1,0] *np.pi /180.)
			if i < 4:
				oCrns[i] = rotX(oCrns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((2,2,4,2))
	pts[0,0] = np.array([ [oCrns[0,1],oCrns[0,0]],
												[oCrns[1,1],oCrns[1,0]],
												[iCrns[6,1],iCrns[6,0]],
												[iCrns[4,1],iCrns[4,0]] ])
	pts[0,1] = np.array([ [oCrns[2,1],oCrns[2,2]],
												[oCrns[0,1],oCrns[0,2]],
												[iCrns[0,1],iCrns[0,2]],
												[iCrns[2,1],iCrns[2,2]] ])
	pts[1,0] = np.array([ [oCrns[3,1],oCrns[3,0]],
												[oCrns[2,1],oCrns[2,0]],
												[iCrns[5,1],iCrns[5,0]],
												[iCrns[7,1],iCrns[7,0]] ])
	pts[1,1] = np.array([ [oCrns[3,1],oCrns[3,2]],
												[oCrns[1,1],oCrns[1,2]],
												[iCrns[1,1],iCrns[1,2]],
												[iCrns[3,1],iCrns[3,2]] ])
	return pts

# Box points for vertically symmetric dual apertures
def dep(ep,reg): 
	oCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	iCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[5]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[6]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[3]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[5]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[6]],
										[reg[0,0] -ep[2],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] -ep[2],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	if reg[1,0] != 0:
		for i in xrange(12):
			iCrns[i] = rotX(iCrns[i],reg[0],reg[1,0] *np.pi /180.)
			if i < 4:
				oCrns[i] = rotX(oCrns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((5,4,2))
	pts[0] = np.array([ [oCrns[2,1],oCrns[2,0]],
											[oCrns[3,1],oCrns[3,0]],
											[iCrns[11,1],iCrns[11,0]],
											[iCrns[9,1],iCrns[9,0]] ])
	pts[1] = np.array([ [oCrns[0,1],oCrns[0,0]],
											[oCrns[1,1],oCrns[1,0]],
											[iCrns[10,1],iCrns[10,0]],
											[iCrns[8,1],iCrns[8,0]] ])
	pts[2] = np.array([ [oCrns[0,1],oCrns[0,2]],
											[oCrns[2,1],oCrns[2,2]],
											[iCrns[4,1],iCrns[4,2]],
											[iCrns[0,1],iCrns[0,2]] ])
	pts[3] = np.array([ [iCrns[5,1],iCrns[5,2]],
											[iCrns[1,1],iCrns[1,2]],
											[iCrns[2,1],iCrns[2,2]],
											[iCrns[6,1],iCrns[6,2]] ])
	pts[4] = np.array([ [oCrns[1,1],oCrns[1,2]],
											[oCrns[3,1],oCrns[3,2]],
											[iCrns[7,1],iCrns[7,2]],
											[iCrns[3,1],iCrns[3,2]] ])
	return pts
											
# Box points for horizontally symmetric dual apertures
def dep2(ep,reg):
	oCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +reg[3,0] +10,reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	iCrns = np.array([[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] -ep[5]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] -ep[4]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +reg[2,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[5]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] -ep[5]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] -ep[4]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[4]],
										[reg[0,0] +reg[3,0] -10,reg[0,1] +ep[1],reg[0,2] +ep[5]],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[3],reg[0,1] +ep[1],reg[0,2] +reg[2,2] -10],
										[reg[0,0] +ep[2],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10],
										[reg[0,0] +ep[3],reg[0,1] +ep[1],reg[0,2] +reg[3,2] +10]])
	if reg[1,0] != 0:
		for i in xrange(12):
			iCrns[i] = rotX(iCrns[i],reg[0],reg[1,0] *np.pi /180.)
			if i < 4:
				oCrns[i] = rotX(oCrns[i],reg[0],reg[1,0] *np.pi /180.)
	pts = np.zeros((5,4,2))
	pts[0] = np.array([ [oCrns[2,1],oCrns[2,0]],
											[oCrns[3,1],oCrns[3,0]],
											[iCrns[11,1],iCrns[11,0]],
											[iCrns[9,1],iCrns[9,0]] ])
	pts[1] = np.array([ [oCrns[0,1],oCrns[0,0]],
											[oCrns[1,1],oCrns[1,0]],
											[iCrns[10,1],iCrns[10,0]],
											[iCrns[8,1],iCrns[8,0]] ])
	pts[2] = np.array([ [oCrns[0,1],oCrns[0,2]],
											[oCrns[2,1],oCrns[2,2]],
											[iCrns[4,1],iCrns[4,2]],
											[iCrns[0,1],iCrns[0,2]] ])
	pts[3] = np.array([ [iCrns[5,1],iCrns[5,2]],
											[iCrns[1,1],iCrns[1,2]],
											[iCrns[2,1],iCrns[2,2]],
											[iCrns[6,1],iCrns[6,2]] ])
	pts[4] = np.array([ [oCrns[1,1],oCrns[1,2]],
											[oCrns[3,1],oCrns[3,2]],
											[iCrns[7,1],iCrns[7,2]],
											[iCrns[3,1],iCrns[3,2]] ])
	return pts
												
										
			
# Read and interpret directive file
lines = [line.strip() for line in open('d.dat')]

region = 0
names = []
regions = []
eps = []
boxes = []

for i in xrange(4,len(lines)):
	if len(lines[i])==40 and lines[i][39]=='*':
		tBox = np.zeros((4,3))
		tBox[0] = map(float,lines[i+3].split(','))
		tBox[1] = map(float,lines[i+4].split(','))
		tBox[2] = map(float,lines[i+5].split(','))
		tBox[3] = map(float,lines[i+6].split(','))
		names.append(lines[i+1])
		regions.append(tBox)
		boxes.append(boxPts(tBox))
	if lines[i][0]=='y':
		eps.append(np.array([region,
										float(lines[i].split(',')[1]),
										float(lines[i].split(',')[4]),
										float(lines[i].split(',')[5]),
										float(lines[i].split(',')[6]),
										float(lines[i].split(',')[7]),
										float(lines[i].split(',')[8]),
										float(lines[i].split(',')[9])]))
	if lines[i][0:3]=='eol':
		region += 1


# Sort endplane apertures
oEps = []
bEps = []
cEps = []
tEps = []
dzEps = []
dxEps = []
for ep in eps:

	# Open Endplane
	if all(ep[2:8]==0):
		oEps.append(oep(ep,regions[int(ep[0])]))

	# Blocking Aperture
	elif ep[7]==999:
		bEps.append(bep(ep,regions[int(ep[0])]))

	# Dual Aperture
	elif ep[7]==888:
		dzEps.append(dep(ep,regions[int(ep[0])])[0])
		dzEps.append(dep(ep,regions[int(ep[0])])[1])
		dxEps.append(dep(ep,regions[int(ep[0])])[2])
		dxEps.append(dep(ep,regions[int(ep[0])])[3])
		dxEps.append(dep(ep,regions[int(ep[0])])[4])
	elif ep[7]==777:
		dzEps.append(dep2(ep,regions[int(ep[0])])[0])
		dzEps.append(dep2(ep,regions[int(ep[0])])[1])
		dxEps.append(dep2(ep,regions[int(ep[0])])[2])
		dxEps.append(dep2(ep,regions[int(ep[0])])[3])
		dxEps.append(dep2(ep,regions[int(ep[0])])[4])

	# Circular Aperture
	elif ep[7]==0:
		cEps.append(cep(ep,regions[int(ep[0])])[0])
		cEps.append(cep(ep,regions[int(ep[0])])[1])

	# Trapezoidal Aperture
	else:
		tEps.append(tep(ep,regions[int(ep[0])])[0])
		tEps.append(tep(ep,regions[int(ep[0])])[1])


# Sort trajectory points, units->mm
fDat = np.loadtxt('dFPts.txt')
fse = np.zeros((max(fDat[:,0]),max(fDat[:,1]),3))
for e in fDat:
  fse[e[0]-1,e[1]-1] = e[2:5]
for i in xrange(len(fse)):
  for j in xrange(1,len(fse[0])):
    if all(fse[i,j]==0.0) and j!=0:
      fse[i,j] = fse[i,j-1]
rDat = np.loadtxt('dRPts.txt')
rse = np.zeros((max(rDat[:,0]),max(rDat[:,1]),3))
for e in rDat:
  rse[e[0]-1,e[1]-1] = e[2:5]
for i in xrange(len(rse)):
  for j in xrange(1,len(rse[0])):
    if all(rse[i,j]==0.0) and j!=0:
      rse[i,j] = rse[i,j-1]

# Generate figures
fig = plt.figure()

ax1 = fig.add_subplot(212)
for i in xrange(len(boxes)):
	ax1.add_patch(pat.Polygon(boxZ(boxes[i]),fill=False))
for i in xrange(len(oEps)):
	ax1.add_patch(pat.Polygon(oEps[i][0],color='#00FF00',lw=2))
for i in xrange(len(bEps)):
	ax1.add_patch(pat.Polygon(bEps[i][0],color='r',lw=2))
for i in xrange(len(cEps)):
	ax1.add_patch(pat.Polygon(cEps[i][0],color='r',lw=2))
for i in xrange(len(tEps)):
	ax1.add_patch(pat.Polygon(tEps[i][0],color='r',lw=2))
for i in xrange(len(dzEps)):
	ax1.add_patch(pat.Polygon(dzEps[i],color='r',lw=2))
for traj in fse:
	ax1.plot(traj[:,1],traj[:,0],'b-',alpha=1)
for traj in rse:
	ax1.plot(traj[:,1],traj[:,0],'b-',alpha=1)
ax1.set_xlim((-900,7000))
ax1.set_ylim((-700,400))
ax1.set_xlabel('y (mm)',fontsize=24)
ax1.set_ylabel('x (mm)',fontsize=24)

ax2 = fig.add_subplot(211, sharex=ax1)
for i in xrange(len(boxes)):
	ax2.add_patch(pat.Polygon(boxX(boxes[i]),fill=False))
for i in xrange(len(oEps)):
	ax2.add_patch(pat.Polygon(oEps[i][1],color='#00FF00',lw=2))
for i in xrange(len(bEps)):
	ax2.add_patch(pat.Polygon(bEps[i][1],color='r',lw=2))
for i in xrange(len(cEps)):
	ax2.add_patch(pat.Polygon(cEps[i][1],color='r',lw=2))
for i in xrange(len(tEps)):
	ax2.add_patch(pat.Polygon(tEps[i][1],color='r',lw=2))
for i in xrange(len(dxEps)):
	ax2.add_patch(pat.Polygon(dxEps[i],color='r',lw=2))
for traj in fse:
	ax2.plot(traj[:,1],traj[:,2],'b-')
for traj in rse:
	ax2.plot(traj[:,1],traj[:,2],'b-')
ax2.set_ylim((-300,300))
ax2.set_ylabel('z (mm)',fontsize=24)

plt.title('Hall-A Moller Polarimeter with SNAKE',fontsize=26)
plt.xlim((-900,7000))
plt.show()
