import pylab
import numpy as np
import matplotlib.pyplot as plt

#import data

thetc = np.loadtxt('run_thetc.dat').T
anpower= (np.sin(thetc)**2*(7+np.cos(thetc)**2))/((3+np.cos(thetc)**2)**2)

phic = np.loadtxt('run_phic.dat').T

thetl = np.loadtxt('run_thetl.dat').T

phil = np.loadtxt('run_phil.dat').T

thetr = np.loadtxt('run_thetr.dat').T

phir = np.loadtxt('run_phir.dat').T

apxlc = np.loadtxt('run_xlc.dat').T

apxrc = np.loadtxt('run_xrc.dat').T

apylc = np.loadtxt('run_ylc.dat').T

apyrc = np.loadtxt('run_yrc.dat').T

siglc = np.loadtxt('run_siglc.dat').T

sigrc = np.loadtxt('run_sigrc.dat').T

#do some math

counts=len(thetc)
ancut= []
anpower= (np.sin(thetc)**2*(7+np.cos(thetc)**2))/((3+np.cos(thetc)**2)**2)

for i in range(1,len(sigrc)):
	if sigrc[i]>5250 and sigrc[i]<5750:
		ancut.append(anpower[i])


thetc = thetc*180/3.14159
phic = phic*180/3.14159
thetl = thetl*180/3.14159
phil = phil*180/3.14159
thetr = thetr*180/3.14159
phir = phir*180/3.14159


philp=-1*phil
phirp=-1*phir
thetlp=thetl-2*(thetl-90)
thetrp=thetr-2*(thetr-90)


#throw down some plots

#thet phi scatter
plt.figure(1)
plt.plot(phil, thetl, 'bx', markersize=5, label="Left Arm")
plt.plot(phir, thetr, 'rx', markersize=5, label="Right Arm")
plt.plot(phic, thetc, 'kx', markersize=5, label="Coincidence")
plt.axis([-40, 40, 70, 110])
plt.legend()

#thet phi rotated
plt.figure(2)
plt.plot(phil, thetl, 'bx', markersize=5, label="Left Arm")
plt.plot(phirp, thetrp, 'rx', markersize=5, label="Right Arm")
plt.axis([-40, 40, 70, 110])
plt.legend()

#xy
plt.figure(3)

ax=plt.subplot(121)
ax.scatter(apxlc, apylc, s=5, c=anpower, marker='x')
ax.axis([-6.6,-4.7,-42.0,-30.0])
ax.set_title('Left Arm')


ax=plt.subplot(122)
im = ax.scatter(apxrc, apyrc, s=5, c=anpower, marker='x')
ax.axis([4.7,6.6,-42.0,-30.0])
ax.set_title('Right Arm')
thebar = plt.colorbar(im)

#energy anpower scatter
plt.figure(4)

hist, xedges, yedges=np.histogram2d(anpower, siglc, bins=60, range=[[.76,.78],[4000,6500]])
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
plt.imshow(hist, extent=extent, aspect=100000, interpolation='nearest', origin='lower')
plt.colorbar()

#anpower hist
plt.figure(5)

plt.hist(anpower, 100)
plt.xlim(0.75, 0.78)
plt.ylim(0,7000)
plt.title("%s counts, avg = %s"%(len(anpower), np.mean(anpower)))
#anpower cut hist
plt.figure(6)

plt.hist(ancut, 100)
plt.xlim(0.75, 0.78)
plt.ylim(0,7000)
plt.title("%s counts, avg = %s"%(len(ancut), np.mean(ancut)))

#theta hist
plt.figure(7)
plt.hist(thetc, 100)

#theta hist
plt.figure(8)
plt.hist(phic, 100)

plt.show()

#pylab.savefig('scatters.png')

