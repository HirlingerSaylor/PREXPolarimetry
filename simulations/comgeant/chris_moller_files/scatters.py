import pylab
import numpy as np
import matplotlib.pyplot as plt

#import data

thetc = np.loadtxt('run_thetc.dat').T
thetc = thetc*180/3.14159

phic = np.loadtxt('run_phic.dat').T
phic = phic*180/3.14159

thetl = np.loadtxt('run_thetl.dat').T
thetl = thetl*180/3.14159

phil = np.loadtxt('run_phil.dat').T
phil = phil*180/3.14159

thetr = np.loadtxt('run_thetr.dat').T
thetr = thetr*180/3.14159

phir = np.loadtxt('run_phir.dat').T
phir = phir*180/3.14159



#do some math

philp=-1*phil
phirp=-1*phir
thetlp=thetl-2*(thetl-90)
thetrp=thetr-2*(thetr-90)



#throw down some plots

plt.plot(phil, thetl, 'bo', markersize=10, alpha=0.1)
plt.plot(phir, thetr, 'yo', markersize=10, alpha=0.1)


plt.show()
#pylab.savefig('scatters.png')

