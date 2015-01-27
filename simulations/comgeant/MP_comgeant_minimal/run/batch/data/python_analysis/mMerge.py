#!/usr/bin/python
# Merges files dumped from comgeant ntuple

import numpy as np

# Merge angle data from initial events into single file
thet = np.loadtxt('thetcm_i.dat')
phi = np.loadtxt('phicm_i.dat')
initAng = open('dInitialAngles.txt','w')
for i in xrange(len(thet)):
	initAng.write(str(thet[i])+' '+str(phi[i])+'\n')
initAng.close()

# Merge angle and position data from coincident events into single file
fileList = ['thetcm_f.dat',
						'phicm_f.dat',
						'anp.dat',
						'x_l.dat',
						'x_r.dat',
						'y_l.dat',
						'y_r.dat',
						'z_l.dat',
						'z_r.dat']
data = []
for fil in fileList:
	data.append(np.loadtxt(fil))
coinData = open('dCoincidentData.txt','w')
for i in xrange(len(data[0])):
	dataList = ''
	for j in xrange(len(fileList)):
		dataList += str(data[j][i])
		if j < len(fileList) -1:
			dataList += ' '
		else:
			dataList += '\n'
	coinData.write(dataList)
coinData.close()

# Merge energy by detector channel files into sigle file
ene = []
for i in xrange(8):
	ene.append(np.loadtxt('ene_'+str(i+1)+'.dat'))
eneCha = open('dEnergyByChannels.txt','w')
for i in xrange(len(ene[0])):
	eneList = ''
	for j in xrange(8):
		eneList += str(ene[j][i])
		if j < 7:
			eneList += ' '
		else:
			eneList += '\n'
	eneCha.write(eneList)
eneCha.close()
