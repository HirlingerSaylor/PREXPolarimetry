#!/usr/bin/python

from subprocess import call
import numpy as np


# Forward scattered electron cm->lab
def fse(thet,ene):
	the = thet *np.pi /180.
	ec = np.sqrt(ene*0.000511/2.)
	v = np.sqrt( (ene-0.000511)/(ene+0.000511) )
	g = np.sqrt( (ene+0.000511)/(2.*0.000511) )
	en = ec*g*(1+v*np.cos(the))
	ang = np.arctan(np.sin(the) /(g*(v+np.cos(the))))
	return np.array([en,ang])

# Reverse scattered electron cm->lab
def rse(thet,ene):
	the = thet *np.pi /180.
	ec = np.sqrt(ene*0.000511/2.)
	v = np.sqrt( (ene-0.000511)/(ene+0.000511) )
	g = np.sqrt( (ene+0.000511)/(2.*0.000511) )
	en = ec*g*(1-v*np.cos(the))
	ang = np.arctan(np.sin(the) /(g*(np.cos(the)-v)))
	return np.array([en,ang])

# Lab angles to SNAKE angles
def angCon(lThet,phi):
	ph = phi *np.pi /180.
	return np.array([lThet *np.sin(ph),-lThet *np.cos(ph)])

# Moller scattering cross section
#  with longitudinal beam polarization p (i.e. p_targ=1)
def cs(thet,p,ene):
	th = thet*np.pi /180.
	v = np.sqrt( (ene-0.000511)/(ene+0.000511) )
	g = np.sqrt( (ene+0.000511)/(2.*0.000511) )
	lcs = ( (v*np.cos(th)+1.)*( (1.+np.cos(th))*(3.+np.cos(th)**2.) )**2. 
			/( g**2. *np.sin(th)**3. *(v+np.cos(th))**3.
			*(np.sin(th)**2. /(g**2. *(v+np.cos(th))**2.) +1.)**(3./2.)))
	asy = (1-p*np.sin(th)**2. *(7.+np.cos(th)**2.)
			/( (3.+np.cos(th)**2.)**2. ))
	return lcs*asy


# Spectrometer Settings
#spec = np.array([0.,0.180542,0.0407379,0.0,0.0445758,0.18])#PREX settings
spec = np.array([0.,0.0327421188133,0.179795628736,0.0,0.0629126565596,0.18*2.2/1.063])#CREX settings

# Magnets
#  (center of box, box 1/2 widths)
mag = np.array([[[0.,0.,0.],[200.,800.,200.]],
                [[0.,682.9,0.],[50.8,182.9,50.8]],
                [[0.,1335.6,0.],[50.8,223.0,50.8]],
                [[0.,2021.8,0.],[50.8,183.7,50.8]],
                [[0.,2676.9,0.],[50.8,183.7,50.8]],
                [[0.,4165.,0.],[300.,822.5,80.]]])

# Detector
#  (center of detector face, rotation)
det = np.array([[449.,6926.,0.],[-7.,0.,0.]])


## Create Directive File ##
d = open('d.dat','w')

# Particle and Integration Parameters
d.write('20.,40.,20.\n'                   # MinStep, MaxStep, StepSize
        +'false,false\n'                  # Spin RayTracing, Synchrotron Radiation
        +'-1.,0.000511\n'                 # Particle Charge, Mass
        +'1.,1.7588e11,1.,0.000511\n')    # ?(magneton), q/m (rad s-1 T-1), - , m(GeV)

# Target Field Box
d.write('Helmholtz*******************************\n'    # Name
        +'HELMHOLT\n'                                   # Label
        +'absolute\n'                                   # Box Placement
        +str(mag[0,0,0])+','
          +str(mag[0,0,1])+','
          +str(mag[0,0,2])+'\n'                         # Local Origin in Global Coords
        +'0.,0.,0.\n'                                   # Local Coord Rotation
        +str(-mag[0,1,0])+','
					+str(-mag[0,1,1])+','
					+str(-mag[0,1,2])+'\n'                        # Box Minimums
        +str(mag[0,1,0])+','
					+str(mag[0,1,1])+','
					+str(mag[0,1,2])+'\n'                         # Box Maximums
        +'3\n'                                          # Integration 3->Runge-Kutta
        +'2\n'                                          # Field Type 2->analytic
        +'19\n'                                         # Magnet Type 19->HH-Coil
        +'none\n'                                       # Field Map or Dipole Input
        +'r\n'                                          # r->rectangular box(?)
        +str(spec[0])+'\n'                              # Field Multiplier
        +'2\n'                                          # Additional Data (see manual)
        +'350000.\n'                                    #  coil radius?
        +'0.09\n'                                       #  current?
        +'list\n'                                       # Endplane Definition Type
        +'y,0.,0.,none,0.,0.,0.,0.,0.,0.\n'             # Helmholtz Endplanes

				# Q1 Beginning Fringe Aperture
				+'y,'+str(mag[1,0,1]-mag[1,1,1]-100.)+',0.,none,47.8,47.8,0.,0.,0.,0.\n'
        +'eol\n')


# Quad 1
#  Beginning
d.write('QUAD_01B********************************\n'
				+'QUAD_01B\nabsolute\n'
				+str(mag[1,0,0])+','
					+str(mag[1,0,1]-mag[1,1,1])+','
					+str(mag[1,0,2])+'\n'
				+'0.,0.,0.\n'
				+str(-mag[1,1,0])+','
					+'-110.,'
					+str(-mag[1,1,2])+'\n'
				+str(mag[1,1,0])+','
					+str(mag[1,1,1]+10.)+','
					+str(mag[1,1,2])+'\n'
				+'3\n2\n7\nnone\nr\n1.\n6\n'
				+str(spec[1])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+'list\n'

				# Q1 Beginning EFB Aperture
				+'y,0.,0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q1 Midpoint
				+'y,'+str(mag[1,1,1])+',0.,none,47.8,47.8,0.,0.,0.,0.\n'
				+'eol\n')

# Q1/Q2 Overlap
d.write('QUAD0102********************************\n'
				+'QUAD0102\nabsolute\n'
				+'0.,'+str(mag[1,0,1]+mag[1,1,1])+',0.\n'
				+'0.,0.,0.\n'
				+str(-mag[1,1,0])+','
					+str(-mag[1,1,1]-10.)+','
					+str(-mag[1,1,2])+'\n'
				+str(mag[1,1,0])+','
					+str(mag[2,0,1]-mag[1,0,1]-mag[1,1,1]+10.)+','
					+str(mag[1,1,2])+'\n'
				+'3\n2\n10\nnone\nr\n1.\n13\n'
				+str(spec[1])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(spec[2])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(mag[2,0,1]-mag[2,1,1]-mag[1,0,1]-mag[1,1,1])+'\n'
				+'list\n'

				# Q1 End EFB Aperture
				+'y,0.,0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q2 Begin EFB Aperture
				+'y,'+str(mag[2,0,1]-mag[2,1,1]-mag[1,0,1]-mag[1,1,1])
					+',0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q2 Midpoint
				+'y,'+str(mag[2,0,1]-mag[1,0,1]-mag[1,1,1])
					+',0.,none,0.,0.,0.,0.,0.,0.\n'
				+'eol\n')

# Q2/Q3 Overlap
d.write('QUAD0203********************************\n'
				+'QUAD0203\nabsolute\n'
				+'0.,'+str(mag[2,0,1]+mag[2,1,1])+',0.\n'
				+'0.,0.,0.\n'
				+str(-mag[2,1,0])+','
					+str(-mag[2,1,1]-10.)+','
					+str(-mag[2,1,2])+'\n'
				+str(mag[2,1,0])+','
					+str(mag[3,0,1]-mag[2,0,1]-mag[2,1,1]+10.)+','
					+str(mag[2,1,2])+'\n'
				+'3\n2\n10\nnone\nr\n1.\n13\n'
				+str(spec[2])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(spec[3])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(mag[3,0,1]-mag[3,1,1]-mag[2,0,1]-mag[2,1,1])+'\n'
				+'list\n'

				# Q2 End EFB Aperture
				+'y,0.,0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q3 Begin EFB Aperture
				+'y,'+str(mag[3,0,1]-mag[3,1,1]-mag[2,0,1]-mag[2,1,1])
					+',0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q3 Midpoint
				+'y,'+str(mag[3,0,1]-mag[2,0,1]-mag[2,1,1])
					+',0.,none,0.,0.,0.,0.,0.,0.\n'
				+'eol\n')

# Q3/Q4 Overlap
d.write('QUAD0304********************************\n'
				+'QUAD0304\nabsolute\n'
				+'0.,'+str(mag[3,0,1]+mag[3,1,1])+',0.\n'
				+'0.,0.,0.\n'
				+str(-mag[3,1,0])+','
					+str(-mag[3,1,1]-10.)+','
					+str(-mag[3,1,2])+'\n'
				+str(mag[3,1,0])+','
					+str(mag[4,0,1]-mag[3,0,1]-mag[3,1,1]+10.)+','
					+str(mag[3,1,2])+'\n'
				+'3\n2\n10\nnone\nr\n1.\n13\n'
				+str(spec[3])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(spec[4])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+str(mag[4,0,1]-mag[4,1,1]-mag[3,0,1]-mag[3,1,1])+'\n'
				+'list\n'

				# Q3 End EFB Aperture
				+'y,0.,0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q4 Begin EFB Aperture
				+'y,'+str(mag[4,0,1]-mag[4,1,1]-mag[3,0,1]-mag[3,1,1])
					+',0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q4 Midpoint
				+'y,'+str(mag[4,0,1]-mag[3,0,1]-mag[3,1,1])
					+',0.,none,0.,0.,0.,0.,0.,0.\n'
				+'eol\n')

# Quad 4
#  End
d.write('QUAD_04E********************************\n'
				+'QUAD_04E\nabsolute\n'
				+str(mag[4,0,0])+','
					+str(mag[4,0,1]+mag[4,1,1])+','
					+str(mag[4,0,2])+'\n'
				+'0.,0.,0.\n'
				+str(-mag[4,1,0])+','
					+str(-mag[4,1,1]-10.)+','
					+str(-mag[4,1,2])+'\n'
				+str(mag[4,1,0])+','
					+str(110.)+','
					+str(mag[4,1,2])+'\n'
				+'3\n2\n9\nnone\nr\n1.\n6\n'
				+str(spec[4])+'\n'
				+'0.\n0.\n0.\n0.\n50.8\n'
				+'list\n'

				# Q4 End EFB Aperture
				+'y,0.,0.,none,47.8,47.8,0.,0.,0.,0.\n'

				# Q4 End Fringe Aperture
				+'y,100.,0.,none,47.8,47.8,0.,0.,0.,0.\n'
				+'eol\n')


# Drift 4
d.write('DRIFT_04********************************\n'
				+'DRIFT_04\nabsolute\n'
				+'0.,3100.,0.\n'
				+'0.,0.,0.\n'
				+'-47.8,-200.,-47.8\n'
				+'47.8,150.,47.8\n'
				+'1\n0\n0\nnone\nr\n0.\n0\nlist\n'

				# Dipole Box Entrance Aperture
				#  positive and negative apertures respectively
				+'y,74.,0.,none,47.8,47.8,0.,0.,0.,0.\n'
				+'y,97.,0.,none,-47.8,47.8,-30.,30.,0.,999\n'

				# Dipole Box Entrance Collimator (DCOL from GEANT)
				+'y,109.,0.,none,15.,-50.,-30.,30.,50.,888\n'
				+'y,112.5,0.,none,0.,0.,0.,0.,0.,0.\n'
				+'eol\n')


# Dipole
d.write('DIPOLE01********************************\n'
                                +'DIPOLE01\nabsolute\n'
				+'-136.47,4987.5,0.\n'
				+'0.,0.,0.\n'
				+'-163.53,-1785.,-80.\n'
				+'436.47,140.,80.\n'
				+'3\n2\n18\ndip.dat\nr\n'
                                +'2.2/1.063\n' #+'1.\n' for PREX, but scale by 2.2/1.063 for CREX
                                +'0\nlist\n'
                                +'y,-1718.5,0.,none,121.47,151.47,30.,50.,0.,777\n'
				+'y,130.,0.,none,0.,0.,0.,0.,0.,0.\n'
				+'eol\n')


# Drift 5
d.write('DRIFT_05********************************\n'
				+'DRIFT_05\nabsolute\n'
				+'-160.,5150.,0.\n'
				+'0.,0.,0.\n'
				+'-90.,-50.,-60.\n'
				+'230.,50.,60.\n'
				+'1\n0\n0\nnone\nr\n0.\n0\nlist\n'

				# Dipole Box Exit Aperture
				+'y,-29.,0.,none,-80.,220.,-30.,30.,0.,999\n'
				+'y,-6.,0.,none,80.,-52.95,-29.5,29.5,52.95,888\n'
				+'y,6.8,0.,none,80.,-53.1,-29.5,29.5,53.1,888\n'
				+'eol\n')


# Detector
d.write('DETECTOR********************************\n'
				+'DETECTOR\nabsolute\n'
				+'-450.39,6866.42,0.\n'
				+'7.,0.,0.\n'
				+'-180.,-1800.,-70.\n'
				+'180.,10.,70.\n'
				+'1\n0\n0\nnone\nr\n0.\n0\nlist\n'

				# Detector Aperture (HOD1 in GEANT)
				+'y,-6.5,0.,none,155.,-63.,-23.,23.,63.,888\n'
				+'y,6.5,0.,none,155.,-63.,-23.,23.,63.,888\n'
				+'eol\n')


d.close()


## Create Dipole File ##
dip = open('dip.dat','w')

dip.write( '2.,4.,2.,8.,4.\n'
					+'0.,0.,15.,9982.3,-0.1631456798\n'
					+'9.485,0.,9.485\n'
					+'0.,0.,0.,0.\n'
					+'25.,-25.,-25.,25.\n'
					+'0.04725,2.2395,-0.9768,0.7288,-0.1299,0.0222\n'
					+'0.04725,2.2395,-0.9768,0.7288,-0.1299,0.0222\n'
					+'0.,0.,0.,0.,0.,0.\n'
					+'0.,0.,0.,0.\n'
					+'0.,0.,0.,0.,0.,0.,0.\n'
					+'0.,0.,0.,0.,0.,0.,0.\n' )

dip.close()


## Create Trajectory Files ##
t1 = open('t1.dat','w')
t2 = open('t2.dat','w')


# Random Thet/Phy traj generator
#  distributiom from scattering cross section
ang = open('analysis/dAngs.txt','w')
tLim = np.array([70.,110.])
pLim = np.array([-7.,7.])
xLim = 0.1
zLim = 0.1
eLim = 0.0002
beamenergy = 2.2 #2.2 GeV for CREX, 1.063 for PREX
pol = 0.
tAngs = np.zeros((7500,2))
t1.write('i\n'+str(len(tAngs))+'\n')
t2.write('i\n'+str(len(tAngs))+'\n')
i = 0
while i<len(tAngs):
	thet = np.random.uniform(tLim[0],tLim[1])
	e = np.random.uniform(beamenergy-eLim,beamenergy+eLim)
	prob = np.random.uniform(0.,1.)
	if prob <= cs(thet,pol,e)/cs(tLim[0],pol,e):
		phi = np.random.uniform(pLim[0],pLim[1])
		x = np.random.uniform(-xLim,xLim)
		z = np.random.uniform(-zLim,zLim)
		tAngs[i] = np.array([thet,phi])
		ang.write(str(thet)+' '+str(phi)+'\n')
		e1 = fse(thet,e)
		e2 = rse(thet,e)
		a1 = angCon(e1[1],phi)
		a2 = angCon(e2[1],phi)
		t1.write(str(x)+',0.,'+str(z)+','+str(a1[0])+','+str(a1[1])+','+str(e1[0])+'\n')
		t2.write(str(x)+',0.,'+str(z)+','+str(a2[0])+','+str(a2[1])+','+str(e2[0])+'\n')
		i += 1
		

t1.close()
t2.close()
ang.close()


## Run Snak ##
call('cp t1.dat t.dat; ./xSnake',shell=True)
call('cp data.txt analysis/dFse.txt; cp pts.txt dFPts.txt;'
			+'cp t2.dat t.dat; ./xSnake',shell=True)
call('cp data.txt analysis/dRse.txt; cp pts.txt dRPts.txt;'
			+'rm t.dat data.txt pts.txt',shell=True)
