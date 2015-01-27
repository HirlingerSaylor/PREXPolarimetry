Welcome to snake! - Ted (goofyboson@gmail.com) 6/14

Intro
	SNAKE provides the means to integrate particle trajectories through various
	magnet elements.  The focus of this simulation is spectrometer analysis.  As
	of the end of June 2014, the SNAKE model of the Hall-A Moller Polarimeter
	includes an analytic Helmholtz target field, analytic quadrupole magnets with
	fringe fields, an analytic dipole with fringe fields, and appropriate
	apertures up to and including the detector aperture.  The Moller scattering
	process is built into Snake.py, which provides the functionality to
	manipulate the polarization of beam/target in the scattering cross section,
	and generates two distinct trajectory files to track particle pairs through
	the spectrometer.  The geometry of this model was painstakingly extracted
	from the GEANT model, and it should be kept in mind that results from these
	two simulations should closely agree.


	Developments in the near future could include: modification of fringe
	fields, including a field map for the target field, and incorporating
	regions of field overlap.  

	The goal of this simulation is first to understand how the target field
	affects acceptance, and potentially develop methods to calibrate for this if
	the results are detrimental.  From there, this model could be used to
	investigate the spectrometer settings, specifically the quadrupole magnets.
	Another analysis on the horizon would be tracking error through the
	spectrometer based on uncertainties in element positions, magnetic fields,
	and initial beam properties.  It's also possible that this model could allow
	for a simulated measurement of beam polarization.

Generally ...
	s* - SNAKE source code
	x* - SNAKE executable
	d*.txt - data files
	*.dat - files called by executable
	m* - python modules to manipulate data
	p* - python plotting programs

To build SNAKE executable from source
	fort77 sSnake.f sSnakeLib.f -o xSnake

To run simulation
	- Executable xSnake calls (by default) the directive file (d.dat) and the
		trajectory file (t.dat), and for raytrace-style dipole magnets an
		accompanying dipole parameter file (dip.dat) as specified in the
		directive.  With these files in place, a simple ./xSnake should do it.
	- For ease of use, mSnake.py provides a method to control this entire
		process.  It generates the directive, trajectory, and dipole files,
		proceeds to run the executable, then reorganizes the output to obey the
		general rules described above.

Output files
	-	The executable .xSnake generates data files data.txt and pts.txt, which
		contain particle positions at endplanes and at integration steps
		respectively.  The format of data.txt is (trajectory number,
		endplane number, x, y, z), and that of pts.txt is (trajectory number,
		integration step number, x, y, z).  Typically, data.txt is used for
		quantitatively analyzing spectrometer acceptance, whereas pts.txt provides
		the points necessary to plot smooth trajectory paths in a simple
		visualizer.
	- The mSnake.py frontend renames these files from data.txt to
		dFse.txt/dRse.txt which describes the forward and reverse scattered
		electrons respectively, and from pts.txt to dFPts.txt/dRPts.txt again for
		the forward and reverse scattered electrons.  Note that the endplane data
		files are placed directly into snake/analysis/, while the integration point
		files are left in snake/ to be used with pSnake.py.  The mSnake.py program
		also generates the dAngs.txt file that records the cm scattering angles
		for each trajectory pair; this too is placed directly into snake/analysis/.


Analysis
	see readme.txt in snake/analysis
