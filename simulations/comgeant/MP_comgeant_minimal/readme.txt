Welcome to MP_comgeant_minimal! - Ted (goofyboson@gmail.com) 6/14

First of all, the dependencies.  

	I can only confirm that this software runs on Ubuntu-12.04, though similar
	platforms should work as well.

	The crucial dependency is cernlib, which can be challenging to install.  I've
	written some basic installation instructions for Ubuntu; if you can't find
	them contact me.

	From there, I've developed a few python programs to make my life easier.  By
	no means are these necessary to run comgeant, but they're probably a decent
	place to start.  They require python, numpy, and matplotlib.

Background

	This simulation is based on comgeant specific to JLAB, originally adapted
	from GEANT3.  I have worked exclusively with executables handed down from 
	Sasha Glamazdin, and have no experience compiling these executables from
	source code.  From what I've been able to gather in some of the seemingly
	related source code, this model uses the completely unpolarized Moller
	scattering cross section in the CM frame.  It includes a field map for the
	super-conducting target magnet that was developed in physics storage at
	JLab; it has been theorized that there could be differences between this map
	and the field as it stands on the beamline due to surrounding ferrous
	materials.  It also seems to use uniform quadrupole and dipole magnets, by
	which I mean that it does not include fringe fields.  The accuracy of this
	simulation seems to be centered around models of the Moller Polarimeter
	detector, with particular focus on detector apertures, calorimetry, and
	detector response.

	This is what I've gathered from working with the code; I would be surprised
	if it was all correct.

Contents

	include/ contains the basic files necessary to run the comgeant executables

	run/ contains dedicated batch/ and interactive/ simulation modes
		Each simulation mode contains accompanying instructions.  I recommend
		starting with interactive/ as it's a simpler way to make sure the
		necessary dependencies are appropriately configured.

Good luck!
