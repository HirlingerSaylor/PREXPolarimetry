This directory includes the absolute minimum number of files from Sasha's
simulations to run batch and interactive modes, and the incorporated analysis.

comg.exe - batch executable
comg_inter.exe - interactive executable

omg_0.ffr - parameter file stripped of interesting parameters
						used in batch and interactive modes to generate omg.ffr
omg.ffr - backup parameter file with interesting parameters filled
					just for reference

particles.ffr - particle definition file linked as fort.16
materials.ffr - material definition file linked as fort.17
geom_global.ffr - Moller Polarimeter geometry file linked as fort.21
									updated March 2013 (moller_geom_413.ffr)
geom_local.ffr - local geometry file, seems mostly empty but is still called
								 by executable so its still here, linked as fort.22
high_field.map - magnet map for superconducting helmholtz coils
								 linked as magnet_5.map

apnum.f, hitraxy.f, and siglg.f - fortran files used for PAW analysis in batch
