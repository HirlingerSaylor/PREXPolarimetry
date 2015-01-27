Batch Mode: For analysis of large monte-carlo data sets

omg.ffr - previous comgeant run parameters
					used to fill default parameters in batch.py
					if you lose it, there's a backup in include/

analysis.kumac - extensive PAW analysis, called by batch.py

batch.py - This does most of the work to generate large data sets

To run ...

- Run comgeant batch mode to generate data
	in MP_comgeant_minimal/run/batch/ enter "./batch.py" (without quotes)
	follow instructions to modify desired run parameters, then run the sim
	this does several things ...
		- uses comg.exe to generate comgeant ntuple, log, and his files in data/
			directory
		- uses analysis.kumac to generate plots of several interesting results with
			data pulled from the comgeant ntuple 
		- also uses analysis.kumac to dump the ntuple data from the interesting
			variables into individual data files in the data/ directory

- Analysis, which is pretty much up to the user.  The figures in data/figures/
	provide a decent sense of what's going on.  More detailed analysis is done
	using the dumped ntuple variables in data/python_analysis/.  See the
	instructions in data/python_analysis/ for more info.
