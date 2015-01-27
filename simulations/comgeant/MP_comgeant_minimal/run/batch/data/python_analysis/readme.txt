Generally ...
	*.dat files are single variable files dumped from comgeant ntuple
	d* files are 2nd tier data files, build from comgeant *.dat files
	m* files describe python modules, which do some manipulation of data files
	p* files plot some selection of the data

After running batch.py from MP_comgeant_minimal/run/batch/, mMerge.py first
generates dInitialAngles.txt, dCoincidentData.txt, and dEnergyByChannels.txt
which are used by several plotting files.  Then, to examine data based on
channel selections, mChannels.py segregates events based on channels with the
highest energy; it produces dDataChannel*.  If you get errors about these files
not existing, it's probably because you didn't run mMerge.py or mChannels.py.

The various plotting files are fairly self explanatory, but if you have any
questions just run them and see what happens.
