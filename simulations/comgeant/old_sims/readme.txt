original comgeant files from Sasha Glamazdin ~ Spring 2013

differences could include polarimeter geometry and a few PAW analysis files

see sasha_instructions for just that

Subdirectories ...
	code: 32-bit comgeant executables and potential source code
					comg_batch.rhws3u.exe used for batch runs
					comg_inter.rhws3u.exe used for interactive runs
					unknown differences from 3b files
	code: 64-bit comgeant executables and source code
					never used!
	data: comgeant geometry, particle, material, and magnet files
	doc: Eugene's comgeant manual
	example: ?
	ntup_his: ntuples from batch jobs
	paw: comgeant run directory
		to run batch jobs in this directory, call comg_batch.rhws3u.exe,
			then, in PAW enter 'exe run_energy'
			omg_0/ffr contains adjustable parameters
			run_energy.kumac does analysis
		to run interactive jobs, call comg_inter.rhws3u.exe,
			enter 'exe plot_setup' followed by 'exe eve n=1000'
			or change n to however many events to produce
			again omg_0.ffr controls simulation parameters

			files not mentioned are probably called from those that were
				or maybe they're something else entirely
	run: some example interactive and batch runs with various analysis kumacs


