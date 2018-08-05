# Summer Project 2017
Work done in Python for the summer project in the School of Physics in the summer of 2017.

`findSN.py` allows the user to change the time differences between detectors (in source code) and find the corresponding Supernova.  This was the first program written, when I was trying to get to grip with the geometry and mathematics of the problem.

`sampleplotSN.py` is a more complex simulation expanding on previous work by adding uncertainties.  Or get it touch if you have questions about it. 

`sample_dt.py` samples neutrino start times from real detector data.

To run one of the programs as usual type:

`python yourfilename.py <optional input files>`

in the terminal.  `sampleplotSN.py` and `sample_dt.py` will require input files appended to the end of the command as well.  The input file `t1.out` contain the differential rate data and `t1thresh.out` contains the digitised energy threshold data.  Similarly for t2, t3 and t4.  The different numbers indicate different models being used.

The `astropy` (http://www.astropy.org/) library, `ROOT` (https://root.cern.ch/) and `scipy` (https://www.scipy.org/) will be required to run all Python programs in this repository. 
