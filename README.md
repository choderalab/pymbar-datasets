pymbar-datasets
===============

---------------
8proteins

Calculates the heat capacity as a function of temperature for a
simulation of 8 proteins.

Files:
analyzer.py     Calculates heat capacity as a function 
		of temperature from only a few states.
in_mbar		input files
out_mbar        output files

This data set originally caused a crash because the dynamic range of
the free energies was too large when the weights were stored in
non-logarithmic form.

---------------
alchemical-large-datasets

Files to be analyzed using alchemical-gromacs.py

Call is 

python alchemical-gromacs.py -d X -p dhdl -t 298

where X is the relevant directory.

Directory contains:
11states50ns/
    dhdl.*.xvg
    out11states.txt
38states/
    dhdl.*.xvg
    out38states.txt
