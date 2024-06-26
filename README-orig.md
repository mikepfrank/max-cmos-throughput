# Maximum Switching Throughput Density Estimator

The file `IRDS2022_HPHD.m` is a MATLAB script for estimating the
maximum switching throughput density of CMOS logic as a function of
power density based on IRDS roadmap data for various technology
nodes. It covers four basic design scenarios: (1) conventional
switching at standard voltages, (2) conventional switching at
optimized voltages, (3) fully adiabatic switching at standard
voltages, (4) fully adiabatic switching at optimized voltages.

This script was initially developed in July-August, 2022 by Alexander 
J. Edwards (AJE) with support from a DOE SCGSR Fellowship, with additional
contributions by Michael P. Frank (MPF) with support from the DOE/NNSA
ASC Nonconventional Computation project. Additional revisions were
made by MPF in Sep. 2022 and Sep. 2023.

The variable `is_HP` near to start of the script can be set to 1 to
utilize the data for the IRDS "high-performance" (HP) design scenario, 
or to 0 to utilize the data for the IRDS "high-density" (HD) design scenario.
(By default we utilize the HD data, because it yields greater throughput
in the power-constrained setting of interest.)

The preliminary results derived from this script were presented at the 
Texas Symposium on Computing with Emerging Technologies (CoMET) in
Aug. 2023, in an invited talk presented by MPF titled "Limits of CMOS
and Prospects for Adiabatic/Reversible CMOS."

The following presentation materials are publicly available:
* [Slide deck in PDF](https://www.sandia.gov/app/uploads/sites/210/2023/11/Comet23-slides_SAND.pdf)
* [Video of presentation](https://www.youtube.com/watch?v=vALCJJs9Dtw)