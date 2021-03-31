* Multi-lattice-indexing

* Author: Dr.Geng Zhi (Email address: gengz@ihep.ac.cn)

* Multi-lattice indexing is used for indexing diffraction images containing multiple lattices in macromolecular crystallography.

* Please cite the following publication if you use those codes in your work.

* Zhou,Q.,Gao,Z.Q.,Dong,Z.,Jiang,Y.M.,She,Z.,Geng,Z. and Dong,Y.H. A reference-based multi-lattice indexing method integrating prior information

* correction and iterative refinement in protein crystallography.

################################
* Prerequisite:

* gfortran

* openmpi

* gsl-1.15

* fgsl-1.0.0

################################


################################
* Installation:

* Modify the file 'configure' to assign your path to binary mpif90, library of openmpi, gsl and fgsl.

* Then run ./configure. After that, two binary files will be created in ./bin.

################################

################################
* How to use:

* A list of portable python scripts are used to run all steps, including peak search with XDS, run multi-lattice indexing, 

* data integration and merging with XDS, reject overlapped reflections between multiple files and indexing of SFX data.

* The python scripts are placed in ./scripts.

* An example of using these scripts is provided in a file named 'example_command.txt'

* When indexing with above scripts, a spot positions list ('SPOT.TXT') containing peak coordinates in unit of pixels must be provided. 

* Another required file is a configuration file ('params.config').

################################
