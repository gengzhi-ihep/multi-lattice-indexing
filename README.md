MCDPS: A multi-crystal indexing data processing suite in macromolecucar crystallography
======================================================================================
Author: Zhi Geng | Email address: gengz@ihep.ac.cn

Description:
----------------------------------------------------------------------------------
`MCDPS` is a software used for indexing diffraction images containing multiple lattices in macromolecular crystallography. It has several important features:

* It is based on the principle of whole-pattern matching provided cell dimensions and space-group symmetry are known.
* It features a local correction for prior information as well as an iterative refinement of experimental parameters to enhance its robustness to errors.
* Indexing result of the new method can be fed back to softwares _XDS_ and _Crystfel_ to facilitate further data reduction and structural determination.

Please cite the following publication if you use ***MCDPS*** in your work: 😊
>Zhou,Q.,Gao,Z.Q.,Dong,Z.,Jiang,Y.M.,She,Z.,Geng,Z.&Dong,Y.H.(2021). A reference-based multi-lattice indexing method integrating prior information correction and iterative refinement in protein crystallography. Acta Cryst. A77. [https://doi.org/10.1107/S2053273321003521](https://doi.org/10.1107/S2053273321003521)

Prerequsites:
-----------------------------------------------------------------------------------
Before compiling source codes, the following libraries or softwares are required:

* python
* gfortran
* openmpi
* [gsl-1.x](http://www.gnu.org/software/gsl/)
* [fgsl-1.0.0](https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library)
* [XDS](https://xds.mr.mpg.de/) (for oscillation datasets)
* [Crystfel](https://www.desy.de/~twhite/crystfel/) (for single-shot datasets)
* [NanoBragg](https://bl831.als.lbl.gov/~jamesh/fastBragg/) (for simulation)

Installation:
------------------------------------------------------------------------------------

1. Install library: `openmpi`, `gsl` and `fgsl`.

* `fgsl` can be downloaded from: [https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library](https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library)
* `gsl` and `openmpi` can be installed through `yum` or downloaded from source codes.  
*  Note: version of `fgsl` must coincide with version of `gsl` (fgsl-1.0.0-->gsl-1.x)

2. Modify `configure` to specify library path to  `openmpi` and `fgsl`.

3. run `./configure`.  

4. The following 2 binary files will be created in directory `bin`:  
`index`,`prepare-sfx`

5. Finish. Now you can enjoy your test! 😄

Usage:
--------------------------------------------------------------------------------

1. Indexing oscillation (or rotation) based diffraction datasets:

* `Preparation`: Prepare a configuration file `param.config` specifying the following parameters:

```
###information about space group and unit cell
Space_group = 19       # space group number
Crystal_system = Tetragonal  #Triclinic,Monoclinic,Orthorhombic,Trigonal,Tetragonal,Hexagonal,Cubic
Unit_cell_a = 77.0    # unit cell parameter a (unit: angstrom)
Unit_cell_b = 77.0    # unit cell parameter b
Unit_cell_c = 37.2   # unit cell parameter c
Unit_cell_alpha = 90.00  # unit cell parameter alpha (unit:degree)
Unit_cell_beta = 90.00   # unit cell parameter beta
Unit_cell_gama = 90.00   # unit cell parameter gama

###information about diffraction geometry
Wavelength = 1.0      # X-ray wavelength (unit: angstrom)
Distance = 140.0      # sample-to-detector distance (unit: millimeter)
Pixel_size = 0.079    # single pixel size (unit: millimeter)
XCenter = 1034.84     # beam center x (unit: pixel)
YCenter = 1028.31     # beam center y (unit: pixel)

###information about user implementation parameters
Enable_rescut = 1     # 1: perform resolution cutoff for indexing  0: not perform
Rescut = 5.0          # setup resolution cutoff
Enable_refcut = 0     # 1: perform reflection number cutoff for indexing  0: not perform
Refcut = 90           # setup reflection number cutoff
Position_error = 3.5  # positional tolerance to accept reflections as matched (unit: pixel)

###prior information correction parameters
Enable_gridscan = 0   # 1: perform prior information correction   0:not perform
Enable_uc_scan = 0    # 1: perform unit cell correction   0:not perform
CellScan = 1.0        # unit cell scan range (unit: angstrom)
CellScanstep = 0.2    # unit cell scan step (unit: angstrom)
XCScan = 6            # beam center x scan range (unit: pixel)
YCScan = 6            # beam center y scan range (unit: pixel)
XYCScanstep = 2       # beam center scan step (unit: pixel)
Local_angle_range = 5.0 # local orientation search range (unit: degree) 
Local_angle_step = 1.0  # local orientation search step (unit: degree)
```

* `Peak search`: generate `SPOT.TXT` with script `1peak-search-xds.py`, for example:
```
python 1peak-search-xds.py -p "/home/gengzhi/dataset/lyso-3crystals" -f "lyso-3crystals01.????" -n 1 -r "1  0  0"  -m "./results"
```
Options you may set are:
```
Usage: python 1peak-search-xds.py [options]

Options:
  -h, --help              show this help message and exit.
  -p PATH, --path=PATH    Path to where diffraction images saved.  Compulsory.
  -f HEAD, --head=HEAD    Header name of diffraction images.  Compulsory.
  -n NPAT, --npat=NPAT    Individual image number for MCDPS indexing.  Optional, default is 1.
  -r ROT, --rot=ROT       Define rotation axis direction, depending on station setup. Optional, default is '1  0  0'.
  -m RES, --res=RES       Path to saving indexing result. Optional, default is './results'.
```

* `Indexing`: Using both `param.config` and `SPOT.TXT`, ***MCDPS*** can be performed using script `2run-index.py`, for example:
```
python 2run-index.py -k 9 -n 5 -e "/home/gengzhi/codes/index-search/software/bin/index" -d "/usr/lib64/openmpi/bin/mpirun" -m "./results"
```
Options you may set are:
```
Usage: python 2run-index.py [options]

Options:
  -h, --help              show this help message and exit.
  -k NODE, --node=NODe    Number of CPU cores for parallel computing. Optional, default is 9.
  -n NLAT, --nlat=NLAT    Number of target crystal lattices to index. Optional, default is 5.
  -e EXE, --exe=EXE       Path to binary executive file-->index. Compulsory.
  -d MPI, --mpi=MPI       Path to parallel compiler-->mpirun. Compulsory.
  -m RES, --res=RES       Path to saving indexing result. Optional, default is './results'.
```

* `Integration`: Given orientation matrix, integrate all diffraction images with script `3integrate-xds.py`, for example:
```
python 3integrate-xds.py -s 96 -u "77.0 77.0 37.24 90.0 90.0 90.0" -t "1  30" -r "1  0  0"  -o "1034.84 1028.31" -z "50  0" -a "-1.3053 -33.4794 69.3425" -b "47.555 54.1956 27.0615" -c "-29.2227 20.8823 9.5321" -m "./results"
```
Options you may set are:
```
Usage: python 3integrate-xds.py [options]

Options:
  -h, --help              show this help message and exit.
  -t RANG, --rang=RANG    Start and End image number used for integration. Optional, default is '1  30'.
  -s SPGN, --spgn=SPGN    Space group number. Compulsory.
  -u UNIT, --unit=UNIT    Unit Cell parameters. Compulsory.
  -r IROT, --irot=IROT    Rotation axis direction dependent on station setup. Optional, default is '1  0  0'.
  -o ORGN, --orgn=ORGN    Beam center X and Y in unit of pixels. Compulsory.
  -z RESO, --reso=RESO    Resolution range. Optional, default is '50  0'.
  -m RES, --res=RES       Path to saving integration result. Optional, default is './results'.
  -a IOM1, --om1=IOM1     The first line of inverse orientation matrix. Compulsory.
  -b IOM2, --om2=IOM2     The second line of inverse orientation matrix. Compulsory.
  -c IOM3, --om3=IOM3     The third line of inverse orientation matrix. Compulsory.
```

* `Scale and Merge`: Scale and merge integrated intensities for each crystal with script `4correct-xds.py`, for example:
```
python 4correct-xds.py -s 96 -u "77.0 77.0 37.24 90.0 90.0 90.0" -z "50  0" -m "./results"
```
Options you may set are:
```
Usage: python 4correct-xds.py [options]

Options:
  -h, --help              show this help message and exit.
  -s SPGN, --spgn=SPGN    Space group number. Compulsory.
  -u UNIT, --unit=UNIT    Unit Cell parameters. Compulsory.
  -z RESO, --reso=RESO    Resolution range. Optional, default is '50  0'.
  -m RES, --res=RES       Path to saving merged result. Optional, default is './results'.
```

* `Multiple datasets Merge`: Scale and merge multiple integrated intensities to generate more complete dataset with script `5multi-merge.py`, for example:
```
python 5multi-merge.py -f "XDS_ASCII_1.HKL XDS_ASCII_2.HKL" -u "77.0 77.0 37.24 90.0 90.0 90.0" -s 96 -t "TRUE" -m "./results"
```
Options you may set are:
```
Usage: python 5multi-merge.py [options]

Options:
  -h, --help              show this help message and exit.
  -f INPF, --inpf=INPF    A list of all input unmerged XDS_ASCII.HKL filenames. Compulsory.
  -u UNIT, --unit=UNIT    Unit Cell parameters. Compulsory.
  -s SPGN, --spgn=SPGN    Space group number. Compulsory.
  -t FRED, --fred=FRED    Fridel Law: TRUE or False. Optional, default is TRUE.
  -m RES, --res=RES       Path to saving merged result. Optional, default is './results'.
```

* `Remove overlaps`: Remove overlapping reflections between multiple integrated files with script `reject-overlap-xds.py`, for example:
```
python reject-overlap-xds.py -f "INTEGRATE1.HKL  INTEGRATE2.HKL INTEGRATE3.HKL" -t 3.0 -m "./results"
```
Options you may set are:

```
Usage: python reject-overlap-xds.py [options]

Options:
  -h, --help              show this help message and exit.
  -f INPF, --inpf=INPF    A list of all input unmerged XDS_ASCII.HKL files to remove overlaps. Compulsory.
  -t CUT, --cut=CUT       Criterion of overlapping between two reflections in unit of pixels. Optional, default is 3.0.
  -m RES, --res=RES       Path to saving result. Optional, default is './results'.
```

* `show clusters`: Show hierarchy clustering results of match-rate after the first round of whole space search,for example:
```
python3 show_cluster.py -f "class.txt" -c 0.2 -t 30
```
Options you may set are:

```
Usage: python3 show_cluster.py [options]

Options:
  -h, --help              show this help message and exit.
  -f FILE, --file=FILE    a file listing all three orientation angles with corresponding match-rate (in total 4 columns). Compulsory.
  -c MCUT, --cut=MCUT     cutoff value to select high match-rate for display. Compulsory, default is 0.2.
  -t TALL, --tall=TALL    Threshold to cluster neigher points in hierarchy clustering.
```

2. Indexing single-shot (or serial) diffraction images:

* `Preparation`: Prepare a configuration file `param.config` specifying the following parameters:

```
###information about space group and unit cell
Space_group = 96       # space group number
Crystal_system = Tetragonal  #Triclinic,Monoclinic,Orthorhombic,Trigonal,Tetragonal,Hexagonal,Cubic
Unit_cell_a = 79.17    # unit cell parameter a (unit: angstrom)
Unit_cell_b = 79.17    # unit cell parameter b
Unit_cell_c = 38.01   # unit cell parameter c
Unit_cell_alpha = 90.00  # unit cell parameter alpha (unit:degree)
Unit_cell_beta = 90.00   # unit cell parameter beta
Unit_cell_gama = 90.00   # unit cell parameter gama

###information about diffraction geometry
Wavelength = 1.3275      # X-ray wavelength (unit: angstrom)
Distance = 93.0      # sample-to-detector distance (unit: millimeter)
Pixel_size = 0.110    # single pixel size (unit: millimeter)
XCenter = 0.0     # beam center x (unit: pixel)
YCenter = 0.0     # beam center y (unit: pixel)

###information about user implementation parameters
Enable_rescut = 1     # 1: perform resolution cutoff for indexing  0: not perform
Rescut = 3.0          # setup resolution cutoff
Enable_refcut = 1     # 1: perform reflection number cutoff for indexing  0: not perform
Refcut = 90           # setup reflection number cutoff
Position_error = 3.5  # positional tolerance to accept reflections as matched (unit: pixel)

###prior information correction parameters
Enable_gridscan = 0   # 1: perform prior information correction   0:not perform
Enable_uc_scan = 0    # 1: perform unit cell correction   0:not perform
CellScan = 1.0        # unit cell scan range (unit: angstrom)
CellScanstep = 0.2    # unit cell scan step (unit: angstrom)
XCScan = 6            # beam center x scan range (unit: pixel)
YCScan = 6            # beam center y scan range (unit: pixel)
XYCScanstep = 2       # beam center scan step (unit: pixel)
Local_angle_range = 5.0 # local orientation search range (unit: degree) 
Local_angle_step = 1.0  # local orientation search step (unit: degree)
```
* Generate `SPOTS.TXT` from two files: `.stream` from Crystfel-0.7.0 and `.geom` describing detector geometry:
```
prepare-sfc
```

* `Indexing`: Using both `param.config` and `SPOT.TXT`, ***MCDPS*** can be performed using script `run-sfx-index.py`, for example:
```
python run-sfx-index.py -k 9 -n 5 -t 500 -e "/home/gengzhi/codes/index-search/software/bin/index" -d "/usr/lib64/openmpi/bin/mpirun" -m "./results"
```
Options you may set are:

```
Usage: python reject-overlap-xds.py [options]

Options:
  -h, --help              show this help message and exit.
  -k NODE, --node=NODe    Number of CPU cores for parallel computing. Optional, default is 9.
  -n NLAT, --nlat=NLAT    Number of target crystal lattices to index. Optional, default is 5.
  -t NIMG, --nimg=NIMG    Number of images used for indexing. Compulsory, default 500.
  -e EXE, --exe=EXE       Path to binary executive file-->index. Compulsory.
  -d MPI, --mpi=MPI       Path to parallel compiler-->mpirun. Compulsory.
  -m RES, --res=RES       Path to saving indexing result. Optional, default is './results'.
```