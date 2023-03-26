The <Simulate-error> directory is used for testing multi-indexing with artificial prior information errors. The detailed steps are listed as follows:

1. Simulate diffraction images containing multiple randomly distributed lattices.

   1.1 Enter the <fastbragg> directory, modify the target .f90 code to specify the number of lattices and error ranges for simulation. 
   1.2 Run ./configure to generate executable files. 
   1.3 Go back to <simulate_image> directory and modify the run.sh script to assign .mtz data file and executable file. 
   1.4 Modify 1generte-OM-cp.sh to assign unit cell information. 
   1.5 Modify 2sim-run-cp.sh to specify simulation parameters, such as distance, wavelength, pixel-size, et.al.
   1.6 Run run.sh to simulate images and perform peak-search job. This step will output new.h5,peaks.stream,om.txt,parameters.txt.
   1.7 Delete header information and the final two lines in peaks.stream and run ./awk.sh. This step will generate SPOT.TXT.

2. Indexing multi-lattice indexing.

   2.1 Copy SPOT.TXT to <index> directory.
   2.2 Modify param.config to assing necessary information. Refer to 1generate-OM-cp.sh and 2sim-run-cp.sh.
   2.3 Modify test.sh to assign number of lattices to search for.
   2.4 Run ./test.sh.

3. Compare indexing result with simulated ground-truth.

   Compare log and om.txt+parameters.txt to assess the accuracy of indexing.



The <Simulate-free> directory is used for testing multi-indexing with accurate prior information. The detailed steps are as follows:

   1.1 Enter the <fastbragg> directory, modify the target .f90 code to specify the number of lattices and error ranges for simulation. 
   1.2 Run ./configure to generate executable files. 
   1.3 Go back to <simulate_image> directory and modify the run.sh script to assign .mtz data file and executable file. 
   1.4 Modify 1generte-OM.sh to assign unit cell information. 
   1.5 Modify 2sim-run.sh to specify simulation parameters, such as distance, wavelength, pixel-size, et.al.
   1.6 Run run.sh to simulate images and perform peak-search job. This step will output new.h5,peaks.stream,om.txt.
   1.7 Delete header information and the final two lines in peaks.stream and run ./awk.sh. This step will generate SPOT.TXT.

2. Indexing multi-lattice indexing.

   2.1 Copy SPOT.TXT to <index> directory.
   2.2 Modify param.config to assing necessary information. Refer to 1generate-OM.sh and 2sim-run.sh.
   2.3 Modify test.sh to assign number of lattices to search for.
   2.4 Run ./test.sh.

3. Compare indexing result with simulated ground-truth.

   Compare log and om.txt to assess the accuracy of indexing.
