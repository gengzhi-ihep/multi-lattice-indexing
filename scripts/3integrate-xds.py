
"""
 (c) IHEP 2021. All rights reserved.
 Author: Geng Zhi
 Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS).
 If you have any problem with the program, please contact author with the following
 email: gengz@ihep.ac.cn
"""

import subprocess
import os
import sys
import optparse
import numpy as np

if __name__ == "__main__":
   
   usage="python %prog [options]"
   parser = optparse.OptionParser(usage)

   parser.add_option('-n','--npat',dest='npat',default=1,type='int',help='Pattern number for index')
   parser.add_option('-t','--rang',dest='rang',default='1  30',type='string',help='Start and End number for integration')
   parser.add_option('-s','--spgn',dest='spgn',default=0,type='int',help='Space group number')
   parser.add_option('-u','--unit',dest='unit',default='70.0  80.0  90.0  90.0  90.0  90.0',type='string',help='Define Unit Cell parameters')
   parser.add_option('-r','--irot',dest='irot',default='1  0  0',type='string',help='Define rotation axes direction')
   parser.add_option('-o','--orgn',dest='orgn',default='1024.0 1024.0',type='string',help='Define beam center x and y')
   parser.add_option('-z','--reso',dest='reso',default='50  0',type='string',help='Resolution range')
   parser.add_option('-e','--res',dest='res',default='./results',type='string',help='Path to saving result')
   parser.add_option('-m','--om',dest='om',default='1  0  0; 0  1  0; 0  0  1',type='string',help='Input orientation matrix')

   (opts, args) = parser.parse_args()

   if os.path.exists(opts.res) and os.path.isdir(opts.res):
      pass
   else:
      print 'Warning: result saving directory does not exist, a new directory will be created under current path.'
      os.mkdir(opts.res)

   os.chdir(opts.res)

   lines = open('XDS.INP').readlines()
   fp = open('XDS.INP','w')
   for line in lines:
      if "JOB" in line:
          line = "JOB= DEFPIX INTEGRATE \n"
      elif "ORGX" in line:
          [orgx,orgy] = opts.orgn.split()
          line = "ORGX= "+orgx+" ORGY= "+orgy+'\n'
      elif "DATA_RANGE" in line:
          line = "DATA_RANGE= " + opts.rang + '\n'
#      elif "SPOT_RANGE" in line:
#          line = "SPOT_RANGE= 1  1 \n"
      elif "SPACE_GROUP_NUMBER" in line:
          line = "SPACE_GROUP_NUMBER= %d \n" % opts.spgn
      elif "UNIT_CELL_CONSTANTS" in line:
          line = "UNIT_CELL_CONSTANTS= "+opts.unit+'\n'
      elif "ROTATION_AXIS" in line:
          line = "ROTATION_AXIS= "+ opts.irot + '\n'
      elif "INCLUDE_RESOLUTION_RANGE=" in line:
          line = "INCLUDE_RESOLUTION_RANGE= " + opts.reso + '\n'
      elif "STARTING_ANGLE" in line:
          temp = line.split()
          start_angle = temp[1]
      elif "OSCILLATION_RANGE=" in line:
          temp = line.split()
          osc_range = temp[1]
      elif "X-RAY_WAVELENGTH=" in line:
          temp = line.split()
          wavelen = temp[1]
      elif "DETECTOR_DISTANCE=" in line:
          temp = line.split()
          distance = temp[1]
      elif "NX=" in line:
          temp = line.rstrip('\n').split()
          NX = temp[1]
          NY = temp[3]
          QX = temp[5]
          QY = temp[7]

      fp.write(line)

   fp.close()

   om = opts.om.split(';')
   
   m = np.array([float(b) for a in om for b in a.split()]).reshape((3,3))
   iom = np.linalg.inv(m)
   iom1 = "%12.6f   %12.6f   %12.6f"% (iom[0,0],iom[0,1],iom[0,2])
   iom2 = "%12.6f   %12.6f   %12.6f"% (iom[1,0],iom[1,1],iom[1,2])
   iom3 = "%12.6f   %12.6f   %12.6f"% (iom[2,0],iom[2,1],iom[2,2])

   fp = open('XPARM.XDS','w')
   fp.write(" XPARM.XDS    VERSION Feb 5, 2021  BUILT=20210323 \n")
   fp.write('   '+str(opts.npat)+'  '+start_angle+'  '+osc_range+'  '+opts.irot+'\n')
   fp.write('      '+wavelen+'  '+'  0.000    0.000    1.000 \n')
   fp.write('   '+str(opts.spgn)+'  '+opts.unit+'\n')
   fp.write('   '+iom1+'\n')
   fp.write('   '+iom2+'\n')
   fp.write('   '+iom3+'\n')
   fp.write('     1  '+NX+'  '+NY+'  '+QX+'  '+QY+'\n')
   fp.write('   '+opts.orgn+'  '+distance+'\n')
   fp.write('     1.000000     0.000000     0.000000 \n')
   fp.write('     0.000000     1.000000     0.000000 \n')
   fp.write('     0.000000     0.000000     1.000000 \n')
   fp.write('       1       1    '+NX+'    '+' 1   '+NY+'\n')
   fp.write('   0.00   0.00  0.00  1.00000  0.00000  0.00000  0.00000  1.00000  0.00000 \n')
   fp.close()

   command = "xds"
   p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   out,err = p.communicate()
