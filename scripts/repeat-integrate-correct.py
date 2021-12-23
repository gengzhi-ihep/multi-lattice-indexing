
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

if __name__ == "__main__":
   
   usage="python %prog [options]"
   parser = optparse.OptionParser(usage)

   parser.add_option('-p','--path',dest='path',default="/home/gengzhi/dataset/lyso-3crystals",type='string',help='Path to diffraction images')
   parser.add_option('-f','--head',dest='head',default="lyso-3crystals01.????",type='string',help='Header name of images')
   parser.add_option('-t','--rang',dest='rang',default='1  30',type='string',help='Start and End number for integration')
   parser.add_option('-s','--spgn',dest='spgn',default=0,type='int',help='Space group number')
   parser.add_option('-r','--irot',dest='irot',default='1  0  0',type='string',help='Define rotation axes direction')
   parser.add_option('-u','--unit',dest='unit',default='70.0  80.0  90.0  90.0  90.0  90.0',type='string',help='Define Unit Cell parameters')
   parser.add_option('-o','--orgn',dest='orgn',default='1024.0 1024.0',type='string',help='Define beam center x and y')
   parser.add_option('-z','--reso',dest='reso',default='50  0',type='string',help='Resolution range')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')
   parser.add_option('-a','--om1',dest='iom1',default='1  0  0',type='string',help='Define the first line of inverse orientation matrix')
   parser.add_option('-b','--om2',dest='iom2',default='0  1  0',type='string',help='Define the second line of inverse orientation matrix')
   parser.add_option('-c','--om3',dest='iom3',default='0  0  1',type='string',help='Define the third line of inverse orientation matrix')
   parser.add_option('-n','--ntry',dest='ntry',default=1,type='int',help='Iteration number of integration')

   (opts, args) = parser.parse_args()

   if os.path.exists(opts.res) and os.path.isdir(opts.res):
      pass
   else:
      print 'Warning: result saving directory does not exist, a new directory will be created under current path.'
      os.mkdir(opts.res)

   os.chdir(opts.res)

   command = "generate_XDS.INP" 
   image = opts.path+ "/" + opts.head

   p = subprocess.Popen([command,image],stdout=subprocess.PIPE)
   out,err = p.communicate()

   lines = open('XDS.INP').readlines()
   fp = open('XDS.INP','w')
   for line in lines:
      if "JOB" in line:
          line = "JOB= XYCORR INIT DEFPIX INTEGRATE CORRECT \n"
      elif "ORGX" in line:
          [orgx,orgy] = opts.orgn.split()
          line = "ORGX= "+orgx+" ORGY= "+orgy+'\n'
      elif "DATA_RANGE" in line:
          line = "DATA_RANGE= " + opts.rang + '\n'
      elif "SPOT_RANGE" in line:
          line = "SPOT_RANGE= 1  1 \n"
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

   fp = open('XPARM.XDS','w')
   fp.write(" XPARM.XDS    VERSION Jan 31, 2020  BUILT=20200417 \n")
   fp.write('    1 '+start_angle+'  '+osc_range+'  '+opts.irot+'\n')
   fp.write('      '+wavelen+'  '+'  0.000    0.000    1.000 \n')
   fp.write('   '+str(opts.spgn)+'  '+opts.unit+'\n')
   fp.write('   '+opts.iom1+'\n')
   fp.write('   '+opts.iom2+'\n')
   fp.write('   '+opts.iom3+'\n')
   fp.write('     1  '+NX+'  '+NY+'  '+QX+'  '+QY+'\n')
   fp.write('   '+opts.orgn+'  '+distance+'\n')
   fp.write('     1.000000     0.000000     0.000000 \n')
   fp.write('     0.000000     1.000000     0.000000 \n')
   fp.write('     0.000000     0.000000     1.000000 \n')
   fp.write('       1       1    '+NX+'    '+' 1   '+NY+'\n')
   fp.write('   0.00   0.00  0.00  1.00000  0.00000  0.00000  0.00000  1.00000  0.00000 \n')
   fp.close()

   command = "xds"
   for n in range(opts.ntry):
      p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
      out,err = p.communicate()
      os.system('mv GXPARM.XDS XPARM.XDS')
