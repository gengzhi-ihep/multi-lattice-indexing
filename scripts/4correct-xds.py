
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

   parser.add_option('-s','--spgn',dest='spgn',default=0,type='int',help='Space group number')
   parser.add_option('-u','--unit',dest='unit',default='70.0  80.0  90.0  90.0  90.0  90.0',type='string',help='Define Unit Cell parameters')
   parser.add_option('-z','--reso',dest='reso',default='50  0',type='string',help='Resolution range')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')

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
          line = "JOB= CORRECT \n"
      elif "SPACE_GROUP_NUMBER" in line:
          line = "SPACE_GROUP_NUMBER= %d \n" % opts.spgn
      elif "UNIT_CELL_CONSTANTS" in line:
          line = "UNIT_CELL_CONSTANTS= "+opts.unit+'\n'
      elif "INCLUDE_RESOLUTION_RANGE=" in line:
          line = "INCLUDE_RESOLUTION_RANGE= " + opts.reso + '\n'

      fp.write(line)

   fp.close()

   command = "xds"
   p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   out,err = p.communicate()
