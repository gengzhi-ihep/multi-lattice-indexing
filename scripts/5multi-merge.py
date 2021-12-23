
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
   parser.add_option('-f','--inpf',dest='inpf',default="XDS_ASCII_1.HKL XDS_ASCII_2.HKL",type='string',help='A list of input unmerged XDS_ASCII.HKL filenames')
   parser.add_option('-u','--unit',dest='unit',default='70.0  80.0  90.0  90.0  90.0  90.0',type='string',help='Unit Cell parameters')
   parser.add_option('-s','--spgn',dest='spgn',default=0,type='int',help='Space group number')
   parser.add_option('-t','--fred',dest='fred',default='TRUE',type='string',help='Friedel law:TRUE OR FALSE')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')

   (opts, args) = parser.parse_args()

   if os.path.exists(opts.res) and os.path.isdir(opts.res):
      pass
   else:
      print 'Warning: result saving directory does not exist, a new directory will be created under current path.'
      os.mkdir(opts.res)

   files = opts.inpf.split()
   for i in range(len(files)):
      command = "cp "+ files[i]+'  '+opts.res
#      os.system(command)

   os.chdir(opts.res)

   fp = open('XSCALE.INP','w')
   fp.write("SPACE_GROUP_NUMBER= %d \n" % opts.spgn)
   fp.write("UNIT_CELL_CONSTANTS= "+opts.unit+'\n')
   fp.write("OUTPUT_FILE= temp.ahkl \n")
   for i in range(len(files)):
      fp.write("INPUT_FILE= "+files[i]+'\n')
   fp.close()

   command = "xscale"
   p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   out,err = p.communicate()

   fp = open('XDSCONV.INP','w')
   fp.write("SPACE_GROUP_NUMBER= %d \n" % opts.spgn)
   fp.write("UNIT_CELL_CONSTANTS= "+opts.unit+'\n')
#   fp.write("INPUT_FILE= temp.ahkl \n")
   fp.write("INPUT_FILE= XDS_ASCII.HKL \n")
   fp.write("OUTPUT_FILE= temp.hkl CCP4_F \n")
   fp.write("FRIEDEL'S_LAW="+opts.fred+'\n')
   fp.write("MERGE="+opts.fred+'\n')
   fp.close()

   command = "xdsconv"
   p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   out,err = p.communicate()

   command = "f2mtz HKLOUT temp.mtz"
   p = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   params = open("F2MTZ.INP","r").readlines()
   args = ""
   for arg in params:
       args += str(arg)
   p.stdin.write(args)
   out,err = p.communicate()

   command = "cad HKLIN1 temp.mtz HKLOUT output_file_name.mtz"
   p = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   args = "LABIN FILE 1 ALL"
#   args += "DWAVELENGTH FILE 1 1   0.97918"
   p.stdin.write(args)
   out,err = p.communicate()
