
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
   parser.add_option('-n','--npat',dest='npat',default=1,type='int',help='Pattern number to index')
   parser.add_option('-r','--rot',dest='rot',default='1  0  0',type='string',help='Define rotation axes direction')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')

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
          line = "JOB= XYCORR INIT COLSPOT \n"
      elif "DATA_RANGE" in line:
          line = "DATA_RANGE=%d %d \n"%(1,opts.npat)
      elif "SPOT_RANGE" in line:
          line = "SPOT_RANGE=%d %d \n"%(opts.npat,opts.npat)
      elif "ROTATION_AXIS" in line:
          line = "ROTATION_AXIS= "+ opts.rot + '\n'

      fp.write(line)

   fp.close()

   command = "xds"
   p = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   out,err = p.communicate()

   command = "less SPOT.XDS | awk \'{printf \"%s %s \\n \",$1,$2}\' > SPOT.TXT"
   os.system(command)
   
