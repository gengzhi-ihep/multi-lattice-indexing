
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
   
   parser.add_option('-k','--node',dest='node',default="9",type='int',help='Parallel nodes number')
   parser.add_option('-n','--nlat',dest='nlat',default=5,type='int',help='Lattices number')
   parser.add_option('-e','--exe',dest='exe',default='/home/gengzhi/software/multi-indexing/bin/index-mpi',type='string',help='Path to executive index file')
   parser.add_option('-d','--mpi',dest='mpi',default='/usr/lib64/openmpi/bin/mpirun',type='string',help='Path to mpirun')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')

   (opts, args) = parser.parse_args()

   if(not os.path.exists(opts.exe)):
      print 'Error: Binary file index does not exist! Please use --exe=path-to-index-single.'
      sys.exit(1)

   if(not os.path.exists(opts.mpi)):
      print 'Error: mpirun does not exist! Please use --mpi=path-to-mpi.'
      sys.exit(1)

   if opts.res.endswith('/'):
      opts.res = opts.res.rstrip('/')

   if os.path.exists(opts.res) and os.path.isdir(opts.res):
      pass
   else:
      print 'Warning: result saving directory does not exist, a new directory will be created under current path.'
      os.mkdir(opts.res)

   curdir = os.getcwd()
   os.chdir(opts.res)

   command = 'cp '+ curdir+'/SPOT.TXT SPOT-TMP.TXT'
   os.system(command)

#   command = 'cp '+ curdir+'/param.config param.config'
#   os.system(command)

   command = opts.mpi
   args = " -np %d %s" % (opts.node, opts.exe)
   command += args

   with open('log','w') as f:

      for n in range(opts.nlat):

         print "================================================================"
         print "Try to find crystal: ", n+1
         f.write("================================================================\n")
         f.write("Try to find crystal: %d \n"%(n+1))

         p = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

         log = p.stdout.readlines()
         for line in log:
            print line.rstrip('\n')
            f.write(line)

         out, err = p.communicate()

#         cmd = "less best.txt | awk \'{printf \"%s %s %s %s %s %s %s %s \\n \",$1,$2,$3,$4,$5,$6,$7,n+1}\' > best.hkl "

         os.system('mv best.txt crystal-%d.hkl'%(n+1))

   f.close()

   sys.exit(0)

