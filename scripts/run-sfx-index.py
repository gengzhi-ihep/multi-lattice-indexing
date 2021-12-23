
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
   
   parser.add_option('-k','--node',dest='node',default="9",type='int',help='Parallel nodes number')
   parser.add_option('-n','--nlat',dest='nlat',default=10,type='int',help='Maximun number of lattices')
   parser.add_option('-t','--nimg',dest='nimg',default=500,type='int',help='Number of images')
   parser.add_option('-e','--exe',dest='exe',default='/home/gengzhi/codes/index-search/software/bin/index',type='string',help='Path to executive index file')
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

   spots = np.loadtxt('SPOTS.TXT')

   curdir = os.getcwd()
   os.chdir(opts.res)

   command = 'cp '+ curdir+'/param.config param.config'
   os.system(command)

   command = opts.mpi
   args = " -np %d %s" % (opts.node, opts.exe)
   command += args

   with open('log','w') as f:

      for i in range(opts.nimg):

         sel = np.where(spots[:,0]==i+1)[0]
         spot = spots[sel]
         np.savetxt('SPOT-TMP.TXT',spot[:,1:3],fmt='%9.2f')

         print "======================================================"
         f.write("======================================================\n")
         n = 0
         while n < opts.nlat:

            p = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            log = p.stdout.readlines()
            if len(log) == 0 :
               break
            position_err = float(log[103].split()[-1])
            tmp1 = int(log[112].split()[0])
            tmp2 = int(log[112].split()[-1])
            match_rate = float(tmp1)/float(tmp2)

            if tmp1 < 10 and position_err > 1.0:
               break
            if tmp1 >= 10 and position_err > 1.5:
               break
            if tmp1 < 6:
               break

            print "Image: %d  lattice: %d "%(i+1,n+1)
            f.write("Image: %d  lattice: %d \n"%(i+1,n+1))

            f.write(log[107])
            f.write(log[108])
            f.write(log[109])
            f.write(log[110])
            f.write("Position Error: %6.3f  Match rate: %5.3f  %d / %d \n"%(position_err,match_rate,tmp1,tmp2))
            print "Position Error: %6.3f  Match rate: %5.3f  %d / %d"%(position_err,match_rate,tmp1,tmp2)

            out, err = p.communicate()

            if tmp2-tmp1 < 5:
               break

            n = n + 1

#            cmd = "less best.txt | awk \'{printf \"%s %s %s %s %s %s %s %s \\n \",$1,$2,$3,$4,$5,$6,$7,n+1}\' > best.hkl "
#            os.system('mv best.txt crystal-%d.hkl'%(n+1))


   f.close()

   sys.exit(0)
