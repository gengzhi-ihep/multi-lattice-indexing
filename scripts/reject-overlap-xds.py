
"""
 (c) IHEP 2021. All rights reserved.
 Author: Geng Zhi
 Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS).
 If you have any problem with the program, please contact author with the following
 email: gengz@ihep.ac.cn
"""

import numpy as np
import optparse
import math
import os

if __name__ == "__main__":

   usage="python %prog [options]"
   parser = optparse.OptionParser(usage)

   parser.add_option('-f','--inpf',dest='inpf',default="INTEGRATE1.HKL INTEGRATE2.HKL",type='string',help='A list of input unmerged files to remove overlap')
   parser.add_option('-t','--cut',dest='cut',default=3.0,type='float',help='Threshold for spot overlapping')
   parser.add_option('-m','--res',dest='res',default='./results',type='string',help='Path to saving result')
   parser.add_option('-n','--npat',dest='npat',default=360,type='int',help='Number of collected diffraction patterns')
   parser.add_option('-r','--nref',dest='nref',default=1000,type='int',help='Maximum number of reflections for each pattern')
   


   (opts, args) = parser.parse_args()

   if os.path.exists(opts.res) and os.path.isdir(opts.res):
      pass
   else:
      print 'Warning: result saving directory does not exist, a new directory will be created under current path.'
      os.mkdir(opts.res)

   os.chdir(opts.res)

   files = opts.inpf.split()

   num = len(files)
   
   N_images=opts.npat
   max_refls=opts.nref

   nrefls = np.zeros((N_images,num))
   
   xdata = np.zeros((max_refls,N_images,num))
   ydata = np.zeros((max_refls,N_images,num))

   for i in range(num):
      lines = open(files[i]).readlines()
      for line in lines[32:-1]:
         tmp = line.split()
         seq = int(float(tmp[14])+0.5) - 1
         if seq<0: seq = 0
         if seq>N_images-1: seq = N_images-1
         xdata[int(nrefls[seq,i]),seq,i] = float(tmp[12])
         ydata[int(nrefls[seq,i]),seq,i] = float(tmp[13])
         nrefls[seq,i] += 1


   reject = np.zeros((max_refls,N_images,num))

   for i in range(num-1):
      for j in range(i+1, num):

         for n1 in range(N_images):

            if int(nrefls[n1,i]) == 0: continue
            if int(nrefls[n1,j]) == 0: continue

            for n2 in range(int(nrefls[n1,i])):
  
               for n3 in range(int(nrefls[n1,j])):

                  dis = math.sqrt((xdata[n2,n1,i]-xdata[n3,n1,j])**2+(ydata[n2,n1,i]-ydata[n3,n1,j])**2)
                  if dis <= opts.cut:
                     reject[n2,n1,i] = 1
                     reject[n3,n1,j] = 1

   nrefls = np.zeros((N_images,num))
   for i in range(num):
      lines = open(files[i]).readlines()
      fp = open(files[i].split('.')[0]+'_REM.HKL','w')
      for line in lines[0:32]:
         fp.write(line)
      for line in lines[32:-1]:
         tmp = line.split()
         seq = int(float(tmp[14])+0.5) - 1
         if seq<0: seq = 0
         if seq>N_images-1: seq = N_images-1
         if reject[int(nrefls[seq,i]),seq,i]< 0.5:
            fp.write(line)
         nrefls[seq,i] += 1
      fp.write(lines[-1])
      fp.close()
