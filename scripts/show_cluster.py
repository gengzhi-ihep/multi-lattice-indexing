
"""
 (c) IHEP 2021. All rights reserved.
 Author: Geng Zhi
 Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS).
 If you have any problem with the program, please contact author with the following
 email: gengz@ihep.ac.cn
"""

import scipy.cluster.hierarchy as sch
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import optparse

if __name__ == '__main__':

   usage="python %prog [options]"
   parser = optparse.OptionParser(usage)

   parser.add_option('-f','--file',dest='file',default="class.txt",type='string',help='File list containing 3 orientation angles as well as match-rates (4 columns)')
   parser.add_option('-c','--cut',dest='mcut',default=0.2,type='float',help='Cutoff to select higher match-rate')
   parser.add_option('-t','--tall',dest='tall',default=30,type='float',help='Threshold to cluster neighbor points')

   (opts, args) = parser.parse_args()

   data = np.loadtxt(opts.file)

   sel = np.where(data[:,3]> opts.mcut)[0]

   subdata = data[sel]

   x = np.sin(subdata[:,0])*np.cos(subdata[:,1])*subdata[:,2]/np.pi*180
   y = np.sin(subdata[:,0])*np.sin(subdata[:,1])*subdata[:,2]/np.pi*180
   z = np.cos(subdata[:,0])*subdata[:,2]/np.pi*180
   h = subdata[:,3]

   xyz = np.vstack([x,y,z]).T

   fig = plt.figure()

   ax = fig.add_subplot(111,projection='3d')

   dismat = sch.distance.pdist(xyz,'euclidean')

   o = sch.linkage(dismat,method='average')

   r = sch.cut_tree(o, height=opts.tall)

#   np.savetxt('out.tmp',r)

   clusters = []
   for i in np.unique(r):
       clusters.append(xyz[r[:,0] == i].mean(axis=0))

   classes = np.array(clusters).reshape(-1,3)

   np.savetxt('clusters.txt',classes)
 
   opts = {}

   for i in range(np.unique(r).shape[0]):
       rbg = np.random.randint(0,256,size=3)/256.0
       opts[i] = list(rbg)

   labels = []
   for i in range(r.shape[0]):
      labels.append(r[i,0])

   colors = list(map(lambda x: opts[x], labels))

   ax.scatter(x,y,z,s=10,c=colors,marker='.')

   colors2 = list(map(lambda x: opts[x], range(len(opts))))

   ax.scatter(classes[:,0],classes[:,1],classes[:,2], s=120, c=colors2,alpha=0.2,marker='o')

   real = np.array([[22.8,70.31,266.71],[52.39,327.14,194.33],[37.1,235.35,196.76],[27.34,330.1,192.75],[37.47,0.53,201.41]])/180.0*np.pi
   realx = np.sin(real[:,0])*np.cos(real[:,1])*real[:,2]/np.pi*180
   realy = np.sin(real[:,0])*np.sin(real[:,1])*real[:,2]/np.pi*180
   realz = np.cos(real[:,0])*real[:,2]/np.pi*180

   ax.scatter(realx,realy,realz, s=100, c='r',alpha=0.6,marker='x')

   plt.show()

#   d = sch.dendrogram(z)

#   plt.show()

