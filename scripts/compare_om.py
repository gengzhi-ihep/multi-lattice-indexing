import numpy as np

stream_file = 'lysozyme-xgandalf-1000-multi.stream'
files_mcdps = 'files-mcdps.lst'
log_mcdps = 'log-old'

#read in oms from MOSFLM stream file
is_in_image = False
is_indexed = False
oms_mosflm = dict()
npeaks_mosflm = dict()
lines = open(stream_file,'r').readlines()

for line in lines:

    if 'Begin chunk' in line:
        is_in_image = True
        oms = []
    elif 'End chunk' in line:
        is_in_image = False

    if is_in_image:
        if 'Image filename' in line:
            filename = line.split()[-1]
        elif 'num_peaks' in line:
            numpeaks = int(line.split()[-1])
            npeaks_mosflm[filename] = numpeaks
        elif 'Begin crystal' in line:
            is_indexed = True
        elif 'End crystal' in line:
            is_indexed = False

        if is_indexed:

            if 'astar' in line:
                asx = float(line.split()[2])/10
                asy = float(line.split()[3])/10
                asz = float(line.split()[4])/10
            elif 'bstar' in line:
                bsx = float(line.split()[2])/10
                bsy = float(line.split()[3])/10
                bsz = float(line.split()[4])/10
            elif 'cstar' in line:
                csx = float(line.split()[2])/10
                csy = float(line.split()[3])/10
                csz = float(line.split()[4])/10
                om = np.array([[asx,bsx,csx],[asy,bsy,csy],[asz,bsz,csz]])
                oms.append(om)
                oms_mosflm[filename] = oms

print
print "Number of images indexed (MOSFLM): %d"%len(oms_mosflm)
#number of lattices per image
nlatt = []
for key in oms_mosflm.keys():
    nlatt.append(len(oms_mosflm[key]))
print "Number of lattices indexed (MOSFLM): %d"%sum(nlatt)
print "Number of Lattices (MOSFLM):  1: %d, 2: %d, 3: %d, 4: %d, 5: %d\n"%(nlatt.count(1),nlatt.count(2),nlatt.count(3),nlatt.count(4),nlatt.count(5))
        

#read in filenames sequenced from MCDPS indexing
images_seqs = dict()
names=open(files_mcdps,"r").readlines()
for (i,name) in enumerate(names):
    images_seqs[i+1] = name.split()[-1]

#read in oms from MCDPS indexing
oms_mcdps = dict()
pos_mcdps = dict()
lines = open(log_mcdps,'r').readlines()
pos_err = 0
n_latt = 0

for (m,line) in enumerate(lines):

    if "Image" in line:

        oms = []

        num = int(line.split()[1])
        asx = float(lines[m+2].split()[0])
        bsx = float(lines[m+2].split()[1])
        csx = float(lines[m+2].split()[2])
        asy = float(lines[m+3].split()[0])
        bsy = float(lines[m+3].split()[1])
        csy = float(lines[m+3].split()[2])
        asz = float(lines[m+4].split()[0])
        bsz = float(lines[m+4].split()[1])
        csz = float(lines[m+4].split()[2])
        oms.append(np.array([[asx,bsx,csx],[asy,bsy,csy],[asz,bsz,csz]]))
        oms.append(np.array([[bsx,-asx,csx],[bsy,-asy,csy],[bsz,-asz,csz]]))
        oms.append(np.array([[bsx,asx,-csx],[bsy,asy,-csy],[bsz,asz,-csz]]))
        oms.append(np.array([[-asx,bsx,-csx],[-asy,bsy,-csy],[-asz,bsz,-csz]]))
        oms.append(np.array([[-asx,-bsx,csx],[-asy,-bsy,csy],[-asz,-bsz,csz]]))
        oms.append(np.array([[-bsx,-asx,-csx],[-bsy,-asy,-csy],[-bsz,-asz,-csz]]))
        oms.append(np.array([[asx,-bsx,-csx],[asy,-bsy,-csy],[asz,-bsz,-csz]]))
        oms.append(np.array([[-bsx,asx,csx],[-bsy,asy,csy],[-bsz,asz,csz]]))

        if oms_mcdps.has_key(images_seqs[num]):
            oms_tmp = oms_mcdps[images_seqs[num]]
            for om in oms:
                oms_tmp.append(om)
            oms_mcdps[images_seqs[num]] = oms_tmp 
        else:
            oms_mcdps[images_seqs[num]] = oms

    elif "Position Error" in line:
        pos_mcdps[images_seqs[num]] = line
        n_latt += 1
        pos_err += float(line.split()[2])
        
print 'Average positional error (MCDPS): %f \n'%(pos_err/n_latt)
print "Number of images indexed (MCDPS): %d"%len(oms_mcdps)

nlatt = []
for key in oms_mcdps.keys():
    nlatt.append(len(oms_mcdps[key])/8)
print "Number of lattices indexed (MCDPS): %d"%sum(nlatt)
print "Number of Lattices (MCDPS):  1: %d, 2: %d, 3: %d, 4: %d, 5: %d\n"%(nlatt.count(1),nlatt.count(2),nlatt.count(3),nlatt.count(4),nlatt.count(5))


#calculate norm errors between common images
ncommon_imgs = 0
ncommon_latt = 0
nmissed_latt = 0
norms = []
com_files = []
unmatched = []

for key2 in oms_mcdps.keys():

    is_matched = False

    for key1 in oms_mosflm.keys():

        #find common indexed images
        if key1 == key2:
            is_matched = True
            ncommon_imgs += 1
            norm = 1000
            #loop through each lattice for matching
            for om1 in oms_mosflm[key1]:

                for om2 in oms_mcdps[key2]:
                    tmp = np.linalg.norm(np.dot(np.linalg.inv(om1),om2)-np.identity(3))
                    if tmp < norm:
                        norm = tmp
                if norm < 1:
                    norms.append(norm)
                    ncommon_latt += 1
                    com_files.append(key1)

    if not is_matched: 
        unmatched.append("%s %s"%(key2,pos_mcdps[key2]))
#        unmatched.append("%s\n"%(key2))

#write out filenames indexed by mcdps not indexed by mosflm
with open("unmatched_mcdps.lst","w") as fp:
    for line in unmatched:
        fp.write(line)
fp.close()

#write out filenames indexed by mosflm not indexed by mcdps
unmatched = []
for key1 in oms_mosflm.keys():
    is_matched = False
    for key2 in oms_mcdps.keys():
        if key1 == key2:
            is_matched = True
            break
    if not is_matched:
        unmatched.append("%s\n"%key1)
with open("unmatched_mosflm.lst","w") as fp:
    for line in unmatched:
        fp.write(line)
fp.close()

norms = np.array(norms)
#print norms[np.where(norms>0.1)[0]]

print "Number of common images shared by both methods: %d"%ncommon_imgs
print "Number of common lattices shared by both methods: %d"%ncommon_latt
print 'Average OM norm error for all common lattices: %f \n'%np.mean(norms)

#make histogram against peak numbers

hist_mcdps = np.zeros(14)
hist_mosflm = np.zeros(14)

for key in oms_mosflm.keys():
    num = npeaks_mosflm[key]
    if num >= 140:
        hist_mosflm[13] += 1
    elif num < 10:
        hist_mosflm[0] += 1
    else:
        hist_mosflm[num//10] += 1

print "Number of indexed patterns (MOSFLM) against number of peaks with an interval of 10:"
print hist_mosflm

for key in oms_mcdps.keys():
    num = npeaks_mosflm[key]
    if num >= 140:
        hist_mcdps[13] += 1
    elif num < 10:
        hist_mcdps[0] += 1
    else:
        hist_mcdps[num//10] += 1

print "Number of indexed patterns (MCDPS) against number of peaks with an interval of 10:"
print hist_mcdps

hist_norms = np.zeros(14)
hist_count = np.zeros(14)
for (key, norm) in zip(com_files,norms):
    num = npeaks_mosflm[key]
    if num >= 140:
        hist_norms[13] += norm
        hist_count[13] += 1
    elif num < 10:
        hist_norms[0] += norm
        hist_count[0] += 1
    else:
        hist_norms[num//10] += norm
        hist_count[num//10] += 1

avg_norm = np.zeros(14)
n = 0
for (count, norm) in zip(hist_count, hist_norms):
    if count == 0:
        avg_norm[n] = 0
    else:
        avg_norm[n] = norm/count
    n += 1
       
print "Average norm errors against number of peaks with an interval of 10:"
print avg_norm

