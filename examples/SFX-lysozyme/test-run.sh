../../bin/prepare-sfx

less lysozyme-mosflm-1000.stream | grep "Image filename" >files-mcdps.lst

python run-sfx-index.py -k 9 -n 5 -t 1000 -e "/home/gengzhi/software/multi-indexing/bin/index-sfx" -d "/usr/lib64/openmpi/bin/mpirun" -m "."
