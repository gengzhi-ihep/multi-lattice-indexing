python 1peak-search-xds.py -p "/home/gengzhi/software/multi-indexing/examples/SR" -f "lyso-3crystals01.????" -n 1 -r "1  0  0"  -m "./result"

cp ./result/SPOT.TXT .

#python 2run-index-mpi.py -k 9 -n 5 -e "/home/gengzhi/software/multi-indexing/bin/index-mpi" -d "/usr/lib64/openmpi/bin/mpirun" -m "."
python 2run-index-omp.py -n 5 -e "/home/gengzhi/software/multi-indexing/bin/index-omp" -m "."
