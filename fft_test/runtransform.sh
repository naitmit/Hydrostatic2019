#!/bin/bash

#-------------COMMENTS-------------
#For running power spectrum script. Uses HDF5 file in the snapshots folder.
#----------------------------------

# A bash script to run the Dedalus python code

rm *~

time mpiexec -n 1 python3.7 transform.py ../snapshots/snapshots_s1.h5
echo "Done"

exit
