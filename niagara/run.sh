#!/bin/bash

#-------------COMMENTS-------------
# A bash script to run the Dedalus python code
#----------------------------------

rm *~

echo "Running Dedalus script"
time mpirun  python3 simulation.py

echo "Merging outputs"
time mpirun  python3 merge.py snapshots

echo "Done"

exit
