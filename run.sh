#!/bin/bash

#-------------COMMENTS-------------
#runs script, proccesses resultant HDF5 file as well
#----------------------------------

# A bash script to run the Dedalus python code

rm *~
day=$(date +'%m'-'%d') #date in mm-dd format
now=$(date +'%R') #time in 24hr clock

if [ -e ./snapshots ]
then
    echo "Removing previous snapshots"
    rm -r ./snapshots
fi

# #option for storing frames folder in timestamped repository
# if [ ! -e ./simulations ]
# then
#     echo "Making simulation folder"
#     mkdir ./simulations
# fi
#
# if [ ! -e ./simulations/$day ]
# then
#     echo "Making day folder"
#     mkdir ./simulations/$day
# fi
#
# echo "Making time folder"
# mkdir ./simulations/$day/$now


echo "Running Dedalus script"
time mpiexec -n 4 python3.7 simulation.py

echo "Merging outputs"
time mpiexec -n 4 python3.7 merge.py snapshots

echo "Printing figures"
rm -r frames
time mpiexec -n 4 python3.7 plot_2d_series.py snapshots/*.h5

cd frames
ffmpeg -framerate 10 -i write_%06d.png -c:v libx264 -pix_fmt yuv420p movie.mp4

# cd ..
# cp -r frames simulation.py ./simulations/$day/$now

echo "Done"

exit
