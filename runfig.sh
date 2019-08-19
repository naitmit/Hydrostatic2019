#!/bin/bash

#-------------COMMENTS-------------
#Takes desired HDF5 file and runs plotting routine. For
#testing the plotting routine and checking the figures.
#----------------------------------

rm *~
echo "Printing figures"
rm -r frames #be careful, may want to save previous frames folder
time mpiexec -n 4 python3.7 plot_2d_series.py snapshots/*.h5 #--output=frames1

cd frames
ffmpeg -framerate 10 -i write_%06d.png -c:v libx264 -pix_fmt yuv420p movie.mp4

echo "Done"

exit
