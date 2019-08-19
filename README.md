# How to use

There are 3 main parts to this repository: the Niagara directory, fft_test directory, and the scripts outside of these directories. First be familiar with these outside scripts, then look into those directories.

## Scripts outside of directories

These are the main scripts.
1. simulation.py is the Python script for the actual simulation. It outputs the data in the form of HDF5 files in 'snapshots' folder.
2. merge.py merges HDF5 files into one global dataset. This is provided by the Dedalus library
3. plot_2d_series.py is the routine for plotting the HDF5 files based off
https://bitbucket.org/dedalus-project/dedalus/src/tip/examples/ivp/3d_rayleigh_benard/plot_slices.py?at=default
calling upon plot_test.py for plotting tools. The figures are created in a 'frames' folder
4. plot_test.py is an edited version of the plotting tools provided by
https://bitbucket.org/dedalus-project/dedalus/src/tip/dedalus/extras/plot_tools.py?at=default
This allows us to be more specific with how we want to plot our snapshots.
5. run.sh runs the simulation and plots the resultant files
6. runfig.sh only plots the HDF5 file, for testing the figures
