# How to use

There are 3 main parts to this repository: the Niagara directory, fft_test directory, and the scripts outside of these directories. First be familiar with these outside scripts, then look into those directories.
Descriptions of the scripts are below. An example of a typical routine in using these scripts follows.
## Scripts outside of directories

These are the main scripts, usually all run on NICOGWS2015.
1. `simulation.py` is the Python script for the actual simulation. It outputs the data in the form of HDF5 files in `snapshots` folder.
2. `merge.py` merges HDF5 files into one global dataset. This is provided by the Dedalus library
3. `plot_2d_series.py` is the routine for plotting the HDF5 files based off
https://bitbucket.org/dedalus-project/dedalus/src/tip/examples/ivp/3d_rayleigh_benard/plot_slices.py?at=default
,calling upon plot_test.py for plotting tools. The figures are created in a `frames` folder
4. `plot_test.py` is an edited version of the plotting tools provided by
https://bitbucket.org/dedalus-project/dedalus/src/tip/dedalus/extras/plot_tools.py?at=default . This allows us to be more specific with how we want to plot our snapshots.
5. `run.sh` runs the simulation and plots the resultant files
6. `runfig.sh` only plots the HDF5 file, for testing the figures

## fft_test directory
This is used for producing the power spectrum for (total) bouyancy contours in the pycnocline. Some example pictures are provided.
1. `transform.py` takes ../snapshots/\*.hf file in the outside directory and creates power spectrum by calling on `plot_test_fft.py`
2. `plot_test_fft.py` extracts the bouyancy contours of interest in each snapshot. It is essentially plot_test.py without plot formatting
3. `runtransfrom.sh` runs `transform.py`

## Niagara directory
For running on Niagara. No plotting is done on Niagara, we only get the HDF5 files from it and then plot locally.
1. `simulation.py` is the exact same as above, with changes to accompany the use of 40 cores
2. `run.sh` is the same as above but running only the simulations
3. `run.slrm` is submitted to run on Niagara. It runs `run.sh`. See https://github.com/ngrisouard/dedalus-on-niagara for details on this.

## In practice
