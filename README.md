# How to use

There are 3 main parts to this repository: the Niagara directory, fft_test directory, and the scripts outside of these directories. First be familiar with these outside scripts, then look into those directories.
Descriptions of the scripts are below. An example of a typical routine in using these scripts follows.
## Scripts outside of directories

These are the main scripts, usually all run on NICOGWS2015.
1. `simulation.py` is the Python script for the actual simulation. It outputs the data in the form of HDF5 files in 'snapshots' folder.
2. `merge.py` merges HDF5 files into one global dataset. This is provided by the Dedalus library
3. `plot_2d_series.py` is the routine for plotting the HDF5 files based off
https://bitbucket.org/dedalus-project/dedalus/src/tip/examples/ivp/3d_rayleigh_benard/plot_slices.py?at=default
,calling upon plot_test.py for plotting tools. The figures are created in a 'frames' folder
4. `plot_test.py` is an edited version of the plotting tools provided by
https://bitbucket.org/dedalus-project/dedalus/src/tip/dedalus/extras/plot_tools.py?at=default . This allows us to be more specific with how we want to plot our snapshots.
5. `run.sh` runs the simulation and plots the resultant files
6. `runfig.sh` only plots the HDF5 file, for testing the figures

## fft_test directory
This is used for producing the power spectrum for (total) buoyancy contours in the pycnocline. Some example pictures are provided.
1. `transform.py` takes ../snapshots/\*.hf file in the outside directory and creates power spectrum by calling on `plot_test_fft.py`
2. `plot_test_fft.py` extracts the buoyancy contours of interest in each snapshot. It is essentially plot_test.py without plot formatting
3. `runtransfrom.sh` runs `transform.py`

## Niagara directory
For running on Niagara. No plotting is done on Niagara, we only get the HDF5 files from it and then plot locally.
1. `simulation.py` is the exact same as above, with changes to accompany the use of 40 cores
2. `run.sh` is the same as above but running only the simulations
3. `run.slrm` is submitted to run on Niagara. It runs `run.sh`. See https://github.com/ngrisouard/dedalus-on-niagara for details on this.

# In practice
In order to change anything in the simulation (e.g. pycnocline height), I must edit `simple.py`. To specify the variables I want to plot and general plot formatting such as titles, I edit `plot_2d_series.py`.

If I want to plot buoyancy contours, I must edit `plot_test.py`. There is also the option to set the ylim as well as more detailed plot formatting by editing this file. Note that I must hard-code the N^2 stratifaction function in this file in order to integrate and find buoyancy contours at specified levels. Thus any changes in stratification in `simulation.py` must be added to `plot_test.py`.

I submit `run.sh` on the workstation and get the figures. On Niagara, I run `run.slrm` and transfer the HDF5 files to the workstation, where I run `runfig.sh` to get figures. If I want different figures for the same simulation, I edit the plotting files and run `runfig.sh`.

The task of finding buoyancy contours is important to note. In editing `plot_test.py`, I need to specify the buoyancy value range and number of contours at the desired height, which requires calculating. Furthermore, the default setting is that the contour values are equally spaced in the specified range. Since N^2 is not linear, this physical spacing in plot is not linear. In addition, the nonlinearity of N^s also makes specifying the ylim for the plot difficult. I have not come up with a systematic way of doing this yet - it has been more of trial and error in setting the buoyancy range and ylim.

After setting the desired buoyancy contours, I use the fft_test directory for a power spectrum. I must hard code the same stratification and buoyancy interval settings in the contour part of the `plot_test.py` script into that of `plot_test_fft.py`. Then I pick the one specific contour I want by specifying its index, and run `runtransform.py`.

## Possible improvements for the future
Finding the buoyancy contours can be an inconvenient task - one can look for a more systematic routine for this. Furthermore, it may be more convenient to create a global N^2 function in an individual script such that when N^2 is changed, it only need be changed in one script. 
