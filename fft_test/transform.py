"""
Plot for temporal power spectrum of total bouyancy (e.g z = btot(x,t) ) in hydrostatic
and nonhdrostatic simulation. It is important to first plot the bouyancy contours with
runfig.sh to see they are the right ones.
Imports plot_test_fft. Make sure to copy any hardcoded changes in plot_test into
plot_test_fft so that we are using the same contours.
Also can get it to print and analyze the time series intervals, which is currently
not constant.

Usage:
    plot_2d_series.py <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames]

"""
import numpy as np
import h5py
# import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from dedalus.extras import plot_tools
import plot_test_fft as plot_test #edited plotbot, used for setting ylim
from mpl_toolkits.axes_grid1 import make_axes_locatable #for plotting colorbar

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    tasks = ['btot','btot_stat'] #for bouyancy in non-hydrostat and hydrostat cases
    """much of these settings are probably not needed, but I didn't want to break the code
    by deleting them
    """
    scale = 2.5
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_{:06}.png'.format(write)
    # Layout
    nrows, ncols = 1, 1
    image = plot_tools.Box(len(tasks), 1)
    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)
    '------------------------------------------------------'
    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        time = np.zeros(0)
        diff = np.zeros(0)
        for index in range(start, start+count):
            if index == 400: #for testing purposes
                break
            for n, task in enumerate(tasks):
                if n ==0: #zoomed bouyancy on top, would like to plot contours aswell
                    # Build subfigure axes
                    i, j = divmod(n, ncols)
                    axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                    # Call 3D plotting helper, slicing in time
                    dset = file['tasks'][task]
                    ctour = plot_test.plot_bot_3d(dset, 0, index, k= 1, axes=axes, title='b in Pycnocline (Non-Hydrostat)',
                                           even_scale=True, clim=(3,5)) #get contours for Fourier transform
                    # print(task)
                    if index != 0: #bypass first snapshot, since 1st is 0 for all variables
                        tour = ctour.collections[4].get_paths()[0] #gets 4th contour
                        v = tour.vertices
                        if index == 1: #only need to do this once
                            Z = v[:,1] #first set of z values in bouyancy
                            x = v[:,0] #gets x values which are constant
                        if index > 1:
                            z = v[:,1]
                            Z = np.vstack((Z,z)) #build array of z values
                        if index%5 == 0: #just to know where it is at
                            print('n=', n ,'index=', index) #nth variable, ith timestep
                        time= np.append(time,file['scales/sim_time'][index])
                        if len(time) >= 2:
                            difference = time[index-1] - time[index-2] #time series interval
                            diff = np.append(diff, difference) #make an array
                if n == 1: #zoomed bouyancy on top, would like to plot contours aswell
                    # Build subfigure axes
                    i, j = divmod(n, ncols)
                    axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                    # Call 3D plotting helper, slicing in time
                    dset = file['tasks'][task]
                    ctour = plot_test.plot_bot_3d(dset, 0, index, k= 1, axes=axes, title='b in Pycnocline (Non-Hydrostat)',
                                           even_scale=True, clim=(3,5)) #get contours for Fourier transform
                    # print(task)
                    if index != 0: #bypass first snapshot, since 1st is 0 for all variables
                        tour = ctour.collections[4].get_paths()[0] #now we are on 2nd etc snapshots
                        v = tour.vertices
                        if index == 1: #only need to do this once
                            Z_stat = v[:,1]
                        if index > 1:
                            z_stat = v[:,1]
                            Z_stat = np.vstack((Z_stat,z_stat)) #build array of z values
                        if index%1 == 0: #just to know where it is at
                            print('n=', n ,'index=', index)

                        # plt.plot(x,z)
                        # plt.ylim(-0.15-1e-4, -0.15+1e-4)
            #     elif n == 2: #zoomed bouyancy on top
            #         # Build subfigure axes
            #         i, j = divmod(n, ncols)
            #         axes = mfig.add_axes(i, j, [0, 0, 1, 1])
            #         # Call 3D plotting helper, slicing in time
            #         dset = file['tasks'][task]
            #         plot_test.plot_bot_3d(dset, 0, index, k=0, axes=axes, title='Change in b in Pycnocline ( Non-Hydrostat)',
            #                                even_scale=True, clim=(-5e-5,5e-5))
            #     elif n == 3: #zoomed bouyancy on top
            #         # Build subfigure axes
            #         i, j = divmod(n, ncols)
            #         axes = mfig.add_axes(i, j, [0, 0, 1, 1])
            #         # Call 3D plotting helper, slicing in time
            #         dset = file['tasks'][task]
            #         plot_test.plot_bot_3d(dset, 0, index, k=0, axes=axes, title='Change in b_stat in Pycnocline (Hydrostat)',
            #                                even_scale=True, clim=(-5e-5,5e-5))
            # Add time title
            # title = title_func(file['scales/sim_time'][index])
            # title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            # fig.suptitle(title, fontsize='large')#, x=0.48, y=title_height, ha='left')
            # Save figure
            # savename = savename_func(file['scales/write_number'][index])
            # savepath = output.joinpath(savename)
            # fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)
    ##colour plot for x,t,z where z is the dependant variable, just used for testing
    # X, T = np.meshgrid(x,time) #x,t mesh for z plot

    # fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True) #2 plots
    # im = ax0.pcolormesh(X, T, Z)
    # ax0.set_title('Nonhydrostatic')
    # ax0.set_ylabel('Time (s)')
    #
    # im = ax1.pcolormesh(X, T, Z_stat)
    # ax1.set_title('Hydrostatic')
    # ax0.set_xlabel('Horizontal Position (m)')
    # colorbar(im) #to fit in colorbar
    # plt.tight_layout(h_pad=1)
    # plt.savefig('Colorbar_graph2')
    # fig.clear()
    # exit()

    #Power Spectrum
    #nonhydrostatic
    Z = Z[30:, :] #delete first 30 rows, transients
    Z_c = np.fft.fft(Z, axis = 0) #coeff matrix
    Z_p = np.abs(Z_c)**2 #power matrix
    #hydrostatic
    Z_stat = Z_stat[30:, :] #delete first n rows, transients
    Z_c_stat = np.fft.fft(Z, axis = 0) #coeff matrix
    Z_p_stat = np.abs(Z_c)**2 #power matrix

    n = Z.shape[0] #length of time axis
    timestep = time[100]-time[99] #sampling interval, picked random interval
    freq = np.fft.fftfreq(n, d = timestep) #frequency axis

    X, w = np.meshgrid(x,freq)
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True) #2 plots
    im = ax0.pcolormesh(X, w, Z_p, vmin = 1e-10, vmax=1e-8)
    ax0.set_title('Nonhydrostatic')
    ax0.set_ylabel('Freq (Hz)')

    im = ax1.pcolormesh(X, w, Z_p_stat, vmin = 0, vmax=1e-8)
    ax1.set_title('Hydrostatic')
    ax0.set_xlabel('Horizontal Position (m)')
    colorbar(im) #to fit in colorbar
    plt.tight_layout(h_pad=1)
    plt.savefig('Powerspectrum_graph3')

    print(np.std(diff[10:,])*100/diff[100]) #looking at std dev time intervals as fraction of a later interval
    print(len(time)) # look at number of time series
    fig.clear()




if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    #from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path)
