"""
Simulation script for Niagara. Everything is the same as simulation.py except
the dedalus modes are changed for 40 cores, and some parameters are changed.
I believe the ony the amplitude and snapshot intervals are changed here.

Added more variables so that both hydro static and non-hydrostatic
are solved in this srcipt. Variables are then used to find the
difference.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 merge.py snapshots
    $ mpiexec -n 4 python3 plot_2d_series.py snapshots/*.h5

"""

from docopt import docopt
import logging
import time

# import multiprocessing
# from joblib import Parallel, delayed

import numpy as np
from mpi4py import MPI
import scipy as sp
from math import isnan

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.extras.plot_tools import plot_bot_2d
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


nu = 1.e-6 #6 [m2/s] visco
kk = 1.e-9 # [m2/s] kappa

mx = 640  # number of modes in x
mz = mx//2  # number of modes in z  NG reduced
H = .5  # [m]  domain height
lfactor = 10 #factor multiplying domain length
L = lfactor*6.*H  # [m]  domain length
# bigK = 4.  # [rad/m]
num = 12  # no. horizontal wl for axis
theta = -0.95*np.pi/2  # [rad] angle of beam w.r.t. horizontal
wall = 115.  # [minutes] walltime_max

Lx = L/num  # [m] horizontal wavelength
rho0 = 1 # [kg/m^3] background density
N0 = 0.6 # rad/s
g = 9.81  # m/s^2 gravity
kx = 2*np.pi/Lx
# kx = bigK * np.cos(theta)  # individual wave numbers
kz = kx * np.tan(theta)
# L = 5.#horizontal length (if not using kx)
omega = N0*np.cos(theta)  # freq, according to dispersion relation
T = 2*np.pi/omega
A = 1e-4
PolRel = {'u': -A*g*omega*np.tan(theta)/N0**2,  # dictionary of forcing coeffs
          'w': A*g*omega/N0**2,
          'b': -A*g}

for fld in PolRel:
    if fld == 'b':
        units = 'm/s2'
    else:
        units = 'm/s'
    logger.info('amplitude of {0} = {1:.2e} {2}'.format(
        fld, PolRel[fld], units))
logger.info('frequency is {0:.2e} 1/s'.format(omega))
'-----------create basis and domain---------'
x_basis = de.Fourier('x', mx, interval=(0., L), dealias=3/2)

# n1= 2**10#layers' resolution
# n2= 2**7
# n3= 2**6

# zb1 = de.Chebyshev('z1', n1, interval=(-H,-H/3))
# zb2 = de.Chebyshev('z2', n2, interval=(-H/3,-2*H/3))
# zb3 = de.Chebyshev('z3', n3, interval=(-2*H/3,0))
# z_basis = de.Compound('z', (zb1, zb2, zb3))


z_basis = de.Chebyshev('z', mz, interval=(-H, 0.), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
x = domain.grid(0)
z = domain.grid(1)

'-----------PDE parameters and variables---------'
# V = (u,v)
pb = de.IVP(domain, variables=['u', 'uz', 'w',
                               'wz', 'b', 'bz', 'p', 'btot',
                               'u_stat', 'uz_stat', 'w_stat', #exact same variables for hydrostat situation
                               'wz_stat', 'b_stat', 'bz_stat', 'p_stat', 'btot_stat',
                                'u_diff', 'w_diff', #difference between two correspondin variables
                                'b_diff', 'p_diff'])
# NG: something I did not know until recently: all variables are Dirichlet
# by default.
pb.meta['wz', 'bz', 'u', 'p','wz_stat', 'bz_stat', 'u_stat', 'p_stat']['z']['dirichlet'] = False
pb.parameters['nu'] = nu
pb.parameters['kk'] = kk
pb.parameters['kz'] = kz
pb.parameters['kz'] = kz
pb.parameters['lfactor'] = lfactor
pb.parameters['r0'] = rho0

def fill_N2(z):  # non-const Stratification (cts)
    N2 = 0*z
    hwidth = H/80 #half-width of bump
    term = ((z+H/40)/hwidth)**2 #centred at -H/20
    N2[z > -3*H/80] = np.e*4**2* \
        np.exp(-1/(1-term[z > -3*H/80])) + N0**2  # bump function
    N2[z <= -3*H/80] = N0**2
    N2[z >= -H/80] = N0**2
    return N2


# x, z = domain.grids(scales=1)
N2 = domain.new_field()
N2.meta['x']['constant'] = True  # means the NCC is constant along x
N2['g'] =fill_N2(z)
# N2['g'] = N0**2

# from dedalus.extras.plot_tools import plot_bot_2d
# plot_bot_2d(N2);
# plt.savefig('Brunt-Vais freq')
# plt.show()
#
#
# plt.plot(z[0],fill_N2(z[0]))
# plt.title('Brunt-Vais Freq')
# plt.xlabel('z')
# plt.ylabel('N^2')
# plt.savefig('Brunt-Vais Profile')
# exit() #for showing strat


pb.parameters['N2'] = N2  # pass function in as a parameter

# def sponge(z):
#     """sponge layer on the top, bump centred at H DOES NOT WORK DUE TO MEMORY ERROR
#     """
#     y= 0*z
#     bignu = 1e2-nu #peak Viscosity
#     hwidth = H/200 #half-width of bump as fraction of H
#     term = (z)/hwidth
#     y[z <= -hwidth] = nu
#     y[z > -hwidth] = np.e*bignu*np.exp(-1/(1-term[z > -hwidth]**2)) + nu
#     return y
#
# spnu = domain.new_field()
# spnu.meta['x']['constant'] = True  # means the NCC is constant along x
# spnu['g'] = sponge(z)
#
# # plt.plot(z[0],sponge(z[0]))
# # plt.title('Sponge')
# # plt.xlabel('z')
# # plt.ylabel('Viscocity')
# # plt.savefig('sponge_profile')
# # exit() #for showing sponge profile
# pb.parameters['spnu'] = spnu

'---------Boundary forcing----------'


def local1(z):
    """used for localized forcing on bottom, hard
    """
    y = 0*z
    y[z < Lx] = 0
    y[z >= Lx] = 1
    y[z > 2*Lx] = 0
    return y


def local2(z):
    """used for localized forcing on bottom, soft
       centres on 3/2*Lx, window of Lx
    """
    y = 0*z
    term = (2/Lx)*(z-3*Lx/2)
    y[z <= Lx] = 0
    y[z > Lx] = np.e*np.exp(-1/(1-term[z > Lx]**2))
    y[z >= 2*Lx] = 0
    return y

def local3(z):
    """used for localized forcing on bottom, soft
       centres on 3/2*Lx, window of Lx, using sine
    """
    y = 0*z
    y[z <= Lx] = 0
    term = (z-Lx)*np.pi/Lx
    y[z > Lx] = np.sin(term[z > Lx])
    y[z >= 2*Lx] = 0
    return y


def timewindow(t):
    """time windowing: a is start of decreaase, b is when
    it returns 0. It is a C^0 linear fn
    """
    a = 0.2
    b = 0.3
    const = 1 - a/(a-b)
    if t <= a:
        r=1
    elif t > a and t <=b:
        r = t/(a-b) + const
    else:
        r=0
    return r

def timewindowincrease(t):
    """time windowing: start at t=0 with value 0, linearly increase
	to 1 after t=a. Stay at 1 afterwards
	It is a C^0 linear fn
    """
    a = T/5
    if t >= a:
        r=1
    elif t < a:
        r=t/a
    return r

def BoundaryForcingS(*args):
    """this function applies its arguments and returns the forcing, sine
    NG changed name
    """
    t = args[0].value  # this is a scalar; we use .value to get its value
    x = args[1].data  # this is an array; we use .data to get its values
    z = args[2].data
    return np.sin(kx*x+kz*z-omega*t)  * local3(x) *timewindowincrease(t)


def ForcingS(*args, domain=domain, F=BoundaryForcingS):
    """This function takes arguments *args, a function F, and a domain and
    returns a Dedalus GeneralFunction that can be applied.
    NG changed name
    """
    return de.operators.GeneralFunction(domain, layout='g', func=F, args=args)


def BoundaryForcingC(*args):  # duplicate for cos in b forcing
    """this function applies its arguments and returns the forcing, cosine
    NG changed name
    """
    t = args[0].value  # this is a scalar; we use .value to get its value
    x = args[1].data  # this is an array; we use .data to get its values
    z = args[2].data
    return np.cos(kx*x+kz*z-omega*t)   *local3(x) *timewindowincrease(t)


def ForcingC(*args, domain=domain, F=BoundaryForcingC):
    """This function takes arguments *args, a function F, and a domain and
    returns a Dedalus GeneralFunction that can be applied.
    NG changed name
    """
    return de.operators.GeneralFunction(domain, layout='g', func=F, args=args)


for fld in ['u', 'w', 'b']:
    BF = domain.new_field()
    BF['g'] = PolRel[fld]
    pb.parameters['cf' + fld] = BF  # pass function in as a parameter.


# now we make it parseable, so the symbol BF can be used in equations
# and the parser will know what it is.
de.operators.parseables['BFS'] = ForcingS
de.operators.parseables['BFC'] = ForcingC
# '------------PDES----------------'

# '------------NON HYDROSTAT----------------'
pb.substitutions['vdg(A, Az)'] = "- u*dx(A) - w*Az"  # RHS
pb.substitutions['diff(A, Az)'] = "- dx(dx(A)) - dz(Az)"  # LHS

# '------------HYDROSTAT----------------'
pb.substitutions['vdg_stat(A, Az)'] = "- u_stat*dx(A) - w_stat*Az"  # RHS
pb.substitutions['diff_stat(A, Az)'] = "- dx(dx(A)) - dz(Az)"  # LHS

# '------------NON HYDROSTAT----------------'
pb.add_equation("dx(u) + wz = 0")  # continuity
pb.add_equation("dt(b) + N2*w + kk*diff(b, bz)"  # bouyancy
                + "= vdg(b, bz)")
pb.add_equation("dt(u) + dx(p) + nu*diff(u, uz)"  # momentum in x
                + "= vdg(u, uz)")
pb.add_equation("dt(w) + dz(p)/r0 - b + nu*diff(w, wz)"  # momentum in y
                + "= vdg(w, wz)")
pb.add_equation("dz(btot) - bz"  # equation for total bouyancy
                + "= - N2")
pb.add_equation("bz - dz(b) = 0")
pb.add_equation("uz - dz(u) = 0")
pb.add_equation("wz - dz(w) = 0")

# boundary conditions in z
pb.add_bc("left(b) = left(cfb*BFC(t,x,z))")
pb.add_bc("left(uz) = left(kz*cfu*BFC(t,x,z))")
pb.add_bc("left(w) = left(cfw*BFS(t,x,z))")  # forcing
pb.add_bc("right(b) = 0")
pb.add_bc("right(btot) = 0")
pb.add_bc("right(w) = 0", condition="(nx != 0)")
pb.add_bc("left(p) = 0", condition="(nx == 0)")
# pb.add_bc("right(w) = 0")#,condition="(nx != 0)")
# pb.add_bc("left(p) = 0")
pb.add_bc("right(uz) = 0")
# pb.add_bc("right(w) = 0")
# pb.add_bc("right(p) = 0", condition="(nx == 0)")  # not sure about that
# pb.add_bc("integ_z(p) = 0")

# '------------HYDROSTAT----------------'
pb.add_equation("dx(u_stat) + wz_stat = 0")  # continuity
pb.add_equation("dt(b_stat) + N2*w_stat + kk*diff_stat(b_stat, bz_stat)"  # bouyancy
                + "= vdg_stat(b_stat, bz_stat)")
pb.add_equation("dt(u_stat) + dx(p_stat) + nu*diff_stat(u_stat, uz_stat)"  # momentum in x
                + "= vdg_stat(u_stat, uz_stat)")
pb.add_equation("dz(p_stat) - r0*b_stat"  # momentum in y
                + "= 0")
pb.add_equation("dz(btot_stat) - bz_stat"  # equation for total bouyancy
                + "= - N2")
pb.add_equation("bz_stat - dz(b_stat) = 0")
pb.add_equation("uz_stat - dz(u_stat) = 0")
pb.add_equation("wz_stat - dz(w_stat) = 0")

# boundary conditions in z
pb.add_bc("left(b_stat) = left(cfb*BFC(t,x,z))")
pb.add_bc("left(uz_stat) = left(kz*cfu*BFC(t,x,z))")
pb.add_bc("left(w_stat) = left(cfw*BFS(t,x,z))")  # forcing
pb.add_bc("right(b_stat) = 0")
pb.add_bc("right(btot_stat) = 0")
pb.add_bc("right(w_stat) = 0", condition="(nx != 0)")
pb.add_bc("left(p_stat) = 0", condition="(nx == 0)")
# pb.add_bc("right(w) = 0")#,condition="(nx != 0)")
# pb.add_bc("left(p) = 0")
pb.add_bc("right(uz_stat) = 0")
# pb.add_bc("right(w) = 0")
# pb.add_bc("right(p) = 0", condition="(nx == 0)")  # not sure about that
# pb.add_bc("integ_z(p) = 0")
# '------------DIFFERENCE EQUATIONS----------------'
pb.add_equation("u_diff = u - u_stat")
pb.add_equation("w_diff = w - w_stat")
pb.add_equation("p_diff = p - p_stat")
pb.add_equation("b_diff = b - b_stat")
# '---------Build solver ----------'
solver = pb.build_solver(de.timesteppers.RK443)
logger.info('Solver built')
#
'------------Initial Bouyancy------'
# N2_new = domain.new_field()
# # N2_new.meta['x']['constant'] = True  # means the NCC is constant along x
# N2_new['g'] =fill_N2(z)
#
# btot = solver.state['btot']
# btot_stat = solver.state['btot_stat']
# #
# N2_new.integrate('z', out=btot)
# N2_new.integrate('z', out=btot_stat)
# # plot_bot_2d(b)

num_per = 4
# Initial timestep

# NG changed things here
c_px = omega/kx  # [m/s] horizontal phase speed
c_pz = omega/kz  # [m/s] vertical phase speed


# NG changed things here
# make sure timestep is small enough
dt = abs(0.5*min(L/(mx*c_px), H/(mz*c_pz)))/3
# dt = 1e-1*T

# Integration parameters
solver.stop_sim_time = num_per*T
solver.stop_wall_time = wall*60.
solver.stop_iteration = np.inf

# Analysis
# NG changed output frequency
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=T/200,
                                              max_writes=1000)

snapshots.add_system(solver.state)

# CFL
# NG: I changed bits and pieces here and there.
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.8,
                     max_change=1.5, min_change=0.5,  max_dt=T/200,
                     threshold=0.1)
CFL.add_velocities(('u', 'w'))


# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=100)
flow.add_property("4*(bz+N2) - uz**2", name='Ri_red')

start_time = time.time()
while solver.ok:
    dt = CFL.compute_dt()
    dt = solver.step(dt)
    if (solver.iteration) % 5 == 0:
        logger.info('Iteration: {0:d}, Time: {1:e}T, dt: {2:e}T'.format(
            solver.iteration, solver.sim_time/T, dt/T))
        # logger.info('Min Ri_red = {0:f}'.format(flow.min('Ri_red')))
        if np.isnan(flow.min('Ri_red')).any() \
           or np.isinf(flow.min('Ri_red')).any():
            raise NameError('Ri blew up')
end_time = time.time()
# Print statistics
logger.info('Iterations: {0:d}'.format(solver.iteration))
logger.info('Sim end time: {0:f}'.format(solver.sim_time))
logger.info('Periods run: {0:f}'.format(solver.sim_time/T))
logger.info('Run time: {0:.2f} sec'.format(end_time - start_time))
logger.info('Run time: {0:f} cpu-hr'.format((end_time -
                                             start_time)/3600*domain.dist.comm_cart.size))
