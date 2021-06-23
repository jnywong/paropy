'''
convective_power.py

Calculate the convective power as a function of radius.

PARODY-JA power.runid computes Ra_parody/Ekman* 1/(volume of the shell) * int(u_r * r/r_o * (C perturbation) * r^2 dt dtheta sin(theta) dphi)
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import h5py
from numpy.lib.function_base import trapz

from paropy.coreproperties import icb_radius, cmb_radius
from paropy.data_utils import load_dimensionless
from paropy.plot_utils import streamfunction, C_shift, merid_outline
from paropy.routines import convective_power_timeavg, ref_codensity

# matplotlib.use('Agg')  # backend for no display

#%% INPUT PARAMETERS
# run_ID  = 'ref_c'
# run_ID = 'd_0_55a'
# run_ID = 'd_0_6a'
# run_ID = 'd_0_65a'
run_ID = 'c-200a'
# run_ID = 'd_0_75a'
# run_ID = 'd_0_8a'
# path containing runs
directory = '/data/geodynamo/wongj/Work/{}'.format(run_ID)
# directory = '/Volumes/NAS/ipgp/Work/{}'.format(run_ID)

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/convective_power'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/convective_power'

#%%----------------------------------------------------------------------------
# Load data
(NR, Ek, Ra, Pr, Pm, fi, rf) = load_dimensionless(run_ID, directory)
if not os.path.exists('{}/convective_power'.format(directory)):
    (radius, Vr, T) = convective_power_timeavg(run_ID, directory)  # timeavg
else:  # load timeavg data
    print('Loading {}/convective_power'.format(directory))
    f = h5py.File('{}/convective_power'.format(directory), 'r')
    for key in f.keys():
        globals()[key] = np.array(f[key])

# FIX
Tref = ref_codensity(radius, rf, fi)
Tvar = T-Tref # remove reference state

ri = radius[0]
ro = radius[-1]
shell_volume = 4*np.pi*(ro**3 - ri**3)/3
check = Ra*trapz(Vr*T*radius**3, radius)/(Ek*shell_volume*ro)
#%%----------------------------------------------------------------------------
# Plot
w, h = plt.figaspect(fig_aspect)
fig1, (ax1a,ax1b) = plt.subplots(2, 1, figsize=(1.5*w, h),sharex=True)
ax1a.plot(radius, Vr)
ax1b.plot(radius, Tvar)

fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(1.5*w, h), sharex=True)
ax2a.plot(radius, T)
ax2b.plot(radius, Tref)