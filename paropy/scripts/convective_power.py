'''
convective_power.py

Calculate the convective power as a function of radius.

PARODY-JA power.runid computes -Ra_parody/Ekman* 1/(volume of the shell) * int(u_r * r/r_o * (C perturbation) * r^2 dt dtheta sin(theta) dphi)
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import h5py
from numpy.lib.function_base import trapz

from paropy.coreproperties import icb_radius, cmb_radius
from paropy.data_utils import load_dimensionless, load_power
from paropy.plot_utils import streamfunction, C_shift, merid_outline
from paropy.routines import convective_power_timeavg, ref_codensity

# matplotlib.use('Agg')  # backend for no display
plt.close('all')
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
    (radius, I) = convective_power_timeavg(run_ID, directory)  # timeavg
else:  # load timeavg data
    print('Loading {}/convective_power'.format(directory))
    f = h5py.File('{}/convective_power'.format(directory), 'r')
    for key in f.keys():
        globals()[key] = np.array(f[key])
# Check with diagnostics
df_power = load_power(run_ID,directory)
check_diag = df_power["available_convective_power_per_unit_vol"].mean()

ro = radius[-1]
convective_power = Ra*I/(Ek*radius[-1])
#%%----------------------------------------------------------------------------
# Plot
w, h = plt.figaspect(fig_aspect)
# fig1, (ax1a,ax1b) = plt.subplots(2, 1, figsize=(1.5*w, h),sharex=True)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
ax1.plot(radius, convective_power)
# ax1b.plot(radius, Tvar)

# fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(1.5*w, h), sharex=True)
# ax2a.plot(radius, T)
# ax2b.plot(radius, Tref)

plt.show()