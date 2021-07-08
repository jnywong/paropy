'''
gauss_ratios.py

Compare ratios of Gauss coeffcients, with G3 = g30/g10 and G2 = g20/g10
'''

import os
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from paropy.data_utils import surfaceload, load_dimensionless
from paropy.routines import gauss_ratios_timeavg

# matplotlib.use('Agg')  # backend for no display
plt.close('all')
#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = ['chem_200d', 'ref_c', 'd_0_55a', 'd_0_6a', 'd_0_65b',
          'c-200a', 'd_0_75a', 'd_0_8a']  # PARODY simulation tag
# path containing simulation output
# run_ID = ['c-200a']
dirName = '/data/geodynamo/wongj/Work'

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/gauss_ratios'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/gauss_ratios'
#%% Load data
w, h = plt.figaspect(fig_aspect)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
ax1.axhline(0, linestyle='--', color='k')

rf_out=np.zeros(len(run_ID))
G2_out=np.zeros(len(run_ID))
G3_out=np.zeros(len(run_ID))
i = 0
for run in run_ID:
    directory = '/data/geodynamo/wongj/Work/{}'.format(run)
    (_, _, _, _, _, _, rf) = load_dimensionless(run, directory)
    if not os.path.exists('{}/gauss_ratios_timeavg'.format(directory)):
        g10, g20, g30, G2, G3 = gauss_ratios_timeavg(run, directory)
    else:  # load timeavg data
        print('Loading {}/gauss_ratios_timeavg'.format(directory))
        f = h5py.File('{}/gauss_ratios_timeavg'.format(directory), 'r')
        for key in f.keys():
            globals()[key] = np.array(f[key])
    rf_out[i]=rf
    G2_out[i]=G2
    G3_out[i]=G3
    i+=1

ax1.plot(rf_out[0], G2_out[0], marker='s', markersize = '10', markerfacecolor= 'None', color='k')
ax1.plot(rf_out[0], G3_out[0], marker='8', markersize = '10', markerfacecolor= 'None', color='k')
ax1.plot(rf_out[1], G2_out[1], marker='s', markersize='10', markerfacecolor= 'None', color='darkgrey')
ax1.plot(rf_out[1], G3_out[1], marker='8', markersize='10', markerfacecolor= 'None', color='darkgrey')
h1 = ax1.plot(rf_out[2:], G2_out[2:], marker='s', markersize = '10')
h2 = ax1.plot(rf_out[2:], G3_out[2:], marker='8', markersize = '10')

ax1.set_xlabel(r'$r_f$')
ax1.set_ylabel(r'Gauss coefficient ratios')
handles = [h1[-1], h2[-1]]
labels = [r'$G2$', r'$G3$']
ax1.legend(handles, labels)

if saveOn == 1:
    if not os.path.exists('{}'.format(saveDir)):
        os.makedirs('{}'.format(saveDir))
    fig1.savefig('{}/compare_rf.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig1.savefig('{}/compare_rf.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/compare_rf.*'.format(saveDir))
