'''
power_spectrum.py

Compare power spectra of different runs
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from paropy.data_utils import load_dipole, load_spec_l, load_magnetic

# matplotlib.use('Agg')  # backend for no display
plt.close('all')
#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = ['chem_200d','ref_c', 'd_0_55a','d_0_6a','d_0_65b','c-200a','d_0_75a','d_0_8a']  # PARODY simulation tag
# run_ID = ['chem_200d']
l_trunc = 14  # SH degree truncation

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/power_spectrum'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/power_spectrum'
#%% Load data
w, h = plt.figaspect(fig_aspect)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
i = 0
for run in run_ID:
    # path containing simulation output
    directory = '/data/geodynamo/wongj/Work/{}'.format(run)
    df_mag = load_magnetic(run,directory)
    # print('Mean magnetic energy is {:.2f}'.format(df_mag["me_per_unit_vol"].mean()))
    # df_dipole = load_dipole(run,directory)
    df_spec_l = load_spec_l(run,directory)
    # print('Sum over degree l of power spectrum is {:.2f}'.format(df_spec_l["timeavg_field"].sum()))
    spec_l = df_spec_l.to_numpy()
    if i==0:
        ax1.plot(spec_l[1:15,0], spec_l[1:15,4], marker='o', markersize=10, color='k')
    elif i==1:
        ax1.plot(spec_l[1:15,0], spec_l[1:15,4], marker='o', markersize=10, color='darkgrey')
    else:
        ax1.plot(spec_l[1:15,0], spec_l[1:15,4], marker='o', markersize=10)
    i+=1


ax1.set_xticks(np.arange(1,l_trunc+1,1))

ax1.legend(run_ID)
plt.tight_layout()

# Save 
if saveOn==1:
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    fig1.savefig(saveDir+'/compare_timeavg_spec_l.png')
    print('Figures saved in {}/compare_timeavg_spec_l.png'.format(saveDir))
