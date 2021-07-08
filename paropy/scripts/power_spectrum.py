'''
power_spectrum.py

Compare power spectra of different runs
'''

import os
import matplotlib
import matplotlib.pyplot as plt

from paropy.data_utils import load_dipole, load_spec_l

# matplotlib.use('Agg')  # backend for no display
plt.close('all')
#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = ['chem_200d','d_0_55a','d_0_6a','d_0_65b','c-200a','d_0_75a','d_0_8a']  # PARODY simulation tag
l_trunc = 13  # SH degree truncation

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
    df_dipole = load_dipole(run,directory)
    df_spec_l = load_spec_l(run,directory)
    df_spec_l.iloc[1:l_trunc+1].plot("sh_degree", "timeavg_field", marker='o', ax=ax1)

ax1.legend(run_ID)

# Save 
if not os.path.exists(saveDir):
    os.makedirs(saveDir)
fig1.savefig(saveDir+'/compare_timeavg_spec_l.png')
print('Figures saved in {}/compare_timeavg_spec_l.png'.format(saveDir))
