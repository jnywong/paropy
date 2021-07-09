'''
f_b_dip.py

Assess ratios f_dip = dipole rms/ surface rms in 1<degree<12 and b_dip = Elsasser/RMS field at core surface = B_deep/B_surf

'''
#%% Import
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from paropy.data_utils import load_dipole, load_magnetic, load_dimensionless

#%% Input parameters
run_ID = ['chem_200d', 'ref_c', 'd_0_55a','d_0_6a','d_0_65b','c-200a','d_0_75a','d_0_8a']  # PARODY simulation tag
# run_ID = ['chem_200d']
# path containing simulation output
dirName = '/data/geodynamo/wongj/Work'
# directory = '/Volumes/NAS/ipgp/Work//'
# directory = '/Users/wongj/Desktop/data/'

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/f_b_dip'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/f_b_dip'

#%% Pre-allocate
rf = np.zeros(len(run_ID))
Bsurf_rms = np.zeros(len(run_ID))
Bsurf_rms_12 = np.zeros(len(run_ID))
dipole_rms = np.zeros(len(run_ID))
El = np.zeros(len(run_ID))
#%% Load
w, h = plt.figaspect(fig_aspect)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
ax2 = ax1.twinx()

i = 0
for run in run_ID:
    directory = '{}/{}'.format(dirName, run)
    print('Loading {}'.format(directory))
    (_, Ek_out, Ra_out, Pr_out, Pm_out, fi_out,
     rf_out) = load_dimensionless(run, directory)
    magnetic_data = load_magnetic(run, directory)
    dipole_data = load_dipole(run, directory)
    rf[i] = rf_out
    Bsurf_rms[i] = dipole_data["surface_rms_B"].mean()
    Bsurf_rms_12[i] = dipole_data["surface_rms_B_deg12"].mean()
    dipole_rms[i] = dipole_data["dipole_rms_B"].mean()
    El[i] = 2*magnetic_data.me_per_unit_vol.mean(axis=0)*Ek_out*Pm_out

    i+=1

fdip = dipole_rms/Bsurf_rms_12
bdip = El/Bsurf_rms

ax1.plot(rf[0], fdip[0], 'o', color = 'k', markerfacecolor="None", markersize=10)
ax2.plot(rf[0], bdip[0], 's', color = 'k', markerfacecolor="None", markersize=10)
ax1.plot(rf[1], fdip[1], 'o', color = 'darkgrey', markerfacecolor="None", markersize=10)
ax2.plot(rf[1], bdip[1], 's', color = 'darkgrey', markerfacecolor="None", markersize=10)
h1 = ax1.plot(rf[2:], fdip[2:], 'o', markersize = 10, color='tab:blue', label=r'$f_{dip}$')
h2 = ax2.plot(rf[2:], bdip[2:], 's', markersize = 10, color = 'tab:orange', label = r'$b_{dip}$')

handles = h1 + h2
labels = [h.get_label() for h in handles]
ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(r'$r_f$')
ax1.set_ylabel(r'$f_{dip}$')
ax2.set_ylabel(r'$b_{dip}$')

# Save
if saveOn == 1:
    if not os.path.exists('{}'.format(saveDir)):
        os.makedirs('{}'.format(saveDir))
    fig1.savefig('{}/compare_rf.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig1.savefig('{}/compare_rf.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/compare_rf.*'.format(saveDir))
