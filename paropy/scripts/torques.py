'''
torques.py

Analyse gravitational, electromagnetic and viscous torques
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from paropy.coreproperties import icb_radius, cmb_radius
from paropy.data_utils import load_dimensionless, load_mantle, load_innercore

# matplotlib.use('Agg')  # backend for no display
plt.close('all')
#%% INPUT PARAMETERS
run_ID = ['chem_200d', 'ref_c', 'd_0_55a', 'd_0_6a', 'd_0_65b',
          'c-200a', 'd_0_75a', 'd_0_8a']  # PARODY simulation tag
dirName = '/data/geodynamo/wongj/Work'  # path containing simulation output

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/torques'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/torques'

plt.close('all')
# %% LOAD
w, h = plt.figaspect(fig_aspect)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
fig2, ax2 = plt.subplots(1, 1, figsize=(1.5*w, h))
shell_gap = cmb_radius - icb_radius
ri = icb_radius/shell_gap
ax1.axhline(0, linewidth=1, color='k', linestyle='--')

rf = np.zeros(len(run_ID))
Cf = np.zeros(len(run_ID))
Cicb = np.zeros(len(run_ID))
S = np.zeros(len(run_ID))
D = np.zeros(len(run_ID))
i = 0
for run in run_ID:
    fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))
    directory = '{}/{}'.format(dirName, run)
    print('Loading {}'.format(directory))
    (NR, Ek, Ra, Pr, Pm_out, fi, rf_out) = load_dimensionless(run, directory)
    df_mantle = load_mantle(run, directory)
    df_ic = load_innercore(run, directory)

    time = df_mantle.time.to_numpy()
    # shift to start from zero and use magnetic diffusion time
    time = (time - time[0])/Pm_out
    # magnetic_torque_mantle = df_mantle["magnetic_torque_on_mantle"].to_numpy() # NOTE: =0 due to insulating mantle
    gravitational_torque_mantle = df_mantle["gravitational_torque_on_mantle"].to_numpy()
    magnetic_torque_mantle = df_mantle["magnetic_torque_on_mantle"].to_numpy()
    # viscous_torque_ic = df_ic["viscous_torque_ic"].to_numpy() # NOTE: =0 due to stress free boundaries
    magnetic_torque_ic = df_ic["magnetic_torque_ic"].to_numpy()
    gravitational_torque_ic = df_ic["gravity_torque_ic"].to_numpy()
    total_angular_momentum = df_ic["total_angular_momentum_ic+oc+m"]
    print('Time average of total angular momentum = {:.3e}'.format(total_angular_momentum.mean()))

    # Plot
    ax1.plot(time, gravitational_torque_mantle, label = r'$\Gamma_G$')
    ax1.plot(time, magnetic_torque_ic, ls = ':', label = r'$\Gamma_{ICB}$') # NOTE: similar to gravitational torque on mantle
    # ax1.plot(time, gravitational_torque_ic) # NOTE: equal and opposite to gravitational torque on mantle
    ax1.set_xlabel(r'magnetic diffusion time')
    ax1.set_ylabel(r'torque')
    ax1.legend()
    ax1.autoscale(enable=True, axis='x', tight=True)

    if saveOn == 1:
        if not os.path.exists('{}'.format(saveDir)):
            os.makedirs('{}'.format(saveDir))
        fig1.savefig('{}/{}.png'.format(saveDir, run),
                     format='png', dpi=200, bbox_inches='tight')
        fig1.savefig('{}/{}.pdf'.format(saveDir, run),
                     format='pdf', dpi=200, bbox_inches='tight')
        print('Figures saved as {}/{}.*'.format(saveDir, run))

    plt.close(fig1)

    i += 1


    
