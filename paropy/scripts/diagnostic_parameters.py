'''
diagnostic_parameters.py

Calculate diagnostic parameters such as magnetic Reynolds and Elsasser numbers

  Re = np.sqrt(2*kinetic_data.ke_per_unit_vol.mean(axis=0))
  El = 2*magnetic_data.me_per_unit_vol.mean(axis=0)*Ek*Pm
  Rm = Pm*Re
  Ro = Re*Ek

'''
#%% Import
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from paropy.data_utils import load_kinetic, load_magnetic, load_dimensionless

#%% Input parameters
run_ID = ['d_0_55a','d_0_6a','d_0_65b','c-200a','d_0_75a','d_0_8a']  # PARODY simulation tag
# path containing simulation output
dirName = '/data/geodynamo/wongj/Work'
# directory = '/Volumes/NAS/ipgp/Work//'
# directory = '/Users/wongj/Desktop/data/'

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/diagnostic_parameters'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/surface'

#%% Pre-allocate
rf = np.zeros(len(run_ID))
time_mag = np.zeros(len(run_ID))
Re = np.zeros(len(run_ID))
El = np.zeros(len(run_ID))
Rm = np.zeros(len(run_ID))
Ro = np.zeros(len(run_ID))

#%% Load
i=0
for run in run_ID:
    directory = '{}/{}'.format(dirName,run)
    (_, Ek_out, Ra_out, Pr_out, Pm_out, fi_out, rf_out) = load_dimensionless(run,directory)
    kinetic_data = load_kinetic(run,directory)
    magnetic_data = load_magnetic(run, directory)

    #%% Diagnostic parameters
    time_out = (kinetic_data.time.iloc[-1]-kinetic_data.time.iloc[0])/Pm_out
    Re_out = np.sqrt(2*kinetic_data.ke_per_unit_vol.mean(axis=0))
    El_out = 2*magnetic_data.me_per_unit_vol.mean(axis=0)*Ek_out*Pm_out
    Rm_out = Pm_out*Re_out
    Ro_out = Re_out*Ek_out

    rf[i]=rf_out
    time_mag[i] = time_out
    Re[i]=Re_out
    El[i]=El_out
    Rm[i]=Rm_out
    Ro[i]=Ro_out

    i+=1

#%% Plot
w, h = plt.figaspect(fig_aspect)
fig, ax1 = plt.subplots(1, 1, figsize=(1.5*w,h))
ax2 = ax1.twinx()

ax1.plot(rf,Rm,'o',color='tab:blue')
ax2.plot(rf,El,'o',color='tab:orange')

ax1.set_xlabel(r'$r_f$')
ax1.set_ylabel(r'$Rm$')
ax2.set_ylabel(r'$\Lambda$')

# Save
if saveOn == 1:
    if not os.path.exists(saveDir+'/diagnostic_parameters'):
        os.makedirs(saveDir+'/diagnostic_parameters')
    fig.savefig(saveDir+'/diagnostic_parameters.png', format='png',
                dpi=200, bbox_inches='tight')
    fig.savefig(saveDir+'/diagnostic_parameters.pdf', format='pdf',
                dpi=200, bbox_inches='tight')
    print('Figures saved as {}/diagnostic_parameters.*'.format(saveDir))
