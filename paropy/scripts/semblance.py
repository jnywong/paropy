'''
semblance.py

Assess AD/NAD = axial vs non-axial dipole; O/E = degree of equatorial symmetry of the magnetic field; Z/NZ = relative power of axisymmetric components of non-dipole field

'''
#%% Import
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from paropy.data_utils import load_dipole, load_magnetic, load_dimensionless, load_compliance

#%% Input parameters
run_ID = ['chem_200d', 'ref_c', 'd_0_55a', 'd_0_6a', 'd_0_65b',
          'c-200a', 'd_0_75a', 'd_0_8a']  # PARODY simulation tag
# path containing simulation output
dirName = '/data/geodynamo/wongj/Work'

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/semblance'  # path to save files

#%% Pre-allocate
rf = np.zeros(len(run_ID))
ADNAD = np.zeros(len(run_ID))
OE = np.zeros(len(run_ID))
ZNZ = np.zeros(len(run_ID))
#%% Load
w, h = plt.figaspect(fig_aspect)
fig1, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(2.5*w, h))
# ax2 = ax1.twinx()

i = 0
for run in run_ID:
    directory = '{}/{}'.format(dirName, run)
    print('Loading {}'.format(directory))
    (_, Ek_out, Ra_out, Pr_out, Pm_out, fi_out,
     rf_out) = load_dimensionless(run, directory)
    compliance_data = load_compliance(run, directory)
    rf[i] = rf_out
    ADNAD[i] = compliance_data["ADNAD"].mean()
    OE[i] = compliance_data["OE"].mean()
    ZNZ[i] = compliance_data["ZNZ"].mean()

    i += 1

ax1.plot(rf[0], ADNAD[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax1.plot(rf[1], ADNAD[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)
ax2.plot(rf[0], OE[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax2.plot(rf[1], OE[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)
ax3.plot(rf[0], ZNZ[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax3.plot(rf[1], ZNZ[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)
h1 = ax1.plot(rf[2:], ADNAD[2:], 'o', markersize=10,
              color='tab:blue', label=r'$AD/NAD$')
h2 = ax2.plot(rf[2:], OE[2:], 'o', markersize=10,
              color='tab:orange', label=r'$O/E$')
h3 = ax3.plot(rf[2:], ZNZ[2:], 'o', markersize=10,
              color='tab:green', label=r'$Z/NZ$')              

# handles = h1 + h2
# labels = [h.get_label() for h in handles]
# ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(r'$r_f$')
ax1.set_ylabel(r'$AD/NAD$')
ax2.set_xlabel(r'$r_f$')
ax2.set_ylabel(r'$O/E$')
ax3.set_xlabel(r'$r_f$')
ax3.set_ylabel(r'$Z/NZ$')
plt.tight_layout()

# Save
if saveOn == 1:
    if not os.path.exists('{}'.format(saveDir)):
        os.makedirs('{}'.format(saveDir))
    fig1.savefig('{}/compare_rf.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig1.savefig('{}/compare_rf.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/compare_rf.*'.format(saveDir))
