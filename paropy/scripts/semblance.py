'''
semblance.py

Assess AD/NAD = axial vs non-axial dipole; O/E = degree of equatorial symmetry of the magnetic field; Z/NZ = relative power of axisymmetric components of non-dipole field

'''
#%% Import
import os
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from paropy.data_utils import load_dipole, load_magnetic, load_dimensionless, load_compliance
from paropy.routines import fcf

#%% Input parameters
run_ID = ['chem_200d', 'ref_c', 'd_0_55a', 'd_0_6a', 'd_0_65b',
          'c-200a', 'd_0_75a', 'd_0_8a']  # PARODY simulation tag
# run_ID = ['chem_200d']
# path containing simulation output
dirName = '/data/geodynamo/wongj/Work'

fig_aspect = 1  # figure aspect ratio

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/semblance'  # path to save files

#%% Pre-allocate
rf = np.zeros(len(run_ID))
ADNAD = np.zeros(len(run_ID))
OE = np.zeros(len(run_ID))
FCF = np.zeros(len(run_ID))
ZNZ = np.zeros(len(run_ID))
#%% Load
w, h = plt.figaspect(fig_aspect)
fig1, ((ax1,ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(2.5*w, 2*h), sharex=True)
fig2, ax2a = plt.subplots(1, 1, figsize=(1.5*w, h))
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
    if not os.path.exists('{}/fcf'.format(directory)):
        FCF_out = fcf(run, directory)  # timeavg
    else:  # load timeavg data
        print('Loading {}/fcf'.format(directory))
        f = h5py.File('{}/fcf'.format(directory), 'r')
        for key in f.keys():
            globals()[key] = np.array(f[key])
    FCF[i] = FCF_out

    i += 1

# Semblance - rating parameters and errors from Christensen (2010)
mu_ADNAD = 1.4
mu_OE = 1
mu_ZNZ = 0.15
mu_FCF = 1.5
sigma_ADNAD = 2
sigma_OE = 2
sigma_ZNZ = 2.5
sigma_FCF = 1.75
chi_ADNAD = ((np.log(ADNAD)-np.log(mu_ADNAD))/np.log(sigma_ADNAD))**2
chi_OE = ((np.log(OE)-np.log(mu_OE))/np.log(sigma_OE))**2
chi_ZNZ = ((np.log(ZNZ)-np.log(mu_ZNZ))/np.log(sigma_ZNZ))**2
chi_FCF = ((np.log(FCF)-np.log(mu_FCF))/np.log(sigma_FCF))**2
chi = chi_ADNAD + chi_OE + chi_ZNZ + chi_FCF
print('Semblance is {}'.format(chi))

err1_neg = mu_ADNAD - np.log(sigma_ADNAD)
err1_pos = mu_ADNAD + np.log(sigma_ADNAD)
err2_neg = mu_OE - np.log(sigma_OE)
err2_pos = mu_OE + np.log(sigma_OE)
err3_neg = mu_ZNZ - np.log(sigma_ZNZ)
err3_pos = mu_ZNZ + np.log(sigma_ZNZ)
err4_neg = mu_FCF - np.log(sigma_FCF)
err4_pos = mu_FCF + np.log(sigma_FCF)

ax1.axhline(1.4,linestyle='--',color='k')
ax2.axhline(1,linestyle='--',color='k')
ax3.axhline(0.15, linestyle='--', color='k')
ax4.axhline(1.5,linestyle='--',color='k')

ax1.plot(rf[0], ADNAD[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax1.plot(rf[1], ADNAD[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)
ax2.plot(rf[0], OE[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax2.plot(rf[1], OE[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)
ax3.plot(rf[0], ZNZ[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax3.plot(rf[1], ZNZ[1], 'o', color='None',
         markeredgecolor="darkgrey", markersize=10)
ax4.plot(rf[0], FCF[0], 'o', color='None', markeredgecolor="k", markersize=10)
ax4.plot(rf[1], FCF[1], 'o', color='None', markeredgecolor="darkgrey", markersize=10)

h1 = ax1.plot(rf[2:], ADNAD[2:], 'o', markersize=10,
              color='tab:blue', label=r'$AD/NAD$')
h2 = ax2.plot(rf[2:], OE[2:], 'o', markersize=10,
              color='tab:orange', label=r'$O/E$')
h3 = ax4.plot(rf[2:], ZNZ[2:], 'o', markersize=10,
              color='tab:red', label=r'$Z/NZ$')          
h4 = ax3.plot(rf[2:], FCF[2:], 'o', markersize=10,
              color='tab:green', label=r'$FCF$')

# handles = h1 + h2
# labels = [h.get_label() for h in handles]
# ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(r'$r_f$')
ax1.set_ylabel(r'$AD/NAD$')
ax2.set_xlabel(r'$r_f$')
ax2.set_ylabel(r'$O/E$')
ax3.set_xlabel(r'$r_f$')
ax3.set_ylabel(r'$Z/NZ$')
ax4.set_xlabel(r'$r_f$')
ax4.set_ylabel(r'$FCF$')
ax1.set_ylim([0,3])
ax2.set_ylim([0,2.5])
ax3.set_ylim([0,1.5])
ax4.set_ylim([0,2.8])

# Errors
ax1.margins(x=0)
ax2.margins(x=0)
ax3.margins(x=0)
ax4.margins(x=0)
xrange = ax1.get_xlim()
xr = np.linspace(xrange[0],xrange[1])
ax1.fill_between(xr,err1_neg,err1_pos, color = 'lightgray', alpha = 0.3, lw=0)
ax2.fill_between(xr,err2_neg,err2_pos, color = 'lightgray', alpha = 0.3, lw=0)
ax3.fill_between(xr,err3_neg,err3_pos, color = 'lightgray', alpha = 0.3, lw=0)
ax4.fill_between(xr,err4_neg,err4_pos, color = 'lightgray', alpha = 0.3, lw=0)

plt.tight_layout()

ax2a.plot(rf[0], chi[0],'o', color = 'k')
ax2a.plot(rf[1], chi[1],'o', color = 'darkgrey')
ax2a.plot(rf[2:], chi[2:],'o-')
ax2a.set_xlabel(r'$r_f$')
ax2a.set_ylabel(r'$\chi^2$')

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
    fig2.savefig('{}/chi.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig2.savefig('{}/chi.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/chi.*'.format(saveDir))
