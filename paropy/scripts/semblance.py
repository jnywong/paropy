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
from matplotlib.cm import ScalarMappable
from paropy.data_utils import load_dipole, load_magnetic, load_dimensionless, load_compliance
from paropy.routines import fcf
cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(vmin=1.2, vmax=5.2)
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
chi_ADNAD = np.zeros(len(run_ID))
chi_OE = np.zeros(len(run_ID))
chi_FCF = np.zeros(len(run_ID))
chi_ZNZ = np.zeros(len(run_ID))
chi = np.zeros(len(run_ID))
#%% Load
w, h = plt.figaspect(fig_aspect)
fig1, ((ax1,ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(2.5*w, 2*h), sharex=True)

# Semblance - rating parameters and errors from Christensen (2010)
mu_ADNAD = 1.4
mu_OE = 1
mu_ZNZ = 0.15
mu_FCF = 1.5
# NOTE: error = mu*sigma and mu/sigma not standard deviation
sigma_ADNAD = 2
sigma_OE = 2
sigma_ZNZ = 2.5
sigma_FCF = 1.75

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
    chi_ADNAD[i] = ((np.log(ADNAD[i])-np.log(mu_ADNAD))/np.log(sigma_ADNAD))**2
    chi_OE[i] = ((np.log(OE[i])-np.log(mu_OE))/np.log(sigma_OE))**2
    chi_ZNZ[i] = ((np.log(ZNZ[i])-np.log(mu_ZNZ))/np.log(sigma_ZNZ))**2
    chi_FCF[i] = ((np.log(FCF[i])-np.log(mu_FCF))/np.log(sigma_FCF))**2
    chi[i] = chi_ADNAD[i] + chi_OE[i] + chi_ZNZ[i] + chi_FCF[i]
    print('Semblance is {}'.format(chi[i]))

    if i==0:
        ax1.plot(rf[0], ADNAD[0], 'o', color=cmap(norm(chi[i])),
                 markeredgecolor="k", markersize=10)
        ax2.plot(rf[0], OE[0], 'o', color=cmap(norm(chi[i])),
                 markeredgecolor="k", markersize=10)
        ax3.plot(rf[0], ZNZ[0], 'o', color=cmap(norm(chi[i])), markeredgecolor="k", markersize=10)
        ax4.plot(rf[0], FCF[0], 'o', color=cmap(norm(chi[i])),
                 markeredgecolor="k", markersize=10)
    elif i==1:
        ax1.plot(rf[1], ADNAD[1], 'o', color=cmap(norm(chi[i])),
                 markeredgecolor="darkgrey", markersize=10)
        ax2.plot(rf[1], OE[1], 'o', color=cmap(norm(chi[i])),
                markeredgecolor="darkgrey", markersize=10)
        ax3.plot(rf[1], ZNZ[1], 'o', color=cmap(norm(chi[i])),
                markeredgecolor="darkgrey", markersize=10)
        ax4.plot(rf[1], FCF[1], 'o', color=cmap(norm(chi[i])),
                markeredgecolor="darkgrey", markersize=10)
    else:
        h = ax1.plot(rf[i], ADNAD[i], 'o', color=cmap(norm(chi[i])), markersize=10)
        ax2.plot(rf[i], OE[i], 'o', color=cmap(norm(chi[i])), markersize=10)
        ax3.plot(rf[i], ZNZ[i], 'o', color=cmap(norm(chi[i])), markersize=10)
        ax4.plot(rf[i], FCF[i], 'o', color=cmap(norm(chi[i])), markersize=10)
    i += 1

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

# Mean and errors
ax1.axhline(mu_ADNAD, color='k', linestyle='--',lw = 2)
ax2.axhline(mu_OE, color='k', linestyle='--',lw = 2)
ax3.axhline(mu_ZNZ, color='k', linestyle='--',lw = 2)
ax4.axhline(mu_FCF, color='k', linestyle='--',lw = 2)
ax1.margins(x=0)
ax2.margins(x=0)
ax3.margins(x=0)
ax4.margins(x=0)
xrange = ax1.get_xlim()
xr = np.linspace(xrange[0],xrange[1])
ax1.fill_between(xr,mu_ADNAD/sigma_ADNAD,mu_ADNAD*sigma_ADNAD, color = 'lightgray', alpha = 0.3, lw=0)
ax2.fill_between(xr,mu_OE/sigma_OE,mu_OE*sigma_OE, color = 'lightgray', alpha = 0.3, lw=0)
ax3.fill_between(xr,mu_ZNZ/sigma_ZNZ,mu_ZNZ*sigma_ZNZ, color = 'lightgray', alpha = 0.3, lw=0)
ax4.fill_between(xr,mu_FCF/sigma_FCF,mu_FCF*sigma_FCF, color = 'lightgray', alpha = 0.3, lw=0)

sm = ScalarMappable(norm = norm, cmap = cmap)
sm.set_array([])
cbar_ax = fig1.add_axes([0.95, 0.35, 0.015, 0.3])
cbar = fig1.colorbar(sm, cax=cbar_ax, label=r"semblance $(\chi^2)$")
cbar.set_ticks(np.linspace(1.2,5.2,5))

# Save
if saveOn == 1:
    if not os.path.exists('{}'.format(saveDir)):
        os.makedirs('{}'.format(saveDir))
    fig1.savefig('{}/compare_rf.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig1.savefig('{}/compare_rf.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/compare_rf.*'.format(saveDir))
