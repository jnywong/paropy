'''
filter_surface_field.py

Use SHTns library to filter core surface magnetic field
'''

import os
import numpy as np
import shtns
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from paropy.data_utils import surfaceload
from paropy.plot_utils import rad_to_deg, get_Z_lim

# matplotlib.use('Agg')  # backend for no display
#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = 'c-200a'  # PARODY simulation tag
# path containing simulation output
directory = '/data/geodynamo/wongj/Work/{}'.format(run_ID)
# directory = '/Volumes/NAS/ipgp/Work/{}/'.format(run_ID)
# directory = '/Users/wongj/Desktop/data/{}'.format(run_ID)

timestamp = '16.84707134'
l_max = 133 # max. spherical harmonic degree from simulation
l_trunc = 14 # SH degree truncation

fig_aspect = 1  # figure aspect ratio
n_levels = 61  # no. of contour levels

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/filter_surface_field'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/filter_surface_field'
#%% Load data
St_file = 'St={}.{}'.format(timestamp, run_ID)
filename = '{}/{}'.format(directory, St_file)

(version, time, DeltaU, Coriolis, Lorentz, Buoyancy, ForcingU,
 DeltaT, ForcingT, DeltaB, ForcingB, Ek, Ra, Pm, Pr,
 nr, ntheta, nphi, azsym, radius, theta, phi, Vt, Vp,
 Br, dtBr) = surfaceload(filename)

#%% SH transform
m_max = l_trunc
sh = shtns.sht(l_trunc, m_max)
nlat, nlon = sh.set_grid(nphi=nphi, nlat = ntheta)
vr = Br.T.astype('float64') # NOTE: array has to be dtype='float64' and not 'float32'

coeff = sh.analys(vr) # spatial to spectral
out = sh.synth(coeff) # spectral to spatial

#%% Plot
w, h = plt.figaspect(fig_aspect)
fig, ax = plt.subplots(1, 1, figsize=(1.5*w, h),
                       subplot_kw={'projection': ccrs.Mollweide()})
X, Y = rad_to_deg(phi, theta)
Z = out
Z_lim = get_Z_lim(Z)
levels = np.linspace(-Z_lim, Z_lim, n_levels)
c = ax.contourf(X, Y, Z, levels, transform=ccrs.PlateCarree(), cmap='PuOr_r',
                extend='both')
cbar_ax = fig.add_axes([0.2, 0.06, 0.6, 0.04])
cbar = fig.colorbar(c, cax=cbar_ax, orientation='horizontal')
cbar.set_ticks([-Z_lim, -Z_lim/2, 0, Z_lim/2, Z_lim])
cbar.ax.set_xlabel(r'$B_{r}$', fontsize=12)
cbar.ax.tick_params(labelsize=12)
cbar.ax.tick_params(length=6)
ax.gridlines()
ax.set_global()

if saveOn == 1:
    if not os.path.exists('{}/{}'.format(saveDir, run_ID)):
        os.makedirs('{}/{}'.format(saveDir, run_ID))
    fig.savefig('{}/{}/{}.png'.format(saveDir, run_ID, timestamp),
                format='png', dpi=200, bbox_inches='tight')
    fig.savefig('{}/{}/{}.pdf'.format(saveDir, run_ID, timestamp),
                format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/{}/{}.*'.format(saveDir, run_ID, timestamp))
