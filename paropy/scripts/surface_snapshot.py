#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:13:57 2021

@author: wongj

Python version of PARODY-JA4.3 Matlab file 'Matlab/surfaceload_v4.m'.

Loads core surface data and plots core surface field snapshot.

ATTENTION: Folder structure should be of the form `<folder>/<run_ID>/Gt_*.run_ID`
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from paropy.data_utils import surfaceload
from paropy.plot_utils import rad_to_deg, get_Z_lim

#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = 'c-200a' # PARODY simulation tag
directory = '/Volumes/NAS/ipgp/Work/{}/'.format(run_ID) # path containing simulation output
timestamp = '16.84707134' 

fig_aspect = 1 # figure aspect ratio
n_levels = 61 # no. of contour levels
saveOn = 0 # save figures?

#%%----------------------------------------------------------------------------
# Load data
St_file = 'St={}.{}'.format(timestamp,run_ID)
filename = directory + St_file

(version, time, DeltaU, Coriolis, Lorentz, Buoyancy, ForcingU, 
            DeltaT, ForcingT, DeltaB, ForcingB, Ek, Ra, Pm, Pr,
            nr, ntheta, nphi, azsym, radius, theta, phi, Vt, Vp,
            Br, dtBr) = surfaceload(filename)

#%%----------------------------------------------------------------------------
# Plot
w, h = plt.figaspect(fig_aspect)
fig, ax = plt.subplots(1, 1, figsize=(1.5*w,h), 
                       subplot_kw={'projection': ccrs.Mollweide()})
X,Y = rad_to_deg(phi, theta)
Z = Br.T
Z_lim = get_Z_lim(Z)
levels = np.linspace(-Z_lim,Z_lim,n_levels)
c = ax.contourf(X, Y, Z, levels, transform=ccrs.PlateCarree(), cmap='PuOr_r',
                extend='both')
cbar_ax = fig.add_axes([0.2,0.06,0.6,0.04])
cbar = fig.colorbar(c, cax=cbar_ax, orientation='horizontal')
cbar.set_ticks([-Z_lim,-Z_lim/2,0,Z_lim/2,Z_lim])
cbar.ax.set_xlabel(r'$B_{r}$',fontsize=12)
cbar.ax.tick_params(labelsize=12)
cbar.ax.tick_params(length=6)
ax.gridlines()
ax.set_global()

# Save
if saveOn==1:
    if not os.path.exists(directory+'/figures'):
        os.makedirs(directory+'/figures')
    fig.savefig(directory+'/figures/surf_{}.png'.format(timestamp),format='png',
                dpi=200,bbox_inches='tight')
    fig.savefig(directory+'/figures/surf_{}.pdf'.format(timestamp),format='pdf',
                dpi=200,bbox_inches='tight')
