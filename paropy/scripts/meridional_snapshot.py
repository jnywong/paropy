#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:49:39 2021

@author: wongj

Python version of PARODY-JA4.3 Matlab file 'Matlab/parodyloadload_v4.m'.

Loads graphics file and plots snapshots of the azimuthal velocity field,
azimuthal magnetic field, temperature/codensity field in meridional slices.

"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmocean.cm as cmo

from paropy.data_utils import parodyload
from paropy.plot_utils import flayer_outline, streamfunction, T_shift, merid_outline

matplotlib.use('Agg')  # backend for no display

#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = 'd_0_6a'
rf = 0.6
directory = '/data/geodynamo/wongj/Work/{}/'.format(run_ID)  # path containing
# timestamp = '20.28436204'
timestamp ='21.49797360'
# timestamp = '21.17830229'

fig_aspect = 1 # figure aspect ratio
n_levels = 21 # no. of contour levels
Vmax = 250 # max Vp
Bmax = 2.5 # max Bp
Tr_min = 1.23

saveOn = 1 # save figures?
saveDir = '/home/wongj/Work/figures/meridional/'  # path to save files

#%%----------------------------------------------------------------------------
# Load data
Gt_file = 'Gt={}.{}'.format(timestamp,run_ID)
filename = directory + Gt_file

(version, time, DeltaU, Coriolis, Lorentz, Buoyancy, ForcingU,
            DeltaT, ForcingT, DeltaB, ForcingB, Ek, Ra, Pm, Pr,
            nr, ntheta, nphi, azsym, radius, theta, phi, Vr, Vt, Vp,
            Br, Bt, Bp, T) = parodyload(filename)

#%%----------------------------------------------------------------------------
# Plot
w, h = plt.figaspect(fig_aspect)
fig = plt.figure(constrained_layout=True, figsize = (2*w,h))
spec = gridspec.GridSpec(ncols = 3, nrows = 1, figure=fig)

# Velocity
ax = fig.add_subplot(spec[0,0])
X = np.outer(np.sin(theta),radius)
Y = np.outer(np.cos(theta),radius)
# azimuthal
Z = np.mean(Vp,0)
# Z_lim = get_Z_lim(Z)
Z_lim = Vmax
levels = np.linspace(-Z_lim,Z_lim,n_levels)
c = ax.contourf(X,Y,Z,levels,cmap='RdYlBu_r',extend='both')
cbar=plt.colorbar(c,ax=ax, aspect = 50, ticks=levels[::2])
cbar.ax.set_title(r'$\mathbf{u}$')
# streamfunction
Vr_m = np.mean(Vr,0)
Vt_m = np.mean(Vt,0)
Z = streamfunction(radius, theta, Vr_m, Vt_m)
c = ax.contour(X,Y,Z, 9 , colors='grey', alpha = 0.5)
# x,y = generate_semicircle(0,0,radius[0], 1e-4)
# ax.plot(x, y, 'k')
# x,y = generate_semicircle(0,0,radius[-1], 1e-5)
# ax.plot(x, y, 'k')
# ax.vlines(0,radius[0],radius[-1],'k')
# ax.vlines(0,-radius[0],-radius[-1],'k')
merid_outline(ax,radius)
ax.axis('off')

# Field
ax = fig.add_subplot(spec[0,1])
# azimuthal
Z = np.mean(Bp,0)
Z_lim = Bmax
levels = np.linspace(-Z_lim,Z_lim,n_levels)
c = ax.contourf(X,Y,Z,levels,cmap='PuOr_r',extend='both')
cbar=plt.colorbar(c,ax=ax, aspect = 50, ticks=levels[::2])
cbar.ax.set_title(r'$\mathbf{B}$')
# poloidal
Br_m = np.mean(Br,0)
Bt_m = np.mean(Bt,0)
Z = streamfunction(radius, theta, Br_m, Bt_m)
c = ax.contour(X,Y,Z, 9 , colors='grey', alpha = 0.5)
merid_outline(ax,radius)
ax.axis('off')

ax = fig.add_subplot(spec[0,2])
# FIX: C shift
# Z = np.mean(T,0)
Z0 = np.mean(T,0)
idx = np.argwhere(radius>rf)[0][0]
# Ts_max = np.mean(Z0[:,idx+1]) # max T near top of F-layer
# Ts_min = np.mean(Z0[:,-1]) # min T near CMB 
Ts_max = np.max(Z0)
Ts_min = np.min(Z0)
h = Ts_max-Ts_min
Z1 = (Z0 - Ts_max)/h
lev_max = 0.5
# lev_max = np.round(max(abs(np.max(Z)),abs(np.min(Z))),1)
Z = Z1 + lev_max
# Z = Z1 + lev_max
lev_min = -lev_max
levels=np.linspace(lev_min,lev_max,n_levels)
c = ax.contourf(X,Y,Z,levels,cmap='inferno',extend='both')
# c = ax.contourf(X,Y,Z,cmap='inferno',extend='both')
cbar=plt.colorbar(c,ax=ax,aspect = 50, ticks=levels[::2])
cbar.ax.set_title(r'$C$')
merid_outline(ax,radius)
flayer_outline(ax,rf)
ax.axis('off')

# Save
if saveOn==1:
    if not os.path.exists(saveDir+'{}'.format(run_ID)):
        os.makedirs(saveDir+'{}'.format(run_ID))
    fig.savefig(saveDir+'{}/{}.png'.format(run_ID, timestamp),format='png',
                dpi=200,bbox_inches='tight')
    fig.savefig(saveDir+'{}/{}.pdf'.format(run_ID, timestamp),format='pdf',
                dpi=200,bbox_inches='tight')
    print('Figures saved as {}{}/{}.*'.format(saveDir,run_ID,timestamp))
