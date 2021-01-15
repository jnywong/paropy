#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 09:29:58 2021

@author: wongj

Loads diagnostic outputs from PARODY-JA4.3:
    - kinetic and magnetic energies
    - Nusselt number (not sure this is computed in PARODY for all HeatingModes)
    - dipole e.g. g10, g11, B_rms
    - power e.g. convective power per unit volume
    - scales e.g. mean l for V, B, T (Christensen & Aubert, 2006)
    - spec_l and spec_m 
    - inner core and mantle rotation and torques (if Coupled Earth run)

"""
import os
import pandas as pd

from paropy.data_utils import load_kinetic,load_magnetic,load_nusselt, \
load_dipole,load_power,load_scales,load_spec_l,load_spec_m, load_mantle, \
load_innercore
from paropy.routines import sim_time, grav_torque

#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
run_ID = 'c-200a' # PARODY simulation tag
directory = '/Volumes/NAS/ipgp/Work/{}/'.format(run_ID) # path containing simulation output

#%%----------------------------------------------------------------------------
# Load data
kinetic_data=load_kinetic(run_ID,directory)
magnetic_data=load_magnetic(run_ID,directory)
# nusselt_data=load_nusselt(run_ID,directory)
# dipole_data=load_dipole(run_ID,directory)
# power_data=load_power(run_ID,directory)
# scales_data=load_scales(run_ID,directory)
# spec_l_data=load_spec_l(run_ID,directory)
# spec_m_data=load_spec_m(run_ID,directory)
try:
    mantle_data=load_mantle(run_ID,directory)
except FileNotFoundError:
    mantle_data=pd.DataFrame({'A' : []})
# try:
#     ic_data=load_innercore(run_ID,directory)         
# except:         
#     ic_data=pd.DataFrame({'A' : []})   

#%%----------------------------------------------------------------------------
# Output data
print('run_ID: {}'.format(run_ID))
time = sim_time(kinetic_data)
print('Simulation time: {:.3f}'.format(time))
(gamma,gamma_max)=grav_torque(mantle_data)
print('Mean of gravitational torque on mantle: {:.2f} ({:.4e} of the maximum absolute value)'.format(
                gamma, gamma/gamma_max))

# Plot
ax1=kinetic_data.plot("time","ke_per_unit_vol")
kinetic_data.plot("time","poloidal_ke",ax=ax1)
kinetic_data.plot("time","toroidal_ke",ax=ax1)

ax2=magnetic_data.plot("time","me_per_unit_vol")
magnetic_data.plot("time","poloidal_me",ax=ax2)
magnetic_data.plot("time","toroidal_me",ax=ax2)

ax3=mantle_data.plot("time","mantle_rotation_rate")

ax4=mantle_data.plot("time","gravitational_torque_on_mantle")