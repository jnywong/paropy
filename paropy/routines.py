#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:23:22 2021

@author: wongj
"""
import numpy as np
import shtns
import h5py

from scipy.integrate import trapz
from paropy.coreproperties import icb_radius, cmb_radius
from paropy.data_utils import parodyload, list_Gt_files, list_St_files, surfaceload

def sim_time(data):
    '''Total simulation time (viscous) from diagnostic data'''
    time = data.time.iloc[-1]-data.time.iloc[0]
    
    return time

def grav_torque(mantle_data):
    '''Mean and maximum gravitational torque'''
    gamma=np.mean(mantle_data.gravitational_torque_on_mantle)
    gamma_max = max(np.abs(mantle_data.gravitational_torque_on_mantle))
    
    return gamma,gamma_max

def meridional_timeavg(run_ID, directory):
    '''Time average of azimuthally averaged Gt files'''
    Gt_file = list_Gt_files(run_ID,directory) # Find all Gt_no in folder
    n = len(Gt_file)

    # Loop over time
    Vt_m, Vp_m, Vr_m, Br_m, Bt_m, Bp_m, T_m = [
        [] for _ in range(7)]
    i=0
    for file in Gt_file:
        print('Loading {} ({}/{})'.format(file, i+1, n))
        filename = '{}/{}'.format(directory,file)
        (_, _, _, _, _, _, _,
            _, _, _, _, _, _, _, _,
            _, _, _, _, radius, theta, phi, Vr, Vt, Vp,
            Br, Bt, Bp, T) = parodyload(filename)
        # Zonal averages
        Vr_m.append(np.mean(Vr,axis=0))
        # print('1')
        Vt_m.append(np.mean(Vt, axis=0))
        # print('2')
        Vp_m.append(np.mean(Vp,axis=0))
        # print('3')
        Br_m.append(np.mean(Br,axis=0))
        # print('4')
        Bt_m.append(np.mean(Bt, axis=0))
        # print('5')
        Bp_m.append(np.mean(Bp,axis=0))
        # print('6')
        T_m.append(np.mean(T,axis=0))
        # print('7')
        i += 1
    # Time average (should really divide by dt but n is good enough according to JA)
    Vr_out = sum(Vr_m)/n
    Vt_out = sum(Vt_m)/n
    Vp_out = sum(Vp_m)/n
    Br_out = sum(Br_m)/n
    Bt_out = sum(Bt_m)/n
    Bp_out = sum(Bp_m)/n
    T_out = sum(T_m)/n
    
    # Save
    with h5py.File('{}/meridional_timeavg'.format(directory), 'a') as f:
        f.create_dataset('radius', data=radius)
        f.create_dataset('theta', data=theta)
        f.create_dataset('phi', data=phi)
        f.create_dataset('Vr', data=Vr_out)
        f.create_dataset('Vt', data=Vt_out)
        f.create_dataset('Vp', data=Vp_out)
        f.create_dataset('Br', data=Br_out)
        f.create_dataset('Bt', data=Bt_out)
        f.create_dataset('Bp', data=Bp_out)
        f.create_dataset('T', data=T_out)

    return (radius, theta, phi, Vr_out, Vt_out, Vp_out, Br_out, Bt_out, Bp_out, T_out)

def surface_timeavg(run_ID, directory):
    '''
    Time average surface fields
    '''
    St_file = list_St_files(run_ID,directory) # Find all Gt_no in folder
    n = len(St_file)

    # Loop over time
    Vt_s, Vp_s, Br_s, dtBr_s = [[] for _ in range(4)]
    i=0
    for file in St_file:
        print('Loading {} ({}/{})'.format(file, i+1, n))
        filename = '{}/{}'.format(directory,file)

        (_, _, _, _, _, _, _,
          _, _, _, _, _, _, _, _,
            _, _, _, _, _, theta, phi, Vt, Vp, Br,
            dtBr) = surfaceload(filename)
        # Append
        Vt_s.append(Vt)
        Vp_s.append(Vp)
        Br_s.append(Br)
        dtBr_s.append(dtBr)
        i += 1
    # Time average (should really divide by dt but n is good enough according to JA)
    Vt_out = sum(Vt_s)/n
    Vp_out = sum(Vp_s)/n
    Br_out = sum(Br_s)/n
    dtBr_out = sum(dtBr_s)/n

    # Save 
    with h5py.File('{}/surface_timeavg'.format(directory), 'a') as f:
        f.create_dataset('theta', data=theta)
        f.create_dataset('phi', data=phi)
        f.create_dataset('Vt', data=Vt_out)
        f.create_dataset('Vp', data=Vp_out)
        f.create_dataset('Br', data=Br_out)
        f.create_dataset('dtBr', data=dtBr_out)
    print('{}/surface_timeavg saved'.format(directory))

    return (theta, phi, Vt_out, Vp_out, Br_out, dtBr_out)


def surface_phiavg_timeavg(run_ID, directory):
    '''
    Time average phi averaged surface fields
    '''
    St_file = list_St_files(run_ID, directory)  # Find all St_no in folder
    n = len(St_file)

    # Loop over time
    Vt_s, Vp_s, Br_s, dtBr_s = [[] for _ in range(4)]
    i = 0
    for file in St_file:
        print('Loading {} ({}/{})'.format(file, i+1, n))
        filename = '{}/{}'.format(directory, file)

        (_, _, _, _, _, _, _,
         _, _, _, _, _, _, _, _,
            _, _, _, _, _, theta, phi, Vt, Vp, Br,
            dtBr) = surfaceload(filename)
        # Append
        Vt_s.append(np.mean(Vt,axis=0))
        Vp_s.append(np.mean(Vp,axis=0))
        Br_s.append(np.mean(Br,axis=0))
        dtBr_s.append(np.mean(dtBr,axis=0))
        i += 1
    # Time average (should really divide by dt but n is good enough according to JA)
    Vt_out = sum(Vt_s)/n
    Vp_out = sum(Vp_s)/n
    Br_out = sum(Br_s)/n
    dtBr_out = sum(dtBr_s)/n

    # Save
    with h5py.File('{}/surface_phiavg_timeavg'.format(directory), 'w') as f:
        f.create_dataset('theta', data=theta)
        f.create_dataset('phi', data=phi)
        f.create_dataset('Vt', data=Vt_out)
        f.create_dataset('Vp', data=Vp_out)
        f.create_dataset('Br', data=Br_out)
        f.create_dataset('dtBr', data=dtBr_out)
    print('{}/surface_phiavg_timeavg saved'.format(directory))

    return (theta, phi, Vt_out, Vp_out, Br_out, dtBr_out)

def convective_power_timeavg(run_ID, directory):
    Gt_file = list_Gt_files(run_ID, directory)  # Find all Gt_no in folder
    n = len(Gt_file)

    # Loop over time
    Vr_avg, T_avg = [[] for _ in range(2)]
    i = 0
    for file in Gt_file:
        print('Loading {} ({}/{})'.format(file, i+1, n))
        filename = '{}/{}'.format(directory, file)
        (_, _, _, _, _, _, _,
            _, _, _, _, _, _, _, _,
            _, _, _, _, radius, theta, phi, Vr, _, _,
            _, _, _, T) = parodyload(filename)
        # Integrate over phi and theta
        Vr_avg.append(trapz(trapz(Vr,phi,axis=0)*np.sin(theta),theta,axis=0))
        T_avg.append(trapz(trapz(T, phi, axis=0)*np.sin(theta), theta, axis=0))
        i += 1
    # Time average (should really divide by dt but n is good enough according to JA)
    Vr_out = sum(Vr_avg)/n
    T_out = sum(T_avg)/n

    # Save
    with h5py.File('{}/convective_power'.format(directory), 'a') as f:
        f.create_dataset('radius', data=radius)
        f.create_dataset('Vr', data=Vr_out)
        f.create_dataset('T', data=T_out)

    return radius, Vr_out, T_out

def ref_codensity(r,rf,fi):
    '''
    Compute reference codensity
    '''
    nr = r.size
    shell_gap = cmb_radius-icb_radius
    ro = cmb_radius/shell_gap
    ri = icb_radius/shell_gap
    idx_f = np.where(r<rf)
    idx_o = np.where(r>=rf)

    if fi<0:
        So = -3/(ro**3-rf**3)
        Sf = 3*(1-fi)/(rf**3-ri**3)
        dCodr = So*(ro**3 - r[idx_o]**3)/(3*r[idx_o]**2)
        dCfdr = -Sf*(r[idx_f]**3-ri**3)/(3*r[idx_f]**2) - fi/r[idx_f]**2

        alpha = -(ro**3/rf+rf**2/2)/(ro**3 - rf**3) - \
            (1-fi)*(rf**2/2 + ri**3/rf)/(rf**3-ri**3) + fi/rf
        Co = (ro**3/r[idx_o] + r[idx_o]**2/2)/(ro**3 - rf**3) + alpha
        Cf = - (1-fi)*(r[idx_f]**2/2 + ri**3/r[idx_f])/(rf**3 - ri**3) + fi/r[idx_f]
        C = np.concatenate([Cf, Co])

    return C

def filter_field(Br,nphi,ntheta,l_trunc):
    m_max = l_trunc
    sh = shtns.sht(l_trunc, m_max)
    # default 'flag = sht_quick_init' uses gaussian grid
    nlat, nlon = sh.set_grid(nphi=nphi, nlat=ntheta)
    # NOTE: array has to be dtype='float64' and not 'float32'
    vr = Br.T.astype('float64')
    coeff = sh.analys(vr)  # spatial to spectral
    Br_f = sh.synth(coeff)  # spectral to spatial
    return Br_f
