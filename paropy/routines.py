#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:23:22 2021

@author: wongj
"""
import numpy as np
import h5py

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
            _, _, _, _, _, _, _, Vr, Vt, Vp,
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

    return (Vr_out, Vt_out, Vp_out, Br_out, Bt_out, Bp_out, T_out)

def surface_timeavg(run_ID, directory):
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
    print('hey tobes')

    # Save 
    # data = [theta, phi, Vt_out, Vp_out, Br_out, dtBr_out]
    # np.savez('{}/surface_timeavg.npz'.format(directory), *data)
    # print('{}/surface_timeavg.npz saved'.format(directory))
    with h5py.File('{}/surface_timeavg'.format(directory), 'w') as f:
        f.create_dataset('theta', data=theta)
        print('hey brah')
        f.create_dataset('phi', data=phi)
        print('hey poop')
        f.create_dataset('Vt', data=Vt_out)
        print('hey Gris')
        f.create_dataset('Vp', data=Vp_out)
        print('hey wow')
        f.create_dataset('Br', data=Br_out)
        print('hey poopbutt')
        f.create_dataset('dtBr', data=dtBr_out)
    print('{}/surface_timeavg saved'.format(directory))

    return (theta, phi, Vt_out, Vp_out, Br_out, dtBr_out)

