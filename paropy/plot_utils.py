#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:10:19 2021

@author: wongj
"""

import numpy as np
import math

def rad_to_deg(phi,theta):
    '''Converts radians into longitudinal and latitudinal degrees'''
    i=0
    phi_deg=np.zeros(len(phi))
    
    for val in phi:
        phi_deg[i]=math.degrees(val)-180
        i=i+1
    theta_deg=np.zeros(len(theta))
    i=0
    for val in theta:
        theta_deg[i]=math.degrees(val)-90
        i=i+1    
    # phi, theta = np.meshgrid(phi_degree, theta_degree) # indexing='xy')
    
    return (phi_deg, theta_deg)

def get_Z_lim(Z):
    '''Choose Z limit for plot to the nearest sig. fig. modulo 5'''
    Z_lim = np.max([np.abs(Z.min()),np.abs(Z.max())])
    index = np.ceil(-np.log10(Z_lim))
    modulo = Z_lim % (5*10**(-index))
    if modulo!=Z_lim:
        Z_lim = Z_lim - (Z_lim % (5*10**(-index)))
    else:
        Z_lim = np.round(Z_lim, index)

    return Z_lim        
    
    