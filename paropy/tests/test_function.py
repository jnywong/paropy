#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:51:42 2021

@author: wongj

Simple test code for parody_py
"""

import numpy as np
from paropy.plot_utils import rad_to_deg

def output(x,y):
    lon, lat = rad_to_deg(x, y)
    out = lon[20]
    return out

def test_output():
    # Generate lat/lon grid in radians
    x = np.linspace(0,2*np.pi)
    y = np.linspace(0,np.pi)
    assert output(x,y)==-33.06122448979593
