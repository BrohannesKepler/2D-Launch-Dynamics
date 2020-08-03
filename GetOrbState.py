#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:47:31 2020

@author: haider
"""

import numpy as np


def main(rx,ry,xi,eta):
    
    """
    Compute the final orbital state of the vehicle given terminal state vector
    
    """
    mu = 3.986e14
    # Specific orbital energy and momentum:
    r = (rx**2+ry**2)**0.5
    V = (xi**2 + eta**2)**0.5
    
    
    
    E = V**2/2 - mu/r
    h = r*xi
    
    # Semi major axis & Eccentricity:
    
    a = -mu/(2*E)
    e = (1 + 2*E*(h/mu)**2)**0.5
    
    #True anomaly:
    
    nu = np.arccos(1/e * (h**2/(mu*r) - 1))
    
    
    return np.array([a,e,nu])