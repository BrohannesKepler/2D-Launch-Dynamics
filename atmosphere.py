#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 17:24:24 2020

@author: haider
"""
import numpy as np


def Atmosphere(h):
    
    """Function that computes density and temperature at a geopotential altitude
        h (m) above the ground using a simplified standard atmosphere model."""
        
    # CONSTANTS FOR EACH SEGMENT:
        
    re = 6378e3
    g0 = 9.81
    R = 287
    
    #First linear segment 0-11 km
    
    T1 = 288.16
    P1 = 101325
    rho1 = 1.225   
    a1 = -6.5E-3


    #First isothermal segment 11-25 km
    
    T2 = 216.66
    P2 = 2.2616E4
    rho2 = 0.3636
    
    # Second linear segment 25-47 km
    
    T3 = 216.66
    P3 = 2.484E3
    rho3 = 0.0399
    a3 = 3E-3
    
    #Second isothermal segment 47-53 km
    
    T4 = 282.66
    P4 = 120.0438
    rho4 = 0.0015
    
    #Third linear segment 53-79 km
    
    T5 = 282.66
    P5 = 58.1075
    rho5 = 7.2608E-4
    a5 = -4.5E-3
    
    #Convert geometric altitude to geopotential altitude:
    
    h = re*h/(re + h)
    
    if h <= 11E3:
        T = T1 + a1*h
        rho = rho1*(T/T1)**(-(g0/(a1*R) + 1))
    
    elif h > 11E3 and h <= 25E3:
        T = T2
        rho = rho2 * np.exp(-g0 * (h - 11E3)/(R*T2))
        
    elif h > 25E3 and h <= 47E3:
        
        T = T3 + a3*(h - 25E3)
        rho = rho3 * (T/T3)**(-(g0/(a3*R) + 1))
        
    elif h > 47E3 and h <= 53E3:
        
        T = T4
        rho = rho4 * np.exp(-g0 * (h-47E3)/(R*T4))
        
    elif h > 53E3 and h <= 79E3:
            
        T = T5 + a5*(h - 53E3)
        rho = rho5 * (T/T5)**(-(g0/(a5*R) + 1)); 
        
    elif h > 79E3:
        
        T = 165.66
        rho = 0
        
    return T, rho
        
        
        
        
    
    
    
    
    
    
    
    