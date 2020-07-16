#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 19:15:04 2020

@author: haider
"""

import numpy as np

def main(rx, ry, vx, vy, M, F9S2, t, t0):
    
    """
    Routine that executes the iterative guidance mode algorithm in order 
    to target the desired orbital radius and velocity 
    
    """
    
    # DECLARE INJECTION CONSTANTS:
    
    Re = 6378e3
    mu = 3.986e14
    
    RT = 200e3 + Re
    VT = np.sqrt(mu/RT)
    # In injection coordinates, for circular orbit we want 0 radial velocity
    xidott = VT
    etadott = 0
    etat = RT
    # CALCULATE CONSTANTS:
    
    R1 = (rx**2 + ry**2)**0.5
    V1 = (vx**2 + vy**2)**0.5
    
    phi1 =  np.arctan(rx/(Re+ry))
    
    g1 = -mu/R1**2
    gT = -mu/RT**2
    
    gstar = 0.5 * (g1 + gT)
    
    # VEHICLE ACCELERATION VARIABLES:
    
    tau = M/F9S2.mdot
    vex = F9S2.Isp * 9.81
    
    ##==================== BEGIN MAIN ROUTINE ================================
    
    # COMPUTE ESTIMATE FOR TIME-TO-GO:
    
    T2s = tau * (1 - np.exp((V1-VT)/vex))
    
    # STAGE RANGE ANGLE:
    
    phi12 = 1/ RT * (V1*T2s + vex * (T2s - (tau-T2s) * np.log(tau/(tau-T2s))))
    
    # TOTAL & AVERAGE RANGE ANGLE:
    
    phiT = phi1 + phi12
    phis = (phiT - phi1)/2
    
    # DEFINE ROTATION MATRIX 
    
    rot = np.array([[np.cos(phiT), -np.sin(phiT)], [np.sin(phiT), np.cos(phiT)]])
    pvec = np.array([rx, ry])
    vvec = np.array([vx, vy])
    
    # ROTATE TO GUIDANCE COORDINATES     
    pstateg = np.matmul(rot, pvec)
    vstateg = np.matmul(rot, vvec)
    
    # VELOCITY DEFICIENCY:
    
    dxi = xidott - vstateg[0] - (gstar * T2s * np.sin(phis))
    deta = etadott - vstateg[1] + (gstar * T2s * np.cos(phis))
    dVs = (dxi**2 + deta**2)**0.5

    # CONSTANTS FOR TGO UPDATE:

    lm = gstar * (deta * np.cos(phis) - (dxi * np.sin(phis)))
    L = vex * np.log(tau/(tau - T2s)) 
    k = vex/(tau - T2s)
    a = k**2 - gstar**2
    b = lm - L*k
    c = dVs**2 - L**2
    
    # TGO UPDATE
    
    dT2 = 1/a * (b + np.sqrt(b**2 + a*c))
    
    T2 = T2s + dT2
    
    # COMPUTE STEERING ANGLE IN INJECTION COORDINATES
    
    X_xi = np.arctan2((deta + gstar * dT2 * np.cos(phis)), (dxi - (gstar * dT2 * np.sin(phis))))
    
    # CONSTANTS FOR COMPUTATION IN INERTIAL COORDINATES
    
    A1 = vex * np.log(tau/(tau - T2))
    B1 = vex * (tau * np.log(tau/(tau - T2)) - T2)
    
    A2 = np.cos(X_xi) * vex * (T2 - (tau - T2)*np.log(tau/(tau - T2)))
    B2 = np.cos(X_xi) * -vex * (T2**2/2 - (tau*(T2 - (tau-T2) * np.log(tau/(tau - T2)))))
    
    C2 = pstateg[1] - etat + vstateg[1]*T2 - (0.5*gstar*T2**2 * np.cos(phis)) + \
        np.sin(X_xi) * (vex * (T2 - ((tau - T2) * np.log(tau/(tau - T2)))))
        
    K1 = B1 * C2/(A2*B1 - (A1*B2))
    K2 = A1*K1/B1
    
    # STEERING ANGLE IN THE INERTIAL FRAME 
    
    X = X_xi - K1 - phiT + K2*(t - t0)
    
    return T2, X






















    