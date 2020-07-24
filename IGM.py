#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 19:15:04 2020

@author: haider
"""

import numpy as np
firstrun = True
def main(rx, ry, vx, vy, M, F9S2, t, t0, T2, IGMState):
    global firstrun
    """
    Routine that executes the iterative guidance mode algorithm in order 
    to target the desired orbital radius and velocity 
    
    """
    
    # DECLARE INJECTION CONSTANTS:
    
    Re = 6378e3
    mu = 3.986e14
    
    RT = 300e3 + Re
    VT = np.sqrt(mu/RT)
    # In injection coordinates, for circular orbit we want 0 radial velocity
    xidott = VT
    etadott = 0
    etat = RT
    
    # CALCULATE CONSTANTS:
    
    R1 = (rx**2 + ry**2)**0.5
    V1 = (vx**2 + vy**2)**0.5
    
    phi1 =  np.arctan(rx/ry)
    
    g1 = mu/R1**2
    gT = mu/RT**2
    
    gstar = 0.5 * (g1 + gT)
    
    # VEHICLE ACCELERATION VARIABLES:
    
    tau = M/F9S2.mdot
    vex = F9S2.Isp * 9.81
    
    ##==================== BEGIN MAIN ROUTINE ================================
    
    
    # CHECK TIME TO GO AND BYPASS ROUTINE IF SMALL TO AVOID INDETERMINATE SOLUTION
    
    if T2 < 10:
        dbg = IGMState
        T2 = T2 - 0.1
        return T2, dbg[-1], dbg 
    
    
    # COMPUTE T2*:
    
    T2s = tau * (1 - np.exp((V1 - VT)/vex))
    
    # IF FIRST GUIDANCE CALL ESTIMATE FOR T2 = T2STAR
    # ELSE, DECREMENT T2' BY THE CYCLE TIME 
    
    if t < t0 + 0.1:
        T2p = T2s
    else:
        T2p = T2 - 0.1
    
    
    # STAGE RANGE ANGLE:
    
    phi12 = 1/ RT * (V1*T2s + vex * (T2s - ((tau-T2s) * np.log(tau/(tau-T2s)))))
    
    # TOTAL & AVERAGE RANGE ANGLE:
    
    phiT = phi1 + phi12
    phis = phi12/2 
  
    # DEFINE ROTATION MATRIX 
    
    rot = np.array([[np.cos(phiT), -np.sin(phiT)], [np.sin(phiT), np.cos(phiT)]])
    pvec = np.array([rx, ry])
    vvec = np.array([vx, vy])
    
    # ROTATE TO GUIDANCE COORDINATES     
    pstateg = np.matmul(rot, pvec)
    vstateg = np.matmul(rot, vvec)
    
    # VELOCITY DEFICIENCY:
    
    dxi = xidott - vstateg[0] - (gstar * T2p * np.sin(phis))
    deta = etadott - vstateg[1] + (gstar * T2p * np.cos(phis))
    dVs = (dxi**2 + deta**2)**0.5

    # CONSTANTS FOR TGO UPDATE:

    lm = gstar * (deta * np.cos(phis) - (dxi * np.sin(phis)))
    L = vex * np.log(tau/(tau - T2p)) 
    k = vex/(tau - T2p)
    a = k**2 - gstar**2
    b = lm - L*k
    c = dVs**2 - L**2
    
    # TGO UPDATE
    
    dT2 = 1/a * (b + np.sqrt(b**2 + a*c))
    
    T2 = T2p + dT2
    
    # COMPUTE STEERING ANGLE IN INJECTION COORDINATES
    
    X_xi = np.arctan2((deta + gstar * dT2 * np.cos(phis)), (dxi - (gstar * dT2 * np.sin(phis))))
    
    # CONSTANTS FOR COMPUTATION IN INERTIAL COORDINATES
    
    A1 = vex * np.log(tau/(tau - T2))
    B1 = vex * (tau * np.log(tau/(tau - T2)) - T2)
    
    A2 = np.cos(X_xi) * vex * (T2 - (tau - T2)*np.log(tau/(tau - T2)))
    B2 = np.cos(X_xi) * -1* vex * (((T2**2)/2) - (tau*(T2 - (tau-T2) * np.log(tau/(tau - T2)))))
    
    C2 = pstateg[1] - etat + vstateg[1]*T2 - ((0.5*gstar*T2**2) * np.cos(phis)) + \
        np.sin(X_xi) * (vex * (T2 - ((tau - T2) * np.log(tau/(tau - T2)))))
        
    K1 = B1 * C2/(A2*B1 - (A1*B2))
    K2 = A1*K1/B1
    
    # STEERING ANGLE IN THE INERTIAL FRAME 
    
    X = X_xi - K1 - phiT + K2*0.1
    #print(T2)
    
    # DEBUG ARRAY
    
    dbg = np.array([A1, B1, A2, B2, C2, K1, K2, X_xi, dxi, deta, T2, X])
    
    
    return T2, X, dbg






















    