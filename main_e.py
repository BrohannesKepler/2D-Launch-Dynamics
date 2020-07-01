#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:03:47 2020

@author: haider
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import atmosphere

class Stage():
    
    """ Define a vehicle class that specifies the number of stages, mass of each 
        stage, thrust, specific impulse, drag expression
        
    """
    def __init__(self, Name, Thrust, Cd, A, Isp, Ms, Mf, Mp):
        
        self.Name = Name
        self.Thrust = Thrust
        self.Cd = Cd
        self.A = A
        self.Isp = Isp
        self.Ms = Ms 
        self.Mf = Mf 
        self.Mp = Mp
        
        #Compute the mass flow rate and burn time for the stage
        self.mdot = self.Thrust/(self.Isp * 9.81)
        self.tb = self.Mf/self.mdot
        self.M = self.Ms + self.Mf + self.Mp
        
        
    def GetDrag(self, h, V):
        """ Function to compute the drag force, modelled on an
            Atlas V launch vehicle as a function of Mach no."""

        #CONSTANTS:
        gamma = 1.14
        R = 287
        
        T, rho = atmosphere.Atmosphere(h)
        
        A = np.sqrt(gamma*R*T)
        M = V/A
        
        if M <= 1.25:
            
            Cd = -0.0415*M**3 + 0.3892*M**2 - 0.2614*M + 0.303
            
        elif M > 1.25 and M <= 4:
            
            Cd = -0.049*M**4 + 0.5664*M**3 - 2.3265*M**2 + 3.8512*M - 1.6625

        elif M > 4 and M <= 10:
            
            Cd = -0.0037 * M**3 + 0.0695*M**2 - 0.4105*M + 0.9732
            
        elif M > 10:
            Cd = 0.255
            
        D = 0.5 * rho * V**2 * self.A * Cd

        return D

def EqOfM(STATE, t, stages):
    
    
     """This function defines the 2D cartesian equations of motion 
     for the launch vehicle."""
        
     # DECLARE CONSTANTS
     mu = 3.986E14

    # PARSE STATE VARIABLES FOR USE 

     rx = STATE[0]
     ry = STATE[1]
     vx = STATE[2]
     vy = STATE[3]
     
     R = (rx**2 + ry**2)**0.5
     V = (vx**2 + vy**2)**0.5
     
     gx = -mu * rx/(rx**2 + ry**2)**1.5
     gy = -mu * ry/(rx**2 + ry**2)**1.5
     
     #Setup equations
     F9S1 = stages[0]
     F9S2 = stages[1]
     

     if t <= F9S1.tb:    
         M = F9S1.M - (F9S1.mdot*t)   
     elif t > F9S1.tb:
         M = F9S2.M - (F9S2.mdot*(t-F9S1.tb))
     

     Tx = 0
     Ty = F9S1.Thrust
     Dx = 0
     Dy = F9S1.GetDrag(R - 6378e3, V)
     print(t, M)
     rx = vx
     ry = vy
    
     vx = Tx/M - Dx/M + gx
     vy = Ty/M - Dy/M + gy
    
    #Return numpy array of the equations
     #print("Vertical velocity: ", ry)
     return np.array([rx, ry, vx, vy])
        

def euler(STATE0, t0, dt, Tf, stages):
    
    """
    Euler's method for integrating the system of equations. First order accurate
    """

    #Pre-allocation of state variables
    
    steps = len(np.arange(t0, Tf, dt))
    
    X = np.zeros((steps+1))    
    Y = np.zeros((steps+1))    
    Vx = np.zeros((steps+1))    
    Vy = np.zeros((steps+1))    
    T_array = np.zeros((steps+1))
    
    #Initial condition allocation
    X[0] = STATE0[0]
    Y[0] = STATE0[1]
    Vx[0] = STATE0[2]
    Vy[0] = STATE0[3]
    
    t = t0
    n = 0
    
    #While the time is less than the final time, call the function and step
    
    while t < Tf:
      
        STATE = np.array([X[n], Y[n], Vx[n], Vy[n]])
        S = EqOfM(STATE, t, stages)
        
        X[n+1] = X[n] + dt * S[0]
        Y[n+1] = Y[n] + dt * S[1]
        Vx[n+1] = Vx[n] + dt * S[2]
        Vy[n+1] = Vy[n] + dt * S[3]
        
        t = t + dt
        T_array[n + 1] = t            
        n = n + 1
    
    #Prescribe state vectors to the state matrix 
    STATEF = np.zeros((steps+1, 5))
    STATEF[:, 0] = T_array
    STATEF[:, 1] = X
    STATEF[:, 2] = Y
    STATEF[:, 3] = Vx
    STATEF[:, 4] = Vy

    return STATEF
    
def main():
    
    """
        Main program routine
        
        Initialise vehicle configuration, initial conditions, and call 
        ODEINT routine to integrate the trajectory
        
    """

    #DEFINE FALCON 9:
    
    F9S1 = Stage("Stage 1", 7426e3, 0.3, 10.75, 282, 27.2e3, 411e3, 116e3)
    F9S2 = Stage("Stage 2", 934e3, 0.3, 10.75, 348, 4.5e3, 111.5e3, 0)
    
    
    # DEFINE INITIAL CONDITIONS:
    # UNITS ARE M AND M/S
    
    x0 = 0
    y0 = 6378e3
    vx0 = 0
    vy0 = 0
    STATE_0 = np.array([x0, y0, vx0, vy0])
    
    #Integration time:
    
    Tf = F9S1.tb
    Stages = [F9S1, F9S2]

    # CALL INTEGRATION ROUTINE 
    sol = euler(STATE_0, 0, 0.5, Tf, Stages)
    
    print("End")
    

    t = sol[:, 0]
    rx = sol[:, 1]
    ry = sol[:, 2]
    vx = sol[:, 3]
    vy = sol[:, 4]
    
    h = (rx**2 + ry**2)**0.5 - 6378e3
    
    
    plt.figure()
    plt.plot(t, vy)
    plt.figure()
    plt.plot(t, h)
    plt.title("Altitude")
    
    
    #plt.plot(T_array, vy, t2, vy2, 'rx')
    
    
    
if __name__ == "__main__":

    main()
    
    
    
    
    
    
    
    
    




