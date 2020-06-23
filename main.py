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
        
     #Call function to compute vehicle forces and mass (T, D, & M):
    
     #Determine flight path angle and divide forces into their x/y components:
    
    
     #If conditions for flight stages, initial launch, pitch kick, gravity turn
        
     #IGM ROUTINE GOES HERE
     rx = STATE[0]
     ry = STATE[1]
     vx = STATE[2]
     vy = STATE[3]
     
     R = (rx**2 + ry**2)**0.5
     V = (vx**2 + vy**2)**0.5
     gx = 0
     gy = 9.81
     #Setup equations
     F9S1 = stages[0]
    
     M = F9S1.M - (F9S1.mdot*t)   
     #print("Solution time: ", t)
     
     print("Mass: ", M)
     Tx = 0
     Ty = F9S1.Thrust
     Dx = 0
     #Dy = 0.5 * 1.225 * V**2 * (F9S1.A) * F9S1.Cd
     Dy = F9S1.GetDrag(R-6378e3,V)

     rx = vx
     ry = vy
    
     vx = Tx/M - Dx/M + gx
     vy = Ty/M - Dy/M - gy
    
    #Return numpy array of the equations
     #print("Vertical velocity: ", ry)
     return np.array([rx, ry, vx, vy])
        
       

def main():
    
    """
        Main program routine
        
        Initialise vehicle configuration, initial conditions, and call 
        ODEINT routine to integrate the trajectory
        
    """

    #DEFINE FALCON 9:
    
    F9S1 = Stage("Stage 1", 7426e3, 0.3, 10.75, 282, 27.2e3, 411e3, 116e3)
    F9S2 = Stage("Stage 2", 934e3, 0.3, 10.75, 348, 4.5e3, 111.5e3, 0)
    
    
    #DEFINE INITIAL CONDITIONS:
    # UNITS ARE M AND M/S
    
    x0 = 0
    y0 = 6378e3
    vx0 = 0
    vy0 = 0
    
    STATE_0 = np.array([x0, y0, vx0, vy0])
    
    #Integration time:
    
    Tf = F9S1.tb
    
    T_array = np.arange(0,Tf,1)
    
    Stages = [F9S1, F9S2]
    
    sol = odeint(EqOfM, STATE_0, T_array, args = (Stages,))
    
    print("End")
    
    rx = sol[:,0]
    ry = sol[:,1]
    vx = sol[:,2]
    vy = sol[:,3]
    
    h = (rx**2 + ry**2)**0.5 - 6378e3
    
    plt.figure()
    plt.plot(T_array, vy)
    plt.figure()
    plt.plot(T_array, h)
    plt.title("Altitude")
    
    
    
    
    
if __name__ == "__main__":

    main()
    
    
    
    
    
    
    
    
    




