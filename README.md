# 2D-Launch-Dynamics

Simple script that solves the equations of motion of a launch vehicle, modelled as a SpaceX Falcon 9 v1.2 "Full Thrust" , in an Earth centered cartesian system, with the y axis
passing through the launch site and the x-axis aligned with the launch direction. A simple Euler integrator in use giving first order temporal and spatial accuracy.

Drag model is based on a 3rd/4th order polynomial fit of the drag coefficient for an Atlas V launch vehicle.

Vehicle will ascend to ~ 100 m before initiating a pitch-kick and then executing a gravity turn until burnout of stage 1, by following the velocity vector.

The second stage will then utilise the Iterative Guidance Mode (IGM) routine used on the Saturn V for upper stage targeting into the desired orbit.
