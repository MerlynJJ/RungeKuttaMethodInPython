#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:13:12 2022

@author: merlyn
"""

import numpy as np
#these two integrator functions throw one value at a time

# first order differential equations INTEGRATOR

def RKN4_1fstODE(f,r,ti,tf,h,params):
    # initializing where to store the results 
    R = np.zeros((len(r),len(time)), dtype = 'float') 
    # fst index cartesian coordinate
    # scnd index the element in such coordinate
    R[:,0] = r
    #generate the new values 
    for t in np.arange(ti,tf,h)
        # slope at the beginning of the time step
        k1 = f(r,t,params)
        # use the slope k1 to step halfway through the time step
        r1 = r + k1*h/2
        # an estimate of the slope at the midpoint
        k2 = f(r1,t+(h/2),params)
        # use the slope k2 to step halfway through the time step
        r2 = r + k2*h/2
        # another estimate of the slope at the r2 midpoint
        k3 = f(r2,t+(h/2),params)
        # use the slope, k3, to step all the way across the time step (to t₀+h)
        r3 = r + k3*(h/2)
        # an estimate of the slope at the endpoint
        k4 = f(r3, t + h,params)
        # use a weighted sum of these slopes to get our final estimate of r*(t₀+h)
        R[:, n+1] = r + (h/6)*(k1+2*k2+2*k3+k4)
        # update the r value
        r = R[:, n+1]
    return R


# Second order differential equations INTEGRATOR

def RKN4_2ndODE(f,r,v,ti,tf,h,params):
    # initializing where to store the results
    R = np.zeros((len(r),len(time)), dtype = 'float') 
    V = np.zeros((len(v),len(time)), dtype = 'float') 
    # fst index cartesian coordinate
    # scnd index element in such coordinate
    R[:,0] = r
    V[:,0] = v
    #generate the new values 
    for n in np.arange(ti, tf+h, h):
        # slope at the beginning of the time step
        k1 = f(r,v,t[n+1],params)
        # use the slope k1 to step halfway through the time step
        v1 = v + k1*h/2
        r1 = r + (h/2)*((v+v1)/2)
        # an estimate of the slope at the midpoint
        k2 = f(r1,v1,t[n+1]+(h/2),params)
        # use the slope k2 to step halfway through the time step
        v2 = v + k2*h/2
        r2 = r + (h/2)*((v+v2)/2)
        # another estimate of the slope at the midpoint
        k3 = f(r2,v2,t+(h/2),params)
        # use the slope, k3, to step all the way across the time step (to t₀+h)
        v3 = v + k3*(h/2)
        r3 = r + h*((v+v3)/2)
        # an estimate of the slope at the endpoint
        k4 = f(r3, v3, t[n+1] + h,params)
        # use a weighted sum of these slopes to get our final estimate of y*(t₀+h)
        V[:, n+1] = v + (h/6)*(k1+2*k2+2*k3+k4)
        R[:, n+1] = r + (h/6)*(v1+2*v2+2*v3+V[:,n+1])
        # update the r and v values
        r = R[:, n+1]
        v = V[:, n+1]
    return R, V
