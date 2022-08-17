#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:13:12 2022

@author: merlyn
"""

# first order differential equations
# Second order differential equations
def RKN4_1fstODE(f,y,t,h,params):
    # slope at the beginning of the time step
    k1 = f(y,t,params)
    # use the slope k1 to step halfway through the time step
    y1 = y + k1*h/2
    # an estimate of the slope at the midpoint
    k2 = f(y1,t+(h/2),params)
    # use the slope k2 to step halfway through the time step
    y2 = y + k2*h/2
    # another estimate of the slope at the midpoint
    k3 = f(y2,t+(h/2),params)
    # use the slope, k3, to step all the way across the time step (to t₀+h)
    y3 = y + k3*(h/2)
    # an estimate of the slope at the endpoint
    k4 = f(y3, t + h,params)
    # use a weighted sum of these slopes to get our final estimate of y*(t₀+h)
    V = y + (h/6)*(k1+2*k2+2*k3+k4)
    return V


# Second order differential equations
def RKN4_2ndODE(f,y,v,t,h,params):
    # slope at the beginning of the time step
    k1 = f(y,v,t,params)
    # use the slope k1 to step halfway through the time step
    v1 = v + k1*h/2
    y1 = y + (h/2)*((v+v1)/2)
    # an estimate of the slope at the midpoint
    k2 = f(y1,v1,t+(h/2),params)
    # use the slope k2 to step halfway through the time step
    v2 = v + k2*h/2
    y2 = y + (h/2)*((v+v2)/2)
    # another estimate of the slope at the midpoint
    k3 = f(y2,v2,t+(h/2),params)
    # use the slope, k3, to step all the way across the time step (to t₀+h)
    v3 = v + k3*(h/2)
    y3 = y + h*((v+v3)/2)
    # an estimate of the slope at the endpoint
    k4 = f(y3, v3, t + h,params)
    # use a weighted sum of these slopes to get our final estimate of y*(t₀+h)
    V = v + (h/6)*(k1+2*k2+2*k3+k4)
    Y = y + (h/6)*(v1+2*v2+2*v3+V)
    return k4,Y,V

