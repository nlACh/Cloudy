# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np

dt = .01 #in seconds

# Parameter
epsilon = 1.0E-20

"""
Kohler theory stuff
"""
#Some specific data:
R = 8.314
T = 277 #Kelvin, 4C
rho_sulphate = 1.77
rho_w = 1.00
M_sulphate = 132
M_w = 18
n=0.8 #dissociation co-efficient
sigma_w = 72E-3 #surface tension @ 4C b/w water droplet and vapour
    
    # Kohler hygroscopicity co-efficient
b = (n*M_w*rho_sulphate)/(M_sulphate*rho_w)
a = (2*sigma_w*M_w)/(rho_w*R*T)
    
"""
Consider the starting equations...
rdr/dt = (S-1 - a/r +b/r3)/(K + D)
and 
ds/dt = aV - bdw/dt
w is some function of r. So replace dw/dt with dr/dt
"""

def Kohler(a, b, r):
    S = []
    for i in range(len(r)):
        S.append(a/r[i] - b/(r[i]**3))
    return S


def integrator(S_i, r_i, step, K, D, V):
    ri_t = [0] #need to move somewhere else
    si_t = [0]
    # calculate time derivative of r
    t = 0
    while(abs(si_t[t]-si_t[t-1]) > epsilon):
        _r = S_i/((K+D)*r_i)
        ri_t.append( ri_t[t] + _r * dt)
        
        _s = a*V - b*_r*(ri_t[t+1]**2)
        si_t.append( si_t[t] + _s * dt)
    return ri_t, si_t

def main():
    # consider an uniform distribution of 1000 particles
    R_T = []
    S_T = []
    print("Enter min and max value of seed radius")
    seed_small, seed_high = input().split(',')
    r = np.random.uniform(seed_small, seed_high, size=1000)
    S = Kohler(a,b,r)
    for i in range(1000):
        ri_t, si_t = integrator(S[i], r[i], dt)
        R_T.append(ri_t)
        S_T.append(si_t)
    
def plot():
    # TODO
    """
    Some way to just see the chnage of radius and supersaturation over time for individual droplets...
    """
    