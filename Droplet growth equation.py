# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np

# Time step
dt = .01 #in seconds

# Parameter
epsilon = 1.0E-20

"""
Kohler theory stuff
"""

# Some constants:
R = 8.314
T = 277 #Kelvin, 4C
rho_w = 1.00E3 #K/m^3
M_w = 18.05E-3 #Kg
sigma_w = 72E-3 #surface tension @ 4C b/w water droplet and vapour

# Some declared stuff, need to be user fed
n=02.0 #dissociation co-efficient
M_species=0.0
rho_species=0.0
species = ''
# Kohler hygroscopicity co-efficient
a=0.0 # Kelvin term
b=0.0 # Raoult's law term

"""
Consider the starting equations...
rdr/dt = (S-1 - a/r +b/r3)/(K + D)
and
ds/dt = A1*V - A2*dw/dt
w is some function of r. So replace dw/dt with dr/dt
"""

# Init values of A1 and A2
A1 = 0.0
A2 = 0.0

# -------------------------------------------------------------------------------------------- #
def init():
    # A function to take in the species of molecules considered
    global rho_species, M_species, n, a, b, species
    #print("What kind of aerosol species would you be using?")
    #species = input()
    #print("Enter some info about the species: molecular mass, density and dissociation coefficient")
    #M_species, rho_species, n = [np.float(x) for x in input().split(',')]
    M_species = 58.440E-3
    rho_species = 2.170E3
    # calculate coefficient of Kohler theory
    a = (2*sigma_w*M_w)/(rho_w*R*T)
    b = (n*M_w*rho_species)/(M_species*rho_w)
    print("A", a)
    print("B", b)
    # TODO
    # Add initializers for A1 and A2

def Kohler(a, b, r, r_dry):
    S = [100*(a/rad - b*(r_dry**3)/(rad**3)) for rad in r]
    sat_ratio = [(1+x) for x in S]
    #print(S)
    return S, sat_ratio


def integrator(S_i, r_i, step, K, D, V):
    ri_t = [] #need to move somewhere else
    si_t = []
    # calculate time derivative of r
    t = 0
    while(abs(si_t[t]-si_t[t-1]) > epsilon):
        _r = S_i/((K+D)*r_i)
        ri_t.append( ri_t[t] + _r * dt)

        _s = A1*V - A2*_r*(ri_t[t+1]**2)
        si_t.append( si_t[t] + _s * dt)
    return ri_t, si_t


def main():
    init()
    # consider an uniform distribution of 1000 particles
    R_T = []
    S_T = []
    #print("Enter min and max value of seed radius")
    #seed_small, seed_high = [np.float(x) for x in input().split(',')]
    r = np.random.uniform(1E-7, 1E-5, size=1000)
    # Arrange the array in ascending order of radius
    r.sort()
    #print(r)
    S, sat_ratio = Kohler(a,b,r, r_dry = 50E-9)
    plot(r, S, 'supersat %', [-0.3, 0.3])
    plot(r, sat_ratio, 'sat ratio', [0.8, 1.3])
    """
    for i in range(1000):
        ri_t, si_t = integrator(S[i], r[i], dt, 1, 1, 1)
        R_T.append(ri_t)
        S_T.append(si_t)
    """

def plot(x,y,string, ylim):
    # TODO
    plt.grid()
    plt.ylim(ylim[0], ylim[1])
    #plt.autoscale(enable=True, axis='both', tight=True)
    plt.xlabel('radius')
    plt.ylabel(string)
    plt.semilogx(x,y)
    plt.show()
    """
    Some way to just see the change of radius and supersaturation over time for individual droplets...
    """

main()
