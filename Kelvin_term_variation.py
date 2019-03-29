# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:28:00 2019

@author: nlpl9
"""
"""
Kohler theory stuff
"""
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np
# Some constants:
R = 8.314
T = 277 #Kelvin, 4C
rho_w = 1.00E3 #K/m^3
M_w = 18.05E-3 #Kg
sigma_w = 72E-3 #surface tension @ 4C b/w water droplet and vapour

# Some declared stuff, need to be user fed
n=02.0 #dissociation co-efficient
M_species=58.440E-3 #Sodium Chloride
rho_species=2.17E-3 #Kg/m3
species = 'Sodium chloride'


# Kohler hygroscopicity co-efficient
a = (2*sigma_w*M_w)/(rho_w*R)
b = (n*M_w*rho_species)/(M_species*rho_w)
_T = np.linspace(230, 350, 100).tolist()
r = np.linspace(1.0E-7, 1.0E-5, 1001).tolist()
r.sort()

def plot_kelvin_r():
    a1 = a/T
    s = [np.exp(a1/x) for x in r]
    plt.grid()
    plt.semilogx(r, s)
    plt.xlabel('radius in m')
    plt.ylabel('Saturation ratio (s)')
    plt.show()
    
    
def plot_kelvin_t():
    A = [(a/x) for x in _T]
    r = [5.0E-8, 1.0E-7, 5.0E-7, 1.0E-6]
    s0 = [np.exp(y/r[0]) for y in A]
    s1 = [np.exp(y/r[1]) for y in A]
    s2 = [np.exp(y/r[2]) for y in A]
    s3 = [np.exp(y/r[3]) for y in A]
    _T1 = [(t-273) for t in _T]
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(_T1, s0, '-b', label = 'r=0.5E-7')
    ax.plot(_T1, s1, '--b', label = 'r=1.0E-7')
    ax.plot(_T1, s2, '-r', label = 'r=5.0E-7')
    ax.plot(_T1, s3, '--r', label = 'r=1.0E-6')
    #ax.axis('equal')
    ax.legend(frameon=False)
    plt.xlabel('Temperature (C)')
    plt.ylabel('Saturation ratio (s)')
    fig#.show()

def kelvin_vsTT():
    plt.grid()
    plt.plot([t-273 for t in _T], [1.0E6*a/x for x in _T], '-g')
    plt.xlabel('Temperature (C)')
    plt.ylabel(r'A in $\mu$m')
    plt.show()

    
plot_kelvin_r()
kelvin_vsTT()
plot_kelvin_t()

