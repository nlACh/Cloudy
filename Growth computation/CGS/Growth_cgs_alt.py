# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 19:14:16 2019

@author: nlpl9
"""

import numpy as np
import matplotlib.pyplot as plt

import constants_cgs as c
import coefficients as cf

dt = 0.1

dry_gas_SI = 8.314

def es(x):
    a = [5.4842763E1, 6.76322E3, 4.21, 3.67E-4, 4.15E-2, 2.188E2, 5.3878E1, 1.33122E3, 9.44523, 1.4025E-2]
    tmp = a[0] - a[1]/x - a[2]*np.log(x) + a[3]*x + np.tanh(a[4]*(x-a[5]))*(a[6] - a[7]/x - a[8]*np.log(x) + a[9]*x)
    tmp1 = np.exp(tmp) # tmp1 is in Pascals
    tmp1 /= 1E5 #converts to bar
    return tmp1

def init(base_temp, updraft, base_p, time): # Import all in CGS
    height_profile = [updraft*x for x in time]
    temp_profile = [(base_temp - c.lapse_rate_env*x) for x in height_profile]
    power = 9.81*(dry_gas_SI*c.lapse_rate_env*100) # unitless, but all in SI
    pressure_profile = [(base_p*(x/base_temp)**power) for x in temp_profile] # in bar
    e_s = [es(x) for x in temp_profile]
    return height_profile, temp_profile, pressure_profile, e_s, updraft



"""
Consider the equations...
dS/dt = A*updraft - 4*pi*B*n*rho*r^2*dr/dt
and...
r*dr/dt = 1/C*[S - P/r + Q/r^3]
We will be using Euler method to solve this stuff.
Finite differences might be another way. RK(4) is possible probably
"""    

def integrator(S_seed, r_seed, u, m_dry_parcel, _t, n):
    S = [S_seed]
    r = [r_seed]
    for t in range(len(_t)-1):
        A = cf.A(c.L_vap, tp[t])
        B = cf.B(c.L_vap, tp[t], pp[t], esT[t])
        C = cf.C(c.L_vap, tp[t], esT[t])
        P, Q = cf.PQ(msp = 58.44, rhosp = 2.17, i=2, rs=5E-5) #Put in required species, all cgs units
        _rdt = (1/C*r[t])*(S[t] - P/r[t] - Q/r[t]**3)
        
        
        _sdt = A*u - 4*np.pi*B*n*c.rho_w*(r[t]**2)*_rdt
        
        S.append(S[t] + _sdt * dt)
        r.append(r[t] + _rdt * dt)
    return S, r


_t = np.linspace(0, 1000, int(1000/dt)+1)
hp, tp, pp, esT, u = init(283, 10, 0.9, _t) # 10C, 10 cm/s, 900 mb
st, rt = integrator(0.1, 6E-3, 10, 0, _t, 300)
plt.grid()
plt.plot(_t, rt)
plt.xlabel('time')
plt.ylabel('radius cm')
plt.show()

plt.grid()
plt.plot(_t, [x*100 for x in st])
plt.xlabel('time')
plt.ylabel('supersat percent')
plt.show()

plt.grid()
plt.plot(rt, [x/100 for x in hp])
plt.ylabel('height')
plt.xlabel('radius cm')
plt.show()