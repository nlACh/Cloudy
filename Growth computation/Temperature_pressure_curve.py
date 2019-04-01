# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 17:26:35 2019

@author: nlpl9
"""
import numpy as np
import matplotlib.pyplot as plt

import Kohler as k
import constants as const
#import MurphyKoop_poly_estimation as MK_poly

t_max = 3600 #s ie, 1 hour
dt = 0.01 #s
# Thus,
_t = np.linspace(0, 200, int(200/dt -1))

S_amb = 0.03

u = 0.0

def es(x):
    a = [5.4842763E1, 6.76322E3, 4.21, 3.67E-4, 4.15E-2, 2.188E2, 5.3878E1, 1.33122E3, 9.44523, 1.4025E-2]
    tmp = a[0] - a[1]/x - a[2]*np.log(x) + a[3]*x + np.tanh(a[4]*(x-a[5]))*(a[6] - a[7]/x - a[8]*np.log(x) + a[9]*x)
    tmp1 = np.exp(tmp)
    return tmp1


def atmos_init(base_temp, updraft_velocity, base_pressure):
    global u
    u = updraft_velocity
    height_profile = [updraft_velocity*x for x in _t]
    temp_profile = [(base_temp - const.lapse_rate_env*x) for x in height_profile]
    power = const.g*(const.dry_gas_const*const.lapse_rate_env)
    pressure_profile = [base_pressure*(x/base_temp)**power for x in temp_profile]
    sat_vap_press = [es(x) for x in temp_profile]
    return height_profile, temp_profile, pressure_profile, sat_vap_press

hp, tp, pp, esT = atmos_init(283, 1, 9E4)
A,B,S_c,r_c = k.kohler_init(.05844, 2.17E3, 2, 50E-9)


def euler(S_amb, r_seed, u, m_dry_parcel):
    S = [S_amb]
    r = [r_seed]
    for t in range(len(_t)-1):
        a1 = (1/tp[t])*((const.L_vap*const.g)/(const.R_moist*const.cp*tp[t]) - const.g/const.dry_gas_const)
        a2 = const.rho_dry*((const.R_moist*tp[t])/(const.eps*esT[t]) + (const.eps*(const.L_vap**2))/(pp[t]*tp[t]*(const.cp)))
        # TODO:
        # Expressed sat_vap_press as fn of T. Employ murphy koop
        Fk = (const.rho_w*const.L_vap**2)/(const.kappa * const.R_moist * tp[t]**2)
        
        Fd = (const.rho_w * const.R_moist * tp[t])/(esT[t] * const.D)
        _rdt = (1/r[t])*(S[t] - 1)/(Fk + Fd)
        r.append(r[t] + _rdt*dt)
        
        _sdt = a1*u - ((4*np.pi*const.rho_w*r[t]*2)/m_dry_parcel)*a2
        S.append(S[t] + dt*_sdt)
    plt.grid()
    plt.ylim(-1, 0.05)
    plt.plot(_t, r)
    plt.show()
    plt.plot(_t, S)
    plt.show()
    return S, r
s,r = euler(-10, 5E-8, u, 0.6)    
### Now consider the growth of droplets.
