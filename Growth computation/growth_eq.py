import numpy as np
import matplotlib.pyplot as plt

import constants as c
import coeff as co


dt = 0.05 # s
###
# Data for species... Using NaCl
M_species = 58.44E-3 # Kg/mol
rho_species = 2.17E3 # Kg/m3
i = 2 # Vant Hoff Factor
rs = 100E-9 # m
#
###

def es(x):
    a = [5.4842763E1, 6.76322E3, 4.21, 3.67E-4, 4.15E-2, 2.188E2, 5.3878E1, 1.33122E3, 9.44523, 1.4025E-2]
    tmp = a[0] - a[1]/x - a[2]*np.log(x) + a[3]*x + np.tanh(a[4]*(x-a[5]))*(a[6] - a[7]/x - a[8]*np.log(x) + a[9]*x)
    tmp1 = np.exp(tmp)
    return tmp1


def profiles(base_temp, base_p, time, updraft):
    height = [updraft*x for x in time]
    power = c.g/(c.R_v*c.adia_lapse_rate)
    p = [(1 - c.adia_lapse_rate*updraft*x/base_temp) for x in time]
    pressure = [base_p*(x**power) for x in p]
    temp = [(base_temp - c.adia_lapse_rate*updraft*x) for x in time]
    
    e_s = [es(x) for x in temp]
    
    return height, pressure, temp, e_s

r_ini = co.r0(1E-15, rho_species)
r_1 = co.r0(1E-16, rho_species)
r_2 = co.r0(1E-17, rho_species)
_t = np.linspace(0, 1000, int(1000/dt)+1) 
deltaz, p, T, es1 = profiles(283, 8E4, _t, updraft = 0.1)


def integrator(S0, r0, num, time, m_dry_parcel, updraft, T, es, p):
    S=[S0]
    r=[r0]
    KB = co.KB(M_species, rho_species, i)
    for x in range(0, len(time)-1):
        KA = co.KA(T[x])
        A = co.A(T[x])
        B = co.B(T[x], es[x], p[x])
        C = co.C_(T[x], es[x])
        
        _rdt = (1/r[x])*(S[x]- KA/r[x] + KB*(rs/r[x])**3)*(1/C)
        r.append(r[x] + _rdt*dt)
        
        B1 = B * (4 * np.pi * num * c.rho_w * (r[x]**2))/m_dry_parcel
        
        _sdt = A * updraft - B1 * _rdt
        S.append(S[x] + _sdt)
    
    return S, r
    """    
    plt.grid()
    plt.semilogx(time, S)
    plt.xlabel('time')
    plt.ylabel('Supersat %')
    plt.show()
    plt.grid()
    #plt.ylim(4.95E-8, 5.05E-8 )
    plt.semilogx(time, r)    
    plt.xlabel('time')
    plt.ylabel('radius m')
    plt.show()
    """
rs = r_ini
S, r = integrator(0.005, r_ini, 1E8, _t, 1, 0.1, T, es1, p)
rs = r_1
S1, r1 = integrator(0.005, r_1, 1E8, _t, 1, 0.1, T, es1, p)
rs = r_2
S2, r2 = integrator(0.005, r_2, 1E8, _t, 1, 0.1, T, es1, p)

plt.grid()
plt.semilogx(_t, S, '-k')
plt.semilogx(_t, S1, '-r')
plt.semilogx(_t, S2, '-b')
plt.xlabel('time')
plt.ylabel('Supersat %')
plt.show()

plt.grid()
plt.semilogy(r, _t, '-k')
plt.semilogy(r1, _t, '-r')
plt.semilogy(r2, _t, '-b')
plt.ylabel('time')
plt.xlabel('radius')
plt.show()