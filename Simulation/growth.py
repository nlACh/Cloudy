#Simulation being done for a uniform distribution of aerosol particles
#with starting radius r_s. Thus r0 = r_s (!)
#Species chosen: NaCl
import numpy as np
import matplotlib.pyplot as plt

from constants import *

U = 0.1 # m/s, updraught velocity

def profiler(T0,p0):
    global _t,_h,_T,_p
    _h = [U*t for t in _t] # m, height profile
    _T = [(T0 - malr*h) for h in _h] # K, temperature profile
    pwr = g/(R_d*malr)
    _p = [p0*(1 - (malr*h)/T0)**pwr for h in _h] # Pa, pressure profile

def integrator(r0,S0,t_f,dt,T0,p0,n):
    global _t,_h,_T,_p
    steps = int(t_f/dt)
    _t = np.linspace(0.0,t_f,steps+1)
    profiler(T0,p0)
    _r = np.zeros(steps+1).tolist()
    _S = np.zeros(steps+1).tolist()
    _S1 = np.zeros(steps+1).tolist()
    _drdt = np.zeros(steps+1).tolist()
    _r[0] = r0
    _S[0] = S0
    _S1[0] = S0
    print(B)
    drdt1 = dSdt1 = 0.0
    #Simulation
    for t in range(0,steps):
        T = _T[t]
        p = _p[t]
        #radius growth
        drdt = (_S[t] - A(T)/_r[t] + B/_r[t]**3)/(_r[t]*F(T,p))
        _drdt[t] = drdt
        _r[t+1] = _r[t] + dt*drdt
        #SS growth
        tmp = 4*np.pi*rho_w*n*_r[t]**2/M_d
        dSdt = P(T)*U - tmp*C(T,p)*drdt
        _S[t+1] = _S[t] + dt*dSdt
        _S1[t+1] = _S1[t] + dt*P(T)*U
        if(t==1):
            print(r_c(T),S_c(T),A(T)/_r[t] - B/_r[t]**3,P(T),C(T,p),tmp*C(T,p),tmp*C(T,p)*drdt,drdt,dSdt)
        if((drdt1>0 and drdt<0) or (drdt1<0 and drdt>0)):
            print("drdt flipped sign!")
        if((dSdt1>0 and dSdt<0) or (dSdt1<0 and dSdt>0)):
            print("dSdt flipped sign! S =",_S[t+1])
        drdt1 = drdt
        dSdt1 = dSdt
        
    plt.plot(_t,_r)
    plt.show()
    plt.xlim(1,5000)
    plt.ylim(0.0,1.0e-8)
    plt.plot(_t,_drdt)
    plt.show()
    plt.plot(_t,_S)
    plt.show()
    
integrator(r_s,0.005,5000,0.05,283,8.0e+4,1.0e+8)
        