# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 03:50:49 2019

@author: nlpl9
An estimation of e_s(T) from clausius clayperon equation...
"""

import numpy as np
import matplotlib.pyplot as plt

a = [5.4842763E1, 6.76322E3, 4.21, 3.67E-4, 4.15E-2, 2.188E2, 5.3878E1, 1.33122E3, 9.44523, 1.4025E-2]
e_s = []
Temp = []
TLow = 0.0
THigh = 0.0

step = 0.1 # Kelvin

def init():
    global TLow, THigh
    TLow, THigh = [np.float(x) for x in input("Enter initial and final temperature range").split(',')]

def es(x):
    global e_s
    tmp = a[0] - a[1]/x - a[2]*np.log(x) + a[3]*x + np.tanh(a[4]*(x-a[5]))*(a[6] - a[7]/x - a[8]*np.log(x) + a[9]*x)
    tmp1 = np.exp(tmp)
    e_s.append(tmp1)
"""
def main():
    init()
    T = TLow
    while(T<=THigh):
        es(T)
        Temp.append(T-273)
        T += step
    
    plt.plot(Temp, e_s)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xlabel('Temperature(C)')
    plt.ylabel('Saturation vapor pressure(Pa)')
    plt.show()
    
main()
"""