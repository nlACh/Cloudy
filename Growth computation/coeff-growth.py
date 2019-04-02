# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:43:15 2019

@author: nlpl9
"""

# Show the variation of the terms in growth eq with temperature!
import numpy as np
import matplotlib.pyplot as plt

dt = 0.25
import coeff as co
import constants as c

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

_t = np.linspace(0, 20000, int(20000/dt)+1) 
deltaz, p, T, es1 = profiles(243, 9E4, _t, updraft = 0.5)

T1 = np.linspace(233, 303, (303-233)/0.1)
es2 = [es(x) for x in T1]
def var(p):
    A=[]
    B=[]
    for x in range(len(T1)):
        A.append(co.A(T1[x]))
        B.append(co.B(T1[x], es2[x], p))
    return A, B

A1, B1 = var(.8*101325)
A2, B2 = var(.6*101325)