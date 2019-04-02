# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:27:36 2019

@author: nlpl9
"""

import numpy as np
import matplotlib.pyplot as plt

es0 = 610.7
t = np.linspace(0, 1000, 100000+1)
T = [273 - 0.004*0.1*x for x in t]
delta = T[2] - T[1]
L_vap = 2.5E6
R_v = 461.5

def es(x):
    es = [es0]
    for i in range(len(x)-1):
        _es = L_vap*es[i]/(R_v*x[i]**2)
        es.append(es[i] + _es*delta)
    return es

e = es(T)
plt.grid()
plt.plot(T, e)
plt.show()
