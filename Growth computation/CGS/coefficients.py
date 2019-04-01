# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 19:32:18 2019

@author: nlpl9
"""
import constants_cgs as ccc
def A(L, T):
    tmp = (3.42E-4/T)*(2.6*L/T - 1)
    return tmp

def B(L, T, p, e_s):
    tmp = (37.8*(L/T)**2 + 1.62*p/e_s)
    return tmp

def C(L, T, e_s):
    tmp = (9.05*(L/T)**2/ccc.kappa) + 4.62E3 * (T/ccc.D)/e_s
    return tmp
    
def PQ(msp, rhosp, i, rs):
    import Kohler_cgs as k
    #Expecting all inputs in cgs
    #msp in g, rhosp in g/cc, rs in cm
    P, Q = k.kohler_init(msp*1E-3, rhosp*1E3, i, rs/100)
    # Q is dimensionless. P is returned in cm
    # Check Kohler_cgs.py
    return P, Q
    