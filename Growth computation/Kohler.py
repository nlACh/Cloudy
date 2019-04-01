# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 17:52:25 2019

@author: nlpl9
"""
#import numpy as np
### Defining the critical radii and supersat of Kohler theory.

## ccants...
import constants_cgs as cc

def kohler_init(M_species, rho_species, n, rs):
    P = (2*cc.sigma_w*cc.M_w)/(cc.rho_w*cc.R*cc.T)
    Q = (n*cc.M_w*rho_species)/(M_species*cc.rho_w)
    #S_c = (2*(A/3*rs)**1.5)/(B**0.5)
    #r_c = np.sqrt((3*B*rs**3)/A)
    
    return P, Q

