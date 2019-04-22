import constants as cc
import numpy as np

KA = lambda T: (2*cc.sigma_w*cc.M_w)/(cc.rho_w*8.314*T) # Taking Kg/mol. So R_v is in J/K mol Use 13.369

KB = lambda M_species, rho_species, i: (i*cc.M_w*rho_species)/(M_species*cc.rho_w)

F_k = lambda T: ((cc.L_vap**2)*cc.rho_w)/(cc.kappa*cc.R_v*T)
F_d = lambda T, e_s: (cc.rho_w*cc.R_v*T)/(e_s*cc.D_v)
C_ = lambda T, e_s: F_k(T) + F_d(T, e_s)

A = lambda T: (((cc.epsilon*cc.L_vap)/(cc.cp*T)-1)*cc.g)/(cc.R_dry*T)
B = lambda T, e_s, p: (cc.epsilon*(cc.L_vap**2))/(cc.cp*cc.R_dry*(T**2)) + p/(cc.epsilon*e_s)

r0 = lambda wt, rho: (wt/((4/3)*np.pi*rho))**(1/3)

"""
def B(T, e_s, p):
    b1 = (cc.R_v*T)/(cc.epsilon*e_s)
    b2 = (cc.epsilon*(cc.L_vap**2))/(p*T*cc.cp)
    B_ = cc.rho_air*(b1 + b2)
    return B_
"""