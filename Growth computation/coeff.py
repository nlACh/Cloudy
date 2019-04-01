import constants as cc

def KA(T):
    A = (2*cc.sigma_w*cc.M_w)/(cc.rho_w*cc.R_v*T) # all in SI
    return A

def KB(M_species, rho_species, i):
    B = (i*cc.M_w*rho_species)/(M_species*cc.rho_w)
    return B

def C(T, e_s):
    K_ = ((cc.L_vap**2)*cc.rho_w)/(cc.kappa*cc.R_v*(T**2))
    D_ = (cc.rho_w*cc.R_v*T)/(e_s*cc.D_v)
    C = K_ + D_
    return C

def A(T):
    a1 = (cc.L_vap*cc.g)/(cc.R_v*cc.cp*T)
    a2 = cc.g/cc.R
    A = (1/T)*(a1 - a2)
    return A

def B(T, e_s, p):
    b1 = (cc.R_v*T)/(cc.epsilon*e_s)
    b2 = (cc.epsilon*(cc.L_vap**2))/(p*T*cc.cp)
    B = cc.rho_air*(b1 + b2)
    return B