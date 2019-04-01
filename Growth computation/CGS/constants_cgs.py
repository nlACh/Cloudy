# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 19:20:55 2019

@author: nlpl9
"""
#TODO convert all to cgs
g = 981 # cm/s^2

lapse_rate_env = 0.004/100 # Kelvin per (centi)meter 
dry_gas_const = 8.314 # SI J/Kg K mol
R = dry_gas_const
R_moist = 461.5 # some bizzare calculation
T = 277 #Kelvin, 4C
rho_w = .99986 #K/m^3
#M_w = 18.05 #Kg
#sigma_w = 72E-3 #surface tension @ 4C b/w water droplet and vapour
L_vap = 2.5E6/(4.2*1E3) # cal/g Latent heat of vaporisation
#rho_dry = 1.29 # vapor density of dry air
eps = 0.622 # some mixing ratio thingies
#cp = 1004 # J/(Kg K) Specific heat a const pressure
kappa = (2.4E-2)/(4.2*100) # cal/(cm s K) Thermal conductivity of air
D = 2.21E-5*1E4 # cm2/s Coefficient of diffusion of air molecules