# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 19:33:50 2019

@author: nlpl9
"""

g = 9.81 # m/s^2

lapse_rate_env = 0.004 # Kelvin per meter 
dry_gas_const = 8.314
R = dry_gas_const
R_moist = 461.5 # some bizzare calculation
T = 277 #Kelvin, 4C
rho_w = 1.00E3 #K/m^3
M_w = 18.05E-3 #Kg
sigma_w = 72E-3 #surface tension @ 4C b/w water droplet and vapour
L_vap = 2.5E6 # J/Kg Latent heat of vaporisation
rho_dry = 1.29 # vapor density of dry air
eps = 0.622 # some mixing ratio thingies
cp = 1004 # J/(Kg K) Specific heat a const pressure
kappa = 2.4E-2 # J/(msK) Thermal conductivity of air
D = 2.21E-5 # m2/s Coefficient of diffusion of air molecules