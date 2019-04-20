#Contains the necessary constants of the growth process
#Temperature range: -30 to +30 deg C
import numpy as np

g = 9.8 # m/s^2, accl. due to gravity
T_a = 273.15 # K, 0 deg C
p_a = 1.01325e+5 # Pa, std. atmos. pressure, 1 atm

R_v = 461.5 # J/Kg/K, water vapour gas constant
R_d = 287.1 # J/Kg/K, dry air gas constant
m_w = 18.015e-3 # Kg/mol, molar mass of water
m_d = 28.96e-3 # Kg/mol, eff. mol. mass of dry air
eps = m_w/m_d # mass ratio
st = 75.86e-3 # N/m, surface tension of water
rho_w = 999.860 # Kg/m^3, density of water at 0.30 deg C
rho_d = 1.275 # Kg/m^3, density of air at 0 deg C
c_p = 1.0038e+3 # J/Kg/K, sp. heat cap. at 2 deg C

M_d = 1.0 #Kg, mass of dry air in the parcel

#Data for NaCl
m_s = 58.44e-3 # Kg, molecular weight [PubChem]
rho_s = 2.17e+3 # Kg/m^3, density at ~4 deg C
i = 2 # van't Hoff factor, generic
r_s = 200.0e-9 # m, dry radius
B = i*m_w*rho_s*r_s**3/(m_s*rho_w) # Kohler hygroscopicity coefficient

A = lambda T: 2*st/(rho_w*R_v*T) # Kelvin coefficient

r_c = lambda T: np.sqrt(3*B/A(T))
S_c = lambda T: (2/np.sqrt(B))*(A(T)/3)**1.5

dalr = g/c_p # K/m, dry adiabatic lapse rate
malr = 4.0e-3 # K/m, moist adiabatic lapse rate

e_s = lambda T: 611.2*np.exp(17.67*(T-273)/(T-29.5)) # Pa, sat. vap. pressure, Bolton 1980

L = lambda T: -2403.43*T + 3.1577e+6 # J/Kg, latent heat of vaporisation, fitted

K = lambda T: 4.184e-3*(1.0465 + 0.017*T) # W/m/K, thermal conductivity of air at 0 deg C

D = lambda T,p:21.1e-6*((T/T_a)**1.94)*(p_a/p) # m^2/s, coeff. of diffusion of vap. in air

F_k = lambda T: (rho_w*L(T)**2)/(K(T)*R_v*T) # thermal term in TME
F_d = lambda T,p: (rho_w*R_v*T)/(D(T,p)*e_s(T)) # diffusion term in TME
F = lambda T,p: F_k(T) + F_d(T,p)

P = lambda T: (((eps*L(T))/(c_p*T) - 1)*g)/(R_d*T) # production coefficient
C = lambda T,p: (eps*L(T)**2)/(c_p*R_d*T**2) + p/(eps*e_s(T)) # condensation coefficient