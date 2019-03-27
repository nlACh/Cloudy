N = 6.022e+23 # (Avogadro number)
R = 8.314 # J/K-mol (universal gas constant)
M_w = 18.015e-3 # Kg (molecular weight of water) [PubChem]
T = 273.45 # K, = 0.30 deg C (temperature)
st = 75.86e-3 # N/m (surface tension) [Dortmund Data Bank]
rho_w = 999.860 # Kg/m^3 (density of water at T) [simetric.co.uk/si_water.htm]
A = 2*st*M_w/(rho_w*R*T) # (Kelvin coefficient)
print(A)

#data for NaCl (1)
M_s1 = 58.44e-3 # Kg (molecular weight) [PubChem]
rho_s1 = 2.17e3 # Kg/m^3 (density at ~4 deg C)
dsc = 2 # (degree of ionic dissociation, generic)
r_s1 = 50.0e-9 # m (dry radius)
r_s2 = 100.0e-9 # m (dry radius)

B1 = dsc*M_w*rho_s1/(M_s1*rho_w) # Kohler hygroscopicity coefficient
print("B1 =",B1)

#data for (NH4)2(SO4) (2)
M_s2 = 132.134e-3 # Kg (molecular weight) [PubChem]
rho_s2 = 1.77e3 # Kg/m^3 (density at ~4 deg C)
dsc = 2 # (degree of ionic dissociation, generic)

B2 = dsc*M_w*rho_s2/(M_s2*rho_w) # Kohler hygroscopicity coefficient
print("B2 =",B2)
import numpy as np
import matplotlib.pyplot as plt
r = []
cycles = 2
r_base = 1.0e-7
for i in range(0,cycles):
    r = r + np.linspace(r_base*10**i,r_base*10**(i+1),1001).tolist()
S11 = [100*(A/rad-B1*r_s1**3/rad**3) for rad in r]
S12 = [100*(A/rad-B1*r_s2**3/rad**3) for rad in r]
S21 = [100*(A/rad-B2*r_s1**3/rad**3) for rad in r]
S22 = [100*(A/rad-B2*r_s2**3/rad**3) for rad in r]
#print(r)
"""
S_c1 = 2*(A/(3*r_s1))**1.5/B**0.5
r_c1 = (3*B*r_s1**3/A)**0.5
print(S_c1,r_c1)
S_c2 = 2*(A/(3*r_s2))**1.5/B**0.5
r_c2 = (3*B*r_s2**3/A)**0.5
print(S_c2,r_c2)
"""
plt.grid()
plt.ylim(-0.2,0.3)
plt.axhline(y=0, color='k', linewidth=0.5)
plt.axvline(x=0, color='k', linewidth=0.5)
plt.semilogx(r,S11,label = 'NaCl 50nm')
plt.semilogx(r,S12,label = 'NaCl 100nm')
plt.semilogx(r,S21,linestyle = 'dashed',label = '(NH4)2(SO4) 50nm')
plt.semilogx(r,S22,linestyle = 'dashed',label = '(NH4)2(SO4) 100nm')
plt.legend()
