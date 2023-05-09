from ambiance import Atmosphere
import numpy as np
from matplotlib import pyplot as plt

def gas_density(molar_mass, altitude):
    R = 8.3144598
    atmos = Atmosphere(altitude)
    rho = (atmos.pressure * molar_mass * 0.001) / (R * atmos.temperature)
    return rho

def calculate_required_size(payload, altitude, molar_mass):
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    volume = payload / delta_rho
    return volume
M = {'H2':2.015894, 'He':4.0026022}


m_module = 5
#desired control altitude, low but above planes
alt = 500
atmos = Atmosphere(alt)
vol = calculate_required_size(m_module, alt, M['H2'])
r = ((3*vol)/(4*np.pi))**(1/3)
A_front = np.pi*(r**2)
CD = 0.47 #sphere
rho = atmos.density

#what velocity do you want

#for 3 hours it is 20m/s if you drift of 200km
V=20

D = 0.5*CD*rho*(V**2)*A_front
Thrust = D/9.81
print(D*V*600)
print("thrust in kg:", Thrust)