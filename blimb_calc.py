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

def estimate_area(volume):
    return ((4*np.pi) * ((3*volume)/(4*np.pi))**(2/3)) 

def estimate_weigth(volume, material_density, material_thickness):
    return estimate_area(volume)* material_thickness * material_density


max_weight = 10
payload = 2
alt = np.arange(0, 33000, 100)
M = {'H2':2.015894, 'He':4.0026022}

vols_H2 = calculate_required_size(max_weight, alt, M['H2']) #crs(alt)
vols_He = calculate_required_size(max_weight, alt, M['He']) #crs(alt)

mass_H2 = estimate_weigth(vols_H2, 75, 1 * 0.001)


plt.plot(alt, mass_H2)
plt.plot(alt, np.ones(330)*(max_weight-payload))
#plt.plot(alt, vols_He)
plt.show()