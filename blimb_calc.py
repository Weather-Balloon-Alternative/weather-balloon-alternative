from ambiance import Atmosphere
import numpy as np
from matplotlib import pyplot as plt

def gas_density(molar_mass, altitude, temperature_offset=0):
    """
        Returns the density[kg/m^3] of a gas with given molar mass[g/mol] at a given altitude[m]
    """
    R = 8.3144598
    atmos = Atmosphere(altitude)
    rho = (atmos.pressure * molar_mass * 0.001) / (R * (atmos.temperature + temperature_offset))
    return rho

def calculate_required_volume(total_mass, altitude, molar_mass):
    """
        Returns the required volume[kg^3] to lift a given mass[kg] to a given altitude[m] using a gas with a given molar mass [g/mol]
    """
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    volume = total_mass / delta_rho
    return volume

def area_of_sphere(volume):
    """
        Returns the surface area of a sphere given the volume
    """
    r = ((3*volume)/(4*np.pi))**(1/3)
    return 4*np.pi*(r**2)

def calculate_weight(volume, material_density, material_thickness):
    """
        Calculates the weight[kg] of a sphereical balloon given the volume[m^3] and material density[kg/m^3] and thickness[m]
    """
    return area_of_sphere(volume)* material_thickness * material_density

def balloon_weight_mylar(total_mass, altitude, molar_mass, custom_vol = 0):
    """
        Calculates the mass of the balloon using mylar. a pressure difference of 500Pa, and a safety factor of 2
    """
    pressure_difference = 500
    safety_factor  = 2
    mylar_strength = 100 * (10**6)
    mylar_density = 1390

    vol = calculate_required_volume(total_mass, altitude, molar_mass)
    if custom_vol != 0: 
        vol = custom_vol
    t = ((pressure_difference * ((3*vol)/(4*np.pi))**(1/3))/(safety_factor * mylar_strength))
    mass_balloon = area_of_sphere(vol) * mylar_density * t
    return mass_balloon

def availible_payload(volume, altitude, molar_mass, use_mylar=True):
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    lifting_ability = delta_rho * volume
    weight_balloon = 0
    if use_mylar:
        weight_balloon = balloon_weight_mylar(0, 0, 0, custom_vol=volume)
    else:
        weight_balloon = traditional_balloon_weight(volume)
    return (lifting_ability - weight_balloon)

def traditional_balloon_weight(volume):
    weight = (4.25*volume + 340) * 0.001
    return weight

max_weight = 11
payload = 2
alt = np.arange(0, 33000, 100)
M = {'H2':2.015894, 'He':4.0026022}
mass = np.zeros(np.size(alt))#mylar_weight(max_weight, alt, M['H2'], 2)


alt = 33000
vol = np.arange(10, 3000, 10)
pl = np.zeros(np.shape(vol))

for i in range(np.size(vol)):
    pl[i] = availible_payload(vol[i], alt, M['H2'], use_mylar=False)

plt.plot(vol, 9.65*np.ones(np.size(vol)))
plt.plot(vol, pl)
plt.show()

