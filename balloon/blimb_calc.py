from ambiance import Atmosphere
import numpy as np
from matplotlib import pyplot as plt

def gas_density(molar_mass, altitude, temperature_offset=0):
    '''
         Returns the density of a gas with given molar mass at a given altitude with an optional temperature offset

            Parameters:
                    molar_mass (float): Molar mass of the gas in [g/mol]
                    altitude (float): Altitude in [m]
                    temperature_offset (float): Termperatre offset in degree

            Returns:
                    density (float): density of gas in kg/m^3
    '''
    R = 8.3144598
    atmos = Atmosphere(altitude)
    density = (atmos.pressure * molar_mass * 0.001) / (R * (atmos.temperature + temperature_offset))
    return density

def calculate_required_volume(total_mass, altitude, molar_mass):
    '''
         Returns the required volume[kg^3] to lift a given mass[kg] to a given altitude[m] using a gas with a given molar mass [g/mol]

            Parameters:
                    total_mass (float): total balloon mass + payload [kg]
                    altitude (float): Altitude in [m]
                    molar_mass (float): Molar mass of the gas in [g/mol]

            Returns:
                    volume (float): volume of the balloon in [m^3]
    '''
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

def calculate_mass(volume, material_density, material_thickness):
    """
        Calculates the weight[kg] of a sphereical balloon given the volume[m^3] and material density[kg/m^3] and thickness[m]
    """
    return area_of_sphere(volume)* material_thickness * material_density

def balloon_mass_mylar(total_mass, altitude, molar_mass, custom_vol = 0):
    """
        Calculates the mass of the balloon using mylar. a pressure difference of 500Pa, and a safety factor of 2
    """ 
    pressure_difference = 150
    safety_factor  = 0.75
    mylar_strength = 100 * (10**6)
    mylar_density = 1390

    vol = calculate_required_volume(total_mass, altitude, molar_mass)
    if custom_vol != 0: 
        vol = custom_vol
    t = ((pressure_difference * ((3*vol)/(4*np.pi))**(1/3))/(safety_factor * mylar_strength))
    mass_balloon = area_of_sphere(vol) * mylar_density * t
    return mass_balloon

def availible_payload(volume, altitude, molar_mass, use_mylar=True):
    """
        Calculate the availible payload [kg] for a given baloon volume[m^3] at a given altitude[m] using either mylar or latex balloons
    """
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    lifting_ability = delta_rho * volume
    mass_balloon = 0
    if use_mylar:
        mass_balloon = balloon_mass_mylar(0, 0, 0, custom_vol=volume)
    else:
        mass_balloon = traditional_balloon_mass(volume)
    return (lifting_ability - mass_balloon)

def traditional_balloon_mass(volume):
    """
        Calculates the mass[kg] of a conventional latex balloon for a given burst volume[m^3]
    """
    mass = (4.25*volume + 340) * 0.001
    return mass

def drag(v):
    c_d = 0.5
    rho = 1.225
    S = 620
    D = 0.5 * rho * (v**2) * S * c_d
    return D

def lift_at_sl(altitude, volume, molar_mass):
    rho_gas_alt = gas_density(molar_mass, altitude)
    rho_gas_sl = gas_density(molar_mass, 0)
    rho_air_sl = Atmosphere(0).density
    L = (rho_gas_alt / rho_gas_sl) * (rho_air_sl - rho_gas_sl) * 9.81 * volume
    return L


max_weight = 11
payload = 2
alt = np.arange(0, 33000, 100)
M = {'H2':2.015894, 'He':4.0026022}
mass = np.zeros(np.size(alt))#mylar_weight(max_weight, alt, M['H2'], 2)




vol = np.arange(10, 3000, 10)
pl = np.zeros(np.shape(vol))
pl2 = np.zeros(np.shape(vol))

rho_air = Atmosphere(alt).density
rho_gas = gas_density(M['H2'], alt)

plt.plot(rho_air, rho_gas)
plt.show()
drho  = rho_air - rho_gas
plt.plot(alt, drho)
plt.show()

m_pl = 5
m_bl = 3.825
m_t = m_pl+m_bl
Vbl = 820

m_liftgas = 820 * gas_density(M['H2'], 33000)
vatalt = m_liftgas / gas_density(M['H2'], alt)
lift = drho * vatalt

plt.plot(alt, lift)
plt.show()


alt = 33000
# lift = lift_at_sl(alt, vol, M['H2'])
# plt.plot(vol, lift)
# plt.show()

for i in range(np.size(vol)):
     pl[i] = availible_payload(vol[i], alt, M['He'], use_mylar=False)
     pl2[i] = availible_payload(vol[i], alt, M['H2'], use_mylar=False)


# plt.plot(vol, ((lift/9.81) / pl2))
# plt.show()
# #plt.plot(vol, 9.65*np.ones(np.size(vol)))
plt.plot(vol, pl, label='Helium')
plt.plot(vol, pl2, label='Hydrogen')
plt.plot(vol, 5*np.ones(np.size(vol)), label='4kg')
plt.legend()
plt.show()

# vs = np.arange(0, 31, 1)
# D = drag(vs)

# plt.plot(vs, D)
# #plt.plot(vs, np.)
# plt.show()


# #m_H2 = 
