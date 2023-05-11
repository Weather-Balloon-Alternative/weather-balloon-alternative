import numpy as np
from ambiance import Atmosphere
from matplotlib import pyplot as plt

M = {'H2':2.015894, 'He':4.0026022}

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
    atmos = Atmosphere(float(altitude))
    density = (atmos.pressure * molar_mass * 0.001) / (R * (atmos.temperature + temperature_offset))
    return density[0]

def climb_rate(m_gas_extra, m_gas_basic, molar_mass):
    alt = np.arange(0, 35000, 10.)
    vel = np.zeros_like(alt)
    for i in range(np.size(alt)):
        rho_air = Atmosphere(float(alt[1])).density[0]
        vol = (m_gas_basic + m_gas_extra) / (rho_air - gas_density(molar_mass, alt[i]))
        rad = ((3*vol)/(4*np.pi))**(1/3)
        area = np.pi*(rad**2)
        v = ((2*9.81*m_gas_extra)/(rho_air*0.5*area))**0.5
        vel[i] = v 
    return alt, vel

alt, vel = climb_rate(1.1, 0.9, M['H2'])

plt.plot(alt, vel)
plt.show()
