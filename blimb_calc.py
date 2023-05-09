from ambiance import Atmosphere
import numpy as np
from matplotlib import pyplot as plt

def gas_density(molar_mass, altitude):
    R = 8.3144598
    atmos = Atmosphere(altitude)
    rho = (atmos.pressure * molar_mass * 0.001) / (R * atmos.temperature)
    return rho

def calculate_required_size(total_mass, altitude, molar_mass):
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    volume = total_mass / delta_rho
    return volume

def estimate_area(volume):
    r = ((3*volume)/(4*np.pi))**(1/3)
    return 4*np.pi*(r**2)

def estimate_weigth(volume, material_density, material_thickness):
    return estimate_area(volume)* material_thickness * material_density

def balloon_weight(total_mass, altitude, molar_mass):
    vol = calculate_required_size(total_mass, altitude, molar_mass)
    mass_balloon = estimate_area(vol)*1390*((500*((3*vol)/(4*np.pi))**(1/3))/(200*1000000)) 
    return mass_balloon

max_weight = 50
payload = 2
alt = np.arange(0, 33000, 100)
M = {'H2':2.015894, 'He':4.0026022}
mass = np.zeros(np.size(alt))#mylar_weight(max_weight, alt, M['H2'], 2)





# for i in range(np.size(alt)):
#     m = max_weight
#     for j in range(50):
#         m = mylar_weight(m, alt[i], M['H2'], payload)
#     mass[i] = m

# plt.plot(alt, mass)
# plt.show()


h = 20000
m_guess = np.arange(1,100, 0.1)
m_calc = []
for i in range(np.size(m_guess)):
    m_calc.append(balloon_weight(m_guess[i]+2, h, M['H2']))

plt.plot(m_guess, m_calc)
plt.show()


#vols_H2 = calculate_required_size(max_weight, alt, M['H2']) #crs(alt)
#vols_He = calculate_required_size(max_weight, alt, M['He']) #crs(alt)

#mass_H2 = estimate_weigth(vols_H2, 75, 1 * 0.001)


#print(estimate_area(calculate_required_size(10, 33000, 2)))
#plt.plot(alt, mass_H2)
#plt.plot(alt, np.ones(330)*(max_weight-payload))
#plt.plot(alt, vols_He)
#plt.show()