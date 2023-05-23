#import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import isacalc
import scipy.stats
import scipy.integrate
import scipy.optimize
#import blimb_calc
from ambiance import Atmosphere


M = {'H2':2.015894, 'He':4.0026022}
glider_data = np.genfromtxt('sailplane_data.csv',delimiter=',')[1:,1:].T


def glider_statistics(val, output, inp="OEM"):
		names = ["OEM","MTOM", "Mpl","Mwing","Mfus","Mhstab","b","AR","S","S/W","cbar"]
		inp_indx = names.index(inp)
		output_indx = names.index(output)
		ratio = np.average(glider_data[output_indx])/np.average(glider_data[inp_indx])
		return ratio*val

def balloon_mass(payload):
	return 0.316*payload + 0.261

#this is a slightly different estimation
def traditional_balloon_weight(volume):
    weight = (4.25*volume + 340) * 0.001
    return weight

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
    volume = total_mass / delta_rho[0]
    return volume

def gas_mass_balloon(m_launch, h_max, ascent_rate, Molw_gas):
	m_gas = calculate_required_size(m_launch, h_max, Molw_gas)*gas_density(Molw_gas, h_max)
	m_extra_gas = 0.01
	alt, vel = climb_rate(m_extra_gas, m_gas, Molw_gas)
	v_avg = np.average(vel)
	v_avg_log = [v_avg]
	while v_avg <= ascent_rate:
		m_extra_gas = m_extra_gas + 0.1
		alt, vel = climb_rate(m_extra_gas, m_gas, Molw_gas)
		v_avg = np.average(vel)
		v_avg_log.append(v_avg)
	return m_gas, m_extra_gas, v_avg

def calculate_required_size(m_tot, altitude, molar_mass):
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    volume = m_tot / delta_rho
    return volume

def balloon_mass_update(m_tot, h_max, molar_mass):
	#yeah i dont really know what this function does
    volume = calculate_required_size(m_tot, h_max, molar_mass)
    mass = (4.25*volume + 340) * 0.001
    return mass

def climb_rate(m_gas_extra, m_gas_basic, molar_mass):
    alt = np.arange(0, 33000, 100.)
    vel = np.zeros_like(alt)
    for i in range(np.size(alt)):
        rho_air = Atmosphere(float(alt[i])).density[0]
        vol = (m_gas_basic + m_gas_extra) / (rho_air - gas_density(molar_mass, alt[i]))
        rad = ((3*vol)/(4*np.pi))**(1/3)
        area = np.pi*(rad**2)
        v = ((2*9.81*m_gas_extra)/(rho_air*0.5*area))**0.5
        vel[i] = v
    return alt, vel

def glide_ratio(CD0, e, AR):
	CL = (np.pi * e * AR * CD0) ** 0.5
	CD = CD0 + (CL ** 2) / (np.pi * AR * e)
	return CL / CD, CL

def glide_flight_props(h_start, h_end, WoS, gamma, CL):
	glide_range = (h_start - h_end)*gamma
	flight_time = scipy.integrate.quad(lambda h: gamma/(v_from_S(WoS, isacalc.isa(h)['dens'], CL)), h_end, h_start)[0]
	return glide_range, flight_time

def vel_profile(h_start, WoS, CL):
	alt = np.arange(0, h_start, 100.)
	v = np.zeros_like(alt)
	M = np.zeros_like(alt)
	for i in range(np.size(alt)):
		rho_air = Atmosphere(float(alt[i])).density[0]
		a = Atmosphere(float(alt[i])).speed_of_sound[0]
		v[i] = v_from_S(WoS,rho_air, CL)
		M[i] = v[i]/a

	return np.max(M), np.max(v), np.min(v)

v_from_S = lambda WoS, rho, CL: np.sqrt((WoS)*(2/rho)*(1/CL))

def third_estimation(payload, heights, aero_params, Molw_gas, S):
	h_max, h_level_flight, h_buffer = heights
	CD0, e, AR = aero_params

	required_ascent_rate = 5 #m/s
	#balloon mass:
	
	m_glide = payload/0.45 #kg
	m_balloon_assumed = 3 #kg
	m_balloon_iteration = m_balloon_assumed
	m_balloon_log = []
	for i in range(10):
		m_tot_initial = m_glide + m_balloon_iteration
		m_balloon_iteration = balloon_mass_update(m_tot_initial, h_max, Molw_gas)
		m_balloon_log.append(m_balloon_iteration)
	m_balloon_final = m_balloon_iteration
	m_launch = m_glide+m_balloon_final

	W = m_glide * 9.81

	#extra gas estimation:
	m_gas, m_extra_gas, v_avg = gas_mass_balloon(m_launch, h_max, required_ascent_rate, Molw_gas)


	#lift over drag
	gamma, CLopt = glide_ratio(CD0, e, AR)

	#velocity calculation
	M_max, v_max, v_min = vel_profile(h_level_flight, W/S, CLopt) 

	glide_range, flight_time = glide_flight_props(h_level_flight, h_buffer, W/S, gamma, CLopt)

	print("Third estimation")
	print("m_balloon: {} [kg], m_glider: {} [kg], m_launch: {} [kg], m_gas: {} [kg] , extra: {} [kg]".format(m_balloon_final, m_glide, m_launch, m_gas, m_extra_gas))
	print("ascent rate: {} [m/s], Max mach: {} [-], max velocity: {} [m/s], min_velocity: {} [m/s]".format(v_avg, M_max, v_max, v_min))
	print("C_L: {} [-], Glide ratio: {} [-], Range: {} [km], Flight time: {} [h]".format(CLopt, gamma,glide_range/1000,flight_time/3600))

	return glide_range, flight_time, m_launch, m_glider, m_balloon


def first_estimation(payload, heights, aero_params, Molw_gas, carry_balloon):

	h_max, h_level_flight, h_buffer = heights
	CD0, e, AR = aero_params

	MTOm = 10 #kg, initial mass guess	
	required_ascent_rate = 5 # m/s

	#Calculate the mass based on ratios from competition gliders, this is iterative
	for i in range(100):
		m_balloon = balloon_mass(MTOm)
		OEM = glider_statistics(m_pl+m_balloon*carry_balloon, "OEM", "Mpl")
		#OEM = (m_pl+m_balloon*carry_balloon) * 0.25
		MTOm = m_pl+m_balloon + OEM
	S = glider_statistics(MTOm, "S", "MTOM")
	#weight for velocity estimation and stuff
	W = MTOm*9.81

	#balloon volume and gass weight required
	V_balloon = calculate_required_volume(MTOm, h_max, Molw_gas)
	m_gas, m_extra_gas, v_avg = gas_mass_balloon(MTOm, h_max, required_ascent_rate, Molw_gas)

	#estimate glide ratio based on aerodynamic parameters
	gamma, CLopt = glide_ratio(CD0, e, AR)
	
	
	rho_high =  isacalc.isa(h_max)['dens']
	rho_level_flight = isacalc.isa(h_level_flight)['dens']
	rho_sl = isacalc.isa(0)['dens']

	v = v_from_S(W/S, rho_level_flight, CLopt)
	b = np.sqrt(S*AR)
	glide_range, flight_time = glide_flight_props(h_level_flight, h_buffer, W/S, gamma, CLopt)

	print("Fisrt estimation")
	print("MTOM: {}, OEM: {}, balloon:{}, gas amt: {} [kg]".format(MTOm, OEM, m_balloon, m_gas))
	print("Wing surface: {} [m2], Wing span: {} [m], speed at top altitude: {} [m/s]".format(S, b, v))
	print("C_L: {} [-], Glide ratio: {} [-], Range: {} [km], Flight time: {} [h]".format(CLopt, gamma,glide_range/1000,flight_time/3600))

	return glide_range, flight_time, MTOm, OEM, m_balloon


if __name__ == '__main__':
	

	m_pl = 2
	heights_1 = [33000, 25000, 500]
	aero_params_1 = [0.02, 0.80, 15]
	

	first_estimation(m_pl, heights_1, aero_params_1, M["H2"], True)


	m_electronic = 1.119
	aero_params_2 = [0.02, 0.85, 5]
	S = 0.8

	third_estimation(m_electronic, heights_1, aero_params_2, M["H2"], S)

	







