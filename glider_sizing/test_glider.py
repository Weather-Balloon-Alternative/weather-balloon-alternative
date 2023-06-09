#import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import glider_sizing

import numpy as np

def test_balloon_mass():
	assert 1 < glider_sizing.balloon_mass(1.361+1.2) < 1.3

def test_gas_density():
	#test if density of air at SL is in correct range
	M_air = 28.9647 #avg molar mass of air
	rho_air_sl = 1.225 #avg sl density
	alt = 0
	rho_calc = glider_sizing.gas_density(M_air, alt)
	assert np.isclose(rho_calc, rho_air_sl)

def test_calc_volume():
	m = 1 #kg
	M_gas = 2
	alt = 0
	V_calc = glider_sizing.calculate_required_volume(m, alt, M_gas)
	V_req = 0.8769 # from manual calculation
	assert np.isclose(V_calc, V_req, rtol=1e-3) 

def test_climb_rate():
	#I dont really know what do 
	assert True

def test_glide_ratio_CL():
	e = 1/np.pi
	AR = 10
	CD0 = 0.036
	
	CL_req = 0.6 #based on manual calc
	CL_calc = glider_sizing.glide_ratio(CD0, e, AR)[1]
	assert np.isclose(CL_calc, CL_req)

def test_glide_ratio_gamma():
	e = 1/np.pi
	AR = 10
	CD0 = 0.036
	
	gamma_req = 0.6/0.072 #based on manual calc
	gamma_calc = glider_sizing.glide_ratio(CD0, e, AR)[0]
	assert np.isclose(gamma_calc, gamma_req)

def test_glide_flight_props():
	h_start = 1000
	h_end = 0
	WoS = 100 # wing loading
	gamma = 10
	CL = 0.5 

	props_req = (10000, 200) #range is pretty obvs, flight time is not tested rn
	props_calc = glider_sizing.glide_flight_props(h_start, h_end, WoS, gamma, CL)
	print(props_calc[1])
	assert np.isclose(props_calc[0], props_req[0])

def test_v_S_val():
	WoS = 20
	rho = 1
	CL = 0.7
	v_req = 7.55928946 #manual calculation
	v_calc = glider_sizing.v_from_S(WoS, rho, CL)
	assert np.isclose(v_calc, v_req)
	WoS *= 4 #check that the relation is indeed square root
	v_req *= 2
	v_calc = glider_sizing.v_from_S(WoS, rho, CL)
	assert np.isclose(v_calc, v_req)

def test_vel_profile():
	h_start = 10000
	WoS = 20
	CL = 0.7
	rho_min = 1.225
	M_max, v_max, v_min = glider_sizing.vel_profile(h_start, WoS, CL)
	assert v_min < v_max
	assert np.isclose(v_min, glider_sizing.v_from_S(WoS, rho_min, CL)) #v_from_S is already verified