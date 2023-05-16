import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import isacalc
import scipy.stats
import scipy.integrate
import scipy.optimize
import blimb_calc


glider_data = np.genfromtxt('sailplane_data.csv',delimiter=',')[1:,1:].T
def glider_statistics(val, output, inp="OEM"):
		
		names = ["OEM","MTOM", "Mpl","Mwing","Mfus","Mhstab","b","AR","S","S/W","cbar"]
		inp_indx = names.index(inp)
		output_indx = names.index(output)
		#print(glider_data)		
		#linfit = scipy.stats.linregress(glider_data[inp_indx], glider_data[output_indx])
		#print(linfit.slope, linfit.intercept, linfit.rvalue)
		#OEW_range = np.arange(0, 10, 0.1)

		#plt.plot(glider_data[inp_indx], glider_data[output_indx], marker='o', linestyle="None")
		#plt.plot(OEW_range, OEW_range*linfit.slope+linfit.intercept)
		#plt.show()
		#return linfit.intercept+linfit.slope*val
		ratio = np.average(glider_data[output_indx])/np.average(glider_data[inp_indx])
		return ratio*val

def balloon_weight(payload):
	return 0.316*payload + 0.261


def third_estimation(payload):
	h_max = 33000 #maximum height to be reached by balloon
	h_level_flight = 25000 #height at wich steady gliding flight begins
	h_buffer = 500 #the glider should reach the target with this height remaining, to account for extra go around, inneficiencies etc 

	CD0 = 0.02
	e = 0.85
	AR = 4.5
	S= 10



	M = {'H2':2.015894, 'He':4.0026022} 
	MTOm = 10 #kg


	for i in range(100):
		m_balloon = balloon_weight(MTOm)

	

	W = MTOm*9.81

	v_from_S = lambda S, rho: np.sqrt((W/S)*(2/rho)*(1/CL))

	CL = (np.pi*e*AR*CD0)**0.5
	CD = CD0 + (CL**2)/(np.pi*AR*e)
	glide_ratio = CL/CD
	glide_range = (h_level_flight - h_buffer)*glide_ratio
	flight_time = scipy.integrate.quad(lambda h: glide_ratio/(v_from_S(S, isacalc.isa(h)['dens'])), h_buffer, h_level_flight)[0]


	#print("MTOM: {}, OEM: {}, balloon:{}, gas amt: {} [kg]".format(MTOm, OEM, m_balloon, 1))
	print("C_L: {} [-], Glide ratio: {} [-], Range: {} [km], Flight time: {} [h]".format(CL, glide_ratio,glide_range/1000,flight_time/3600))

def second_estimation(payload):
	h_max = 33000 #maximum height to be reached by balloon
	h_level_flight = 25000 #height at wich steady gliding flight begins
	h_buffer = 500 #the glider should reach the target with this height remaining, to account for extra go around, inneficiencies etc 

	CD0 = 0.02
	e = 0.85
	AR = 4.5
	wing_loading = 30
	wing_frac_ar = 0.035
	fus_total_frac = 0.2


	M = {'H2':2.015894, 'He':4.0026022} 
	MTOm = 10 #kg

	for i in range(100):
		m_balloon = balloon_weight(MTOm)
		m_fus = MTOm*fus_total_frac
		m_wing = wing_frac_ar*AR*MTOm
		MTOm = m_balloon + m_fus + m_wing + m_pl
		OEM = m_fus+m_wing
	

	W = MTOm*9.81
	print("MTOM: {}, OEM: {}, balloon:{}, gas amt: {} [kg]".format(MTOm, OEM, m_balloon, 1))
	

def CD0AR_plot(glide_range, h_start, h_end):
	e = 0.8
	drag_ratio = 1
	glide_ratio = glide_range/(h_start-h_end)
	CLopt = lambda CD0, AR: ((np.pi*e*AR*CD0)**0.5)
	f = lambda CD0, AR: CLopt(CD0, AR)/((CD0 + (CLopt(CD0, AR)**2)/(np.pi*AR*e))*drag_ratio) - glide_ratio
	ARf = lambda CD0: scipy.optimize.fsolve(lambda x: f(CD0, x), 5)[0]

	print(ARf(0.02))
	CD0s = np.arange(0.01, 0.04, 0.001)
	ARs = []
	for i in CD0s:
		ARs.append(ARf(i))

	plt.plot(CD0s, ARs)
	plt.show()


def first_estimation(payload):
	h_max = 33000 #maximum height to be reached by balloon
	h_level_flight = 25000 #height at wich steady gliding flight begins
	h_buffer = 500 #the glider should reach the target with this height remaining, to account for extra go around, inneficiencies etc 

	CD0 = 0.02
	e = 0.85
	AR = 4.5

	M = {'H2':2.015894, 'He':4.0026022} 
	MTOm = 10 #kg
	
	print(glider_statistics(1, "OEM", "Mpl"))
	carry_balloon = True

	for i in range(100):
		#V_balloon = blimb_calc.calculate_required_volume(MTOm, h_max, M["H2"])
		#m_balloon = blimb_calc.traditional_balloon_weight(V_balloon)
		m_balloon = balloon_weight(MTOm)
		OEM = glider_statistics(m_pl+m_balloon*carry_balloon, "OEM", "Mpl")
		#OEM = (m_pl+m_balloon*carry_balloon) * 0.25
		MTOm = m_pl+m_balloon + OEM
	W = MTOm*9.81
	print(glider_statistics(OEM, "Mwing", "OEM"), "mwing")

	V_balloon = blimb_calc.calculate_required_volume(MTOm, h_max, M["H2"])
	m_gas = V_balloon*blimb_calc.gas_density(M["H2"], h_max)

	CL = (np.pi*e*AR*CD0)**0.5
	CD = CD0 + (CL**2)/(np.pi*AR*e)
	
	rho_high =  isacalc.isa(h_max)['dens']
	rho_level_flight = isacalc.isa(h_level_flight)['dens']
	rho_sl = isacalc.isa(0)['dens']

	S_from_v = lambda v, rho: W/(0.5 * rho * (v**2) * CL)
	v_from_S = lambda S, rho: np.sqrt((W/S)*(2/rho)*(1/CL))

	S = glider_statistics(MTOm, "S", "MTOM")
	v = v_from_S(S, rho_level_flight)
	#v = S_from_v(100, rho_high)
	glide_ratio = CL/CD

	b = np.sqrt(S*AR)
	glide_range = (h_level_flight - h_buffer)*glide_ratio

	flight_time = scipy.integrate.quad(lambda h: glide_ratio/(v_from_S(S, isacalc.isa(h)['dens'])), h_buffer, h_level_flight)[0]


	print("MTOM: {}, OEM: {}, balloon:{}, gas amt: {} [kg]".format(MTOm, OEM, m_balloon, m_gas))
	print("Wing surface: {} [m2], Wing span: {} [m], speed at top altitude: {} [m/s]".format(S, b, v))
	print("C_L: {} [-], Glide ratio: {} [-], Range: {} [km], Flight time: {} [h]".format(CL, glide_ratio,glide_range/1000,flight_time/3600))


print(glider_statistics(2, "OEM", "Mpl"))

if __name__ == '__main__':
	m_pl = 2
	#first_estimation(m_pl)
	#CD0AR_plot(300, 25, 0.5)
	#second_estimation(m_pl)
	third_estimation(m_pl)
	






