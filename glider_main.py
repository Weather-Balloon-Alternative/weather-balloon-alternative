import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import isacalc
import scipy.stats


glider_data = np.genfromtxt('sailplane_data.csv',delimiter=',')[1:,2:].T/50
def glider_statistics(val, output, inp="OEW"):
		
		names = ["OEW","MTOW","Mwing","Mfus","Mhstab","b","AR","S","S/W","cbar"]
		inp_indx = names.index(inp)
		output_indx = names.index(output)
		print(glider_data)		
		linfit = scipy.stats.linregress(glider_data[inp_indx], glider_data[output_indx])
		print(linfit.slope, linfit.intercept, linfit.rvalue)
		OEW_range = np.arange(0, 10, 0.1)

		#plt.plot(glider_data[inp_indx], glider_data[output_indx], marker='o', linestyle="None")
		#plt.plot(OEW_range, OEW_range*linfit.slope+linfit.intercept)
		#plt.show()
		return linfit.intercept+linfit.slope*val

print(glider_statistics(3, "MTOW"))


W = 5 #kg
h_max = 30000

CD0 = 0.02
e = 0.85
A = 15


CL = (np.pi*e*A*CD0)*0.5
CD = CD0 + (CL**2)/(np.pi*A*e)
print(CL/CD)

rho_high =  isacalc.isa(h_max)['dens']
rho_sl = isacalc.isa(0)['dens']

S = lambda v, rho: W/(0.5 * rho * (v**2) * CL)
v_from_S = lambda S, rho: np.sqrt((W/S)*(2/rho)*(1/CL))


#vv = np.arange(10, 100, 0.1)
SS = np.arange(0.3, 3, 0.01)

#plt.plot(vv, S(vv, rho_high), vv, S(vv, rho_sl))
plt.plot(SS, v_from_S(SS, rho_high), SS, v_from_S(SS, rho_sl))
plt.ylabel("velocity")
plt.xlabel("Wing surface area")
plt.show()





