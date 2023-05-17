import numpy as np 
import math
import pandas
import ambiance
import matplotlib.pyplot as plt

def climb_characteristics(weights, surface_areas, rho_array, ROC, V_0, C_d_0, AR, e, plots, SAF_energy_density, th_efficiency):
    '''
    args:   
            weights: list with floats,              gives all different weights for which will be analyzed
            suface_areas: list with floats,         gives the different surface areas which are possible to be used. 
            rho_array: array with floats,           gives density at the different altitudes that are selected for the climb
            ROC: float,                             rate of climb that is selected 
            V_0: float,                             minimum speed at which the system is allowed to fly
            C_d_0: float,                           zero lift drag
            AR: float,                              aspect ratio of the aircraft
            e: float between 0 and 1,               oswald efficiency factor
            plots: Bool,                            plotting results if required
            SAF_energy_density: float,              Energy in a liter of fuel
            th_efficiency: float between 0 and 1,   thermal efficiency
    outputs:
            S_return: list of floats,               optimum surface area for the design requiring minimum energy
            l_SAF_return: list of floats,           estimate of kilograms of Sustainable Aviation Fuel (SAF) required
            P_max_return: list of flaots,           maximum power required during ascent
            min_E_tot_loc: int,                     location of the minimum energy required in array
    '''
    #making return lists
    S_return =      []
    l_saf_return =  []
    P_max_return =  []

    #making lists for measurements
    E_tot_array = []
    S_array = []

    #itterating over the inputted weights
    for W in weights:
        #itterating over the different surface areas
        for S in surface_areas:
            C_L = 2*W/(rho_array*S*V_0**2)

            #appending CL in case 
            for j in range(len(C_L)):
                if C_L[j]>C_L_climb:
                    C_L[j]=C_L_climb

            #calculating velocity and expected drag        
            V = np.sqrt(2*W/(rho_array*S*C_L))
            C_D = C_d_0 + C_L**2/(np.pi*AR*e)

            #Power consumption calculations
            P_equilibrium=0.5*rho_array*S*V**3*C_D/prop_efficiency
            P_climb = np.ones(len(C_L))*ROC*W/prop_efficiency
            P_tot = P_equilibrium+P_climb
            
            #creating time array
            T = np.ones(len(C_L))*1/ROC

            #appending energy required to array
            E_tot = P_tot * T
            E_tot_array.append(sum(E_tot))
            S_array.append(S)
        
        #plotting results if required
        if plots:
            plt.plot(S_array,E_tot_array, label=str("Weight="+str(W)+"[N]"))

        #determining at which surface area minimum energy is required
        min_E_tot_loc = int(np.argwhere(E_tot_array==min(E_tot_array)))
        S = S_array[min_E_tot_loc]
        
        #calculating kilograms of saf required
        l_SAF = convert_fuel_volume(E_tot_array[min_E_tot_loc],thermal_efficiency=th_efficiency, energy_density=SAF_energy_density)
        
        #appending data to full return lists
        S_return.append(S)
        l_saf_return.append(l_SAF)
        P_max_return.append(max(P_tot))

        #print statements to summarize the data
        print("minimum energy required : ",E_tot_array[min_E_tot_loc])
        print("Weight: ", W)
        print("optimum surface area: ",S_array[min_E_tot_loc])
        print("kilograms SAF required for ascent: ", l_SAF)
        print("maximum power requried: ", max(P_tot))
        print("")
        
        #clearing the lists
        S_array = []
        E_tot_array =[]
    if plots:
        plt.legend()
        plt.xlabel("S[m^2]")
        plt.ylabel("Energy[J]")
        plt.show()
    return S_return, l_saf_return, P_max_return, min_E_tot_loc

def loitering_fuel(loitering_time, loitering_altitude, C_L_loiter, S, AR, e, SAF_energy_density, th_efficiency):
    """
    Args:
            loitering_time: array of length n with floats,      gives the loitering times required at each altitude.
            loitering_altitude: array of length n with floats,  gives the loitering altitudes at each altitude.
            C_L_loiter: array of length n with floats,          gives the loitering expected lift coefficients
            S: float,                                           surface area of craft
            AR: float,                                          aspect ratio of the aircraft
            e: float between 0 and 1,                           oswald efficiency factor of the aircraft
            SAF_energy_density: float,                          Energy in a liter of fuel
            th_efficiency: float between 0 and 1,               thermal efficiency
    outputs:
            l_SAF_loit: float,                                  estimate of kilograms of Sustainable Aviation Fuel (SAF) required
    """
    E_loit = 0
    for i in range(len(loitering_time)):
        #determine altitude of loiter phase
        alt = loitering_altitude[i]

        #determine properties for power calculations
        rho = ambiance.Atmosphere(alt).density
        V = np.sqrt((2*W)/(rho*C_L_loiter[i]*S))
        C_D = C_d_0 + C_L_loiter[i]**2/(np.pi*AR*e)
        
        #calculate required fuel by power calculation
        P_loit = 1/2*V**3*rho*S*C_D/prop_efficiency
        print("for phase", i, " fuel required: ", convert_fuel_volume(P_loit*loitering_time[i], thermal_efficiency=th_efficiency, energy_density=SAF_energy_density))
        E_loit += P_loit*loitering_time[i]
    
    #caculate total required fuel
    L_SAF_loit = float(convert_fuel_volume(E_loit, thermal_efficiency=th_efficiency, energy_density=SAF_energy_density))
    print("fuel for loitering:", L_SAF_loit )
    return L_SAF_loit


def descending_fuel(ROD, C_L_descent, C_d_0, AR, e, rho_array, SAF_energy_density, th_efficiency):
    """
    Args:
            ROD: float,                 planned rate of descent
            C_L_descent: float,         planned descent lift coefficient
            c_d_0: float,               zero lift drag
            AR: float,                  Aspect ratio of aircraft
            e: float between 0 and 1,   oswald efficiency factor
            rho_array: array of size n, density at different altitudes
    Outputs:
            l_SAF_descent: float,       kilograms of fuel required for descent
    """
    V = np.sqrt((2*W)/(rho_array*S*C_L_descent))

    C_L_descent = np.ones(len(V))*C_L_descent
    for i in range(len(V)):
        if V[i]>68:
            V[i]=68
            C_L_descent[i] = 2*W/(rho_array[i]*V[i]**2*S)


    C_D = C_d_0 + C_L_descent**2/(np.pi*AR*e)
    P_EH = ROD*W
    P = (0.5*V**3*rho_array*S*C_D-P_EH)/prop_efficiency
    for i in range(len(P)):
        if P[i]>2*P[0]:
            P[i] = 2*P[0]
        if P[i]>2*P[-1]:
            P[i] = 2*P[-1]
        if P[i]<0:
            P[i] = 0
    P = P/th_efficiency
    E = sum(P*1/(ROD))
    l_SAF_descent = convert_fuel_volume(E,th_efficiency, SAF_energy_density)
    print("fuel for descent:", l_SAF_descent)
    return l_SAF_descent


def convert_fuel_volume(required_energy, thermal_efficiency = 0.35, energy_density = 40*10**6):
    """
    args:
          required energy: float or array of size n,    the energy required for a certain operation 
          thermal efficiency: float between 0 and 1,    the efficiency of the thermal process, limited by the carnot efficiency
          energy density: float,                        energy per liter of the propellant
    returns:
          l_fuel: float or array of size n,             fuel volume required
    """
    l_fuel = required_energy/(thermal_efficiency*energy_density)
    return l_fuel





req_energy = []

#constant ROC strategy
ROC = 10    #[m/s]

#https://www.omicsonline.org/articles-images/2168-9792-5-161-g009.html
#2412 foil
cl_alpha = 0.11
cl_0 =0.25 
C_d_0 = 0.07
c_d_2 = 0.000675
W=250 #[N]
C_L_climb = 1.1
C_L_descent = 0.15
prop_efficiency = 0.7
electric_sys_eff = 0.7
AR= 8
e =0.9
S = None

#inputs
startaltitude = 0
top_altitude = 30_500
finish_altitude = 0
V_0 = 15
plots = False
SAF_energy_density = 40*10**6 #[j/kg]
loitering_time = [1500,750,750] #[s]
loitering_altitude = [30_500,10_000,5_000] #[m]
C_L_loiter = [1.2, 0.18, 0.12]
th_efficiency = 0.35
total_flight_time = 10_000 #[s]

#________________________________________________________________________________________________________________________________________________________________________________________
#actual model simulations
rho_array = ambiance.Atmosphere(np.arange(startaltitude,top_altitude)).density
weights = np.arange(250,450,500)

S_return, l_saf_return, P_max_return, min_E_tot_loc = climb_characteristics(weights, np.linspace(7.3,10.0,2), rho_array, ROC, V_0, C_d_0, AR, e, plots, SAF_energy_density, th_efficiency)

#taking the surface area of minimum energy, this is part of a design strategy, can be altered
S = S_return[min_E_tot_loc]
l_SAF = l_saf_return[min_E_tot_loc]
P_max = P_max_return[min_E_tot_loc]

#loitering
l_SAF +=loitering_fuel(loitering_time,loitering_altitude,C_L_loiter,S,AR,e,SAF_energy_density, th_efficiency)


#descending:
descent_time = total_flight_time-((top_altitude-startaltitude)/ROC+sum(loitering_time))
ROD = (top_altitude-finish_altitude)/descent_time
l_SAF+=descending_fuel(ROD, C_L_descent, C_d_0, AR, e, rho_array, SAF_energy_density, th_efficiency)
print("total kilograms of SAF required:",l_SAF)