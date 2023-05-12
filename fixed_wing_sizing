import numpy as np 
import math
import pandas
import ambiance
import matplotlib.pyplot as plt

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
approach = "const_ROC" #const_ROC, const_P
plots = False
SAF_energy_density = 40*10**6
loitering_time = [1500,750,750] #[s]
loitering_altitude = [30_500,10_000,5_000] #[m]
C_L_loiter = [1.2, 0.18, 0.12]
th_efficiency = 0.35
total_flight_time = 10_000 #[s]

rho_array = ambiance.Atmosphere(np.arange(startaltitude,top_altitude)).density
E_tot_array = []
S_array = []
l_SAF = 0
slow_start_height = 1000
weights = np.arange(250,450,500)
dist = 0 
#ascent
if approach == "const_ROC":
    for W in weights:
        for i in np.linspace(73,74,2):
            S = i/10
            C_L = 2*W/(rho_array*S*V_0**2)
            for j in range(len(C_L)):
                if C_L[j]>C_L_climb:
                    C_L[j]=C_L_climb
            V = np.sqrt(2*W/(rho_array*S*C_L))
            C_D = C_d_0 + C_L**2/(np.pi*AR*e)
            #Power consumption calculations
            P_equilibrium=0.5*rho_array*S*V**3*C_D/prop_efficiency
            P_climb = np.ones(len(C_L))*ROC*W/prop_efficiency
            P_climb[0:slow_start_height] = (ROC/3)*W/prop_efficiency
            P_tot = P_equilibrium+P_climb
            # print(P_tot)
            T = np.ones(len(C_L))*1/ROC
            T[0:slow_start_height]=1/(ROC/3)
            E_tot = P_tot * T
            # print(np.max(P_tot))
            # print(np.sum(E_tot))
            # plt.plot(P_tot)
            E_tot_array.append(sum(E_tot))
            S_array.append(S)
        if plots:
            plt.plot(S_array,E_tot_array, label=str("Weight="+str(W)+"[N]"))
        min_E_tot_loc = int(np.argwhere(E_tot_array==min(E_tot_array)))
        print("Weight: ", W)
        print("optimum surface area: ",S_array[min_E_tot_loc])
        S = S_array[min_E_tot_loc]
        print("minimum energy required : ",E_tot_array[min_E_tot_loc])
        l_SAF +=E_tot_array[min_E_tot_loc]/(SAF_energy_density*th_efficiency)
        print("liters SAF required for ascent: ", l_SAF)
        print("maximum power requried: ", max(P_tot))
        print("")
        dist +=sum(V*1/ROC)
        S_array = []
        E_tot_array =[]
if plots:
    plt.legend()
    plt.xlabel("S[m^2]")
    plt.ylabel("Energy[J]")
    plt.show()
print("distance: ", dist)

#loitering:
E_loit = 0
for i in range(len(loitering_time)):
    alt = loitering_altitude[i]
    rho = ambiance.Atmosphere(alt).density
    V = np.sqrt((2*W)/(rho*C_L_loiter[i]*S))
    C_D = C_d_0 + C_L_loiter[i]**2/(np.pi*AR*e)
    P_loit = 1/2*V**3*rho*S*C_D/prop_efficiency
    print("for phase", i, " fuel required: ", P_loit*loitering_time[i]/(th_efficiency*SAF_energy_density))
    E_loit += P_loit*loitering_time[i]
    dist +=int(V*loitering_time[i])
print("distance: ", dist)
L_SAF_loit = float(E_loit/(th_efficiency*SAF_energy_density))
print("fuel for loitering:", L_SAF_loit )
print(l_SAF)
l_SAF +=L_SAF_loit
print(l_SAF)
#descending:
descent_time = total_flight_time-((top_altitude-startaltitude)/ROC+sum(loitering_time))
ROD = (top_altitude-finish_altitude)/descent_time
V = np.sqrt((2*W)/(rho_array*S*C_L_descent))

C_L_descent = np.ones(len(V))*C_L_descent
for i in range(len(V)):
    if V[i]>68:
        V[i]=68
        C_L_descent[i] = 2*W/(rho_array[i]*V[i]**2*S)

dist += sum(V*1/ROD)
print("distance total: ", dist)

C_D = C_d_0 + C_L_descent**2/(np.pi*AR*e)
P_EH = ROD*W
P = (0.5*V**3*rho_array*S*C_D-P_EH)/prop_efficiency
for i in range(len(P)):
    if P[i]>2*P[0]:
        P[i] = 2*P[0]
    if P[i]>2*P[-1]:
        P[i] = 2*P[-1]
P = P/th_efficiency
E = sum(P*1/(ROD))
print("fuel for descent:", E/SAF_energy_density)
l_SAF +=E/SAF_energy_density
print("total SAF required: ", l_SAF)

