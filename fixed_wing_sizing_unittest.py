import pytest
import fixed_wing_sizing
import numpy as np
import ambiance

def test_fuel_volume():
    required_energy =10_000
    thermal_efficiency = 0.50
    energy_density = 20_000
    assert fixed_wing_sizing.convert_fuel_volume(required_energy, thermal_efficiency, energy_density) == 1

def test_descent_fuel():
    ROD =0.1
    C_L_descent = 1.0
    C_D_0 = 0.10
    AR  = 10
    e = 0.5
    rho_array = np.ones(10)
    energy_density = 10_000
    th_efficiency = 0.5
    W = 100
    S = 10
    prop_efficiency = 0.7
    #power from energyheight: 10 W
    #speed: (2*W/(rho*CL*S))^0.5 = 4.47
    #Power: 0.5*rho*V**3*S*(C_D_0+CL^2/(pi*AR*e)) = 73.1919
    #energy: (power*t- ROD*W*t)/(th_efficiency*prop_efficiency) = 1805.48
    # print(fixed_wing_sizing.descending_fuel(ROD, C_L_descent, C_D_0, AR, e, rho_array, energy_density, th_efficiency, W, S, prop_efficiency))
    assert abs(fixed_wing_sizing.descending_fuel(ROD, C_L_descent, C_D_0, AR, e, rho_array, energy_density, th_efficiency, W, S, prop_efficiency=prop_efficiency)-18054.8/ (energy_density))<0.01

def test_loiter_fuel():
    loitering_time = [100,100]
    loitering_altitude = [0, 1000]
    C_L_loiter = [1,1]
    S = 10
    W = 100
    AR = 10
    e = 0.5
    energy_density = 10_000
    th_efficiency = 0.5
    prop_efficiency = 0.7
    C_d_0 =0.10
    #C_D = 0.10+ 1^2/(3.14*10*0.5) = 0.16366
    #V_1 = 4.0406
    #P_1 = 94.469
    #E_1 = 9446.9
    
    #V_2 = 4.2416
    #P_2 = 99.168
    #E_2 = 9916.8

    #E_tot = 19363.7
    #l_fuel = 19363.7/(th_efficiency*energy_density) = 3.87274 

    assert abs(fixed_wing_sizing.loitering_fuel(loitering_time,loitering_altitude, C_L_loiter, S, AR, e, energy_density, th_efficiency, W, prop_efficiency, C_d_0) - 3.87274)<0.01

def test_climb_characteristics():
    W = [100,1000]
    S = [1,10]
    rho_array = np.ones(100)
    ROC = 10
    V_0 = 1
    C_d_0 = 0.10
    AR = 10
    e= 0.5
    energy_density = 10_000
    th_efficiency = 0.5
    C_L_climb = 1
    prop_efficiency = 0.7
    # V_1 = 4.4721
    # V_2 = 14.142
    # C_D = 0.16366
    # P_EH_1 = 1_000/0.7 = 1428.57
    # P_EH_2 = 10_000/0.7 = 14285.7
    # P_eq_1 = 104.5585
    # P_eq_2 = 3306.336
    # P_tot_1 = 1533.12
    # P_tot_2 = 17592.19
    # E_min_1 = 15331.2
    # E_min_2 = 175921.9
    # l_saf_1 = 3.0662
    # l_saf_2 = 35.18438
    print("returndata")
    print(fixed_wing_sizing.climb_characteristics(W, S, rho_array, ROC, V_0, C_d_0, AR, e, False, energy_density, th_efficiency, C_L_climb, prop_efficiency))
    actual_data = fixed_wing_sizing.climb_characteristics(W, S, rho_array, ROC, V_0, C_d_0, AR, e, False, energy_density, th_efficiency, C_L_climb, prop_efficiency) 
    expected_data= [[10,10],[3.0662,35.18438],[1533.13,17592.19],[1]]
    print("actual data", actual_data)
    print("expected data", expected_data)
    array = []
    for i in range(len(expected_data)):
        for j in range(len(expected_data[i])):
            if abs(expected_data[i][j]-actual_data[i][j])<0.01:
                array.append(True)
            else:
                array.append(False)
    assert sum(array)==len(array)


if __name__ =="__main__":
    print(ambiance.Atmosphere(0).density)
    print(ambiance.Atmosphere(1000).density)
    print(np.ones(10))
    print(np.array([[0,0],[1,1]])-np.array([[1,1],[2,2]]))