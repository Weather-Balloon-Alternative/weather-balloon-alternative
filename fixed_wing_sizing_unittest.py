import pytest
import fixed_wing_sizing
import numpy as np

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
    assert fixed_wing_sizing.descending_fuel(ROD, C_L_descent, C_D_0, AR, e, rho_array, energy_density, th_efficiency, W, S, prop_efficiency)-1805.48/energy_density<0.1