from ambiance import Atmosphere
import numpy as np
from matplotlib import pyplot as plt

#parameters
M = {'H2':2.015894, 'He':4.0026022}
CD = 0.7 #sphere is 0.47
prop_efficiency = 0.7
motor_efficiency = 0.9

#mission
max_alt = 33000 #target altitude
m_pay = 2 #kg

#input
ascent_time = 1.6 #hours
descent_time = 6 #hours
Prop_power = 250 #W
m_prop = 0.5
m_balloon_assumed = 5 #varies with weight, will go out of hand
m_struct = 1


def gas_density(molar_mass, altitude):
    R = 8.3144598
    atmos = Atmosphere(altitude)
    rho = (atmos.pressure * molar_mass * 0.001) / (R * atmos.temperature)
    return rho

def calculate_required_size(m_tot, altitude, molar_mass):
    atmos = Atmosphere(altitude)
    delta_rho = atmos.density - gas_density(molar_mass, altitude)
    volume = m_tot / delta_rho
    return volume

"""
def calc_drag_force(m_tot, altitude, velocity):
    atmos = atmosphere(altitude)
    vol = calculate_required_size(m_tot, altitude, M['H2'])
    r = ((3 * vol) / (4 * np.pi)) ** (1 / 3)
    A_front = np.pi * (r ** 2)
    rho = atmos.density
    D = 0.5 * CD * rho * (velocity ** 2) * A_front
    return D
"""

def calc_range(descent_rate, m):
    alt_step = 100 #m
    alt = np.arange(0, max_alt, alt_step)
    timestep = alt_step/descent_rate
    range = float(0)
    V_log = []
    for i in alt:
        atmos = Atmosphere(i)
        rho = atmos.density
        volume = calculate_required_size(m, i, M['H2'])
        r = ((3 * volume) / (4 * np.pi)) ** (1 / 3)
        A_front = np.pi * (r ** 2)
        Power_eff = Prop_power*prop_efficiency*motor_efficiency
        V = (Power_eff/(0.5*CD*rho*A_front))**(1/3)
        range = range + V*timestep
        V_log.append(V)
    return range, V_log, alt
def balloon_mass_update(m_tot):
    volume = calculate_required_size(m_tot, max_alt, M['H2'])
    mass = (4.25*volume + 340) * 0.001
    return mass



def plottingpower():
    power = np.arange(100, 1100, 100)
    range_data = np.zeros(np.shape(power))
    V_avg_data = np.zeros(np.shape(power))
    m_tot_data = np.zeros(np.shape(power))
    for i in range(np.size(power)):

        Prop_power = power[i]
        E_req = Prop_power * descent_time
        wh_kg = 210

        m_bat = E_req / wh_kg

        # balloon mass iteration
        m_balloon_log = []
        m_balloon_iteration = m_balloon_assumed
        for j in range(10):
            m_tot_initial = m_pay + m_struct + m_prop + m_bat + m_balloon_iteration
            m_balloon_iteration = balloon_mass_update(m_tot_initial)
            m_balloon_log.append(m_balloon_iteration)

        m_b_final = m_balloon_iteration
        m_tot = m_pay + m_prop + m_struct + m_bat + m_b_final

        descent_rate = max_alt / (descent_time * 60 * 60)
        R, V_log = calc_range(descent_rate, m_tot)
        V_avg = sum(V_log)/len(V_log)
        range_data[i] = R
        V_avg_data[i] = V_avg
        m_tot_data[i] = m_tot


    #plot 2 in one
    fig, ax1 = plt.subplots()

    color1 = 'tab:red'
    ax1.set_xlabel('Propulsion power used [W]')
    ax1.set_ylabel('Range [m]', color=color1)
    ax1.plot(power, range_data, color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color2 = 'tab:blue'
    ax2.set_ylabel('Total mass [kg]', color=color2)  # we already handled the x-label with ax1
    ax2.plot(power, m_tot_data, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

def plottingdescenttime():
    descent_time_data = np.arange(0.5, 12, 0.5)
    range_data = np.zeros(np.shape(descent_time_data))
    V_avg_data = np.zeros(np.shape(descent_time_data))
    m_tot_data = np.zeros(np.shape(descent_time_data))
    for i in range(np.size(descent_time_data)):
        descent_time = descent_time_data[i]
        E_req = Prop_power * descent_time
        wh_kg = 210

        m_bat = E_req / wh_kg

        # balloon mass iteration
        m_balloon_log = []
        m_balloon_iteration = m_balloon_assumed
        for j in range(10):
            m_tot_initial = m_pay + m_struct + m_prop + m_bat + m_balloon_iteration
            m_balloon_iteration = balloon_mass_update(m_tot_initial)
            m_balloon_log.append(m_balloon_iteration)

        m_b_final = m_balloon_iteration
        m_tot = m_pay + m_prop + m_struct + m_bat + m_b_final

        descent_rate = max_alt / (descent_time * 60 * 60)
        R, V_log, alt = calc_range(descent_rate, m_tot)

        V_avg = sum(V_log)/len(V_log)
        range_data[i] = R
        V_avg_data[i] = V_avg
        m_tot_data[i] = m_tot
    #plot 2 in one
    fig, ax1 = plt.subplots()

    color1 = 'tab:red'
    ax1.set_xlabel('Descent time [hrs]')
    ax1.set_ylabel('Range [m]', color=color1)
    ax1.plot(descent_time_data, range_data, color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color2 = 'tab:blue'
    ax2.set_ylabel('Average velocity [m/s]', color=color2)  # we already handled the x-label with ax1
    ax2.plot(descent_time_data, V_avg_data, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

def calc_parameters():
    E_req = Prop_power*descent_time
    wh_kg = 210

    m_bat = E_req/wh_kg

    #balloon mass iteration
    m_balloon_log = []
    m_balloon_iteration = m_balloon_assumed
    for i in range(10):
        m_tot_initial = m_pay + m_struct + m_prop + m_bat + m_balloon_iteration
        m_balloon_iteration = balloon_mass_update(m_tot_initial)
        m_balloon_log.append(m_balloon_iteration)

    m_b_final = m_balloon_iteration
    m_tot = m_pay + m_prop + m_struct + m_bat + m_b_final
    m_module = m_pay + m_prop + m_bat + m_struct

    descent_rate = max_alt/(descent_time*60*60)
    R, V_log, alt = calc_range(descent_rate, m_tot)

    alt_list = alt.tolist()

    #emmisions:
    max_volume = calculate_required_size(m_tot, max_alt, M['H2'])
    gas_density_at_altitude = gas_density(M['H2'], max_alt)
    gas_mass_max = max_volume*gas_density_at_altitude
    min_volume = calculate_required_size(m_tot, 0, M['H2'])
    gas_density_at_sea = gas_density(M['H2'], 0)
    gas_mass_min = min_volume * gas_density_at_sea
    gas_emitted = gas_mass_max - gas_mass_min


    plt.plot(alt_list, V_log, label='Velocity')
    plt.ylabel('Velocity [m/s]')
    plt.xlabel('Altitude [m]')
    plt.gca().invert_xaxis()
    #plt.legend()
    #plt.show()

    V_avg = sum(V_log)/len(V_log)
    V_max = max(V_log)
    V_min = min(V_log)
    print("Range: ", R)
    print("mtot: ", m_tot)
    print("mbat: ", m_bat)
    print("M_balloon initial: ", m_balloon_assumed)
    print("final balloon mass: ",m_b_final)
    print("m_module: ", m_module)
    print("V_avg: ", V_avg)
    print("V_max: ", V_max)
    print("V_min: ", V_min)
    print("gas emitted: ", gas_mass_min)

plottingdescenttime()