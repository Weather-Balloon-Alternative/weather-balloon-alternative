from Balloon_with_motor import *

# Snowflake parafoil data
m_sfdry = 1.95 # [kg] - Snowflake dry mass
C_bat = 2100. # [mAh] - Battery capacity
P_full = 750. # [mA] - Power consumption (controls + downlink)
P_autopilot = 250. # [mA] - Power consumption (controls only)
V_y = 3.66 # m/s - descend rate
V_x = 7.2 # m/s - forward speed
glide_ratio = 2./1. # - Specific to Snowflake, probably closer to 3
R_min = 15.2 # m - Minimum turn radius

# Other parameters
m_pay = 2. # kg - Payload mass
m_parachute = 0.190 # kg - Parachute mass (high altitude sciences 1.5 m diameter)
alt_open = 17000 # m - Parafoil deployment altitude (crude estimation)
H2_GWP = 11 # - Equivalent global warming as compared to same amount of CO2
max_alt = 30500
H2_costperkg = 3. # $/kg - Cost of hydrogen per kg

# Balloon weight estimation
def weight_estimate(m_sfdry,m_pay,m_parachute):
    m_lifted = m_sfdry + m_pay + m_parachute
    m_balloon = balloon_mass_update(m_lifted)
    m_tot = m_lifted + m_balloon
    return m_tot,m_lifted,m_balloon

# Range estimation
def range_and_endurance_estimation(alt_open,glide_ratio, V_y, V_x):
    # Range based on glide ratio
    range1 = alt_open*glide_ratio

    # Range based on vertical and horizontal velocity
    flight_time = alt_open/V_y
    range2 = flight_time*V_x
    return range1, range2, flight_time


# Operational emissions
def operational_emisions(m_lifted, molar_mass):
    H2_emisions = calculate_required_size(m_lifted,0,molar_mass)
    CO2_equivalent = H2_GWP*H2_emisions
    return CO2_equivalent

# Consumable cost
def consumable_cost(m_lifted, molar_mass):
    H2_volume = calculate_required_size(m_lifted,0,molar_mass)
    rho = gas_density(molar_mass,0)
    H2_mass = rho*H2_volume
    H2_cost = H2_mass * H2_costperkg
    return H2_cost
# Operational cost?

# Transportational cost?


print('Mass = ',weight_estimate(m_sfdry, m_pay, m_parachute)[0])
print('Range = ',range_and_endurance_estimation(alt_open,glide_ratio,V_y,V_x)[0])
print('Endurance = ',range_and_endurance_estimation(alt_open,glide_ratio,V_y,V_x)[2])
print('Operational emissions = ',operational_emisions(weight_estimate(m_sfdry,m_pay,m_parachute)[1],M['H2']))
print('Consumable costs = ',consumable_cost(weight_estimate(m_sfdry,m_pay,m_parachute)[1],M['H2']))

print(weight_estimate(m_sfdry,m_pay,m_parachute)[2])