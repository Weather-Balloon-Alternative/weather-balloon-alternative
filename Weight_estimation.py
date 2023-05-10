from Balloon_with_motor import *

# Snowflake parafoil data
M_sfdry = 1.95 # [kg] - Snowflake dry mass
C_bat = 2100 # [mAh] - Battery capacity
P_full = 750 # [mA] - Power consumption (controls + downlink)
P_autopilot = 250 # [mA] - Power consumption (controls only)
V_y = 3.66 # m/s - descend rate
V_x = 7.2 # m/s - forward speed
glide_ratio = 2/1 # -
R_min = 15.2 # m - Minimum turn radius

# Other parameters
M_pay = 2. # kg - Payload mass
M_parachute = 0.190 # kg - Parachute mass

# Balloon weight estimation
m_tot = M_sfdry + M_pay + M_parachute
M_balloon = balloon_mass_update(m_tot)
print(M_balloon)
