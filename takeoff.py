import ambiance
import numpy as np
g = 9.81
S = 2.0
Mach = 0.8
CD_0 = 0.07
W = 1000

def takeoff(max_z_alt, dz, S, Mach, CD_0, W):
    # altitude array
    arr = np.arange(0, max_z_alt, dz)
    # speed of sound at every altitude
    sos = ambiance.Atmosphere(arr).speed_of_sound
    V = sos * Mach
    # climb speed
    rho = ambiance.Atmosphere(arr).density
    # time at altitudes
    t = 1 / V
    # drag at altitude
    D = 0.5 * rho * V ** 2 * S * CD_0

    I_SP = 330
    T_array = []
    for i in range(len(arr)):
        T = D[i] * W * t[i]
        m_dot = T / (I_SP * 9.81)
        W -= m_dot * t[i] * 9.81
        T_array.append(T)
    #      time at altitude, final weight, max avg thrust
    return t, W, max(T_array)







'''import ambiance
import numpy as np
arr = np.arange(0,33000,1)
S = 2.0
sos= ambiance.Atmosphere(arr).speed_of_sound
V = sos*0.8
rho = ambiance.Atmosphere(arr).density
C_D_0 = 0.07
g = 9.81
t = 1/V
W = 1000
print("total time to climb: ", t)
D = 0.5*rho*V**2*S*C_D_0

I_SP = 330
T_array=[]
for i in range(len(arr)):
    T = D[i]*W*t[i]
    m_dot = T/(I_SP*9.81)
    W -= m_dot*t[i]*9.81
    T_array.append(T)
print("final structure+payload weight: ", W,"[N]")
print("maximum thrust required: ",max(T_array),"[N]")



# T = D+ W*t
# print("Thrust required: ",T)
# I_SP = 330
# m_dot = T/(I_SP*9.81)
# print("total mass for climb: ",sum(m_dot*t))
# print("average thrust required: ", np.average(T))'''