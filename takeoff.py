import ambiance
import numpy as np
import matplotlib.pyplot as plt

def integrate_requied_power(D, W, I_SP, TtoW_engine):
    """
    Args:
            D: array of size n,         drag that the vehicle experiences at specific flight level
            W: float,                   initial weight
            I_SP: float,                specific impulse
            TtoW_engine: float,         thrust to weight of the engine
    outputs:
            W: float,                   final weight
            T_max: float,               maximum thrust required
            W_engine: float,            engine weight
    """
    T_array=[]
    #calculate the thrust required for every point in time
    for i in range(len(arr)):
        T = D[i]+W
        m_dot = T/(I_SP*9.81)
        W -= m_dot*t[i]*9.81
        T_array.append(T)
    T_max = max(T_array)
    W_engine = T_max/TtoW_engine
    print("final structure+payload weight: ", W,"[N]")
    print("maximum thrust required: ",T_max,"[N]")
    print("rocket size required: ",W_engine,"[N]")
    return W, T_max, W_engine

arr = np.arange(0,30500,1)
S = 2.0
sos= ambiance.Atmosphere(arr).speed_of_sound
V = sos*0.8
rho = ambiance.Atmosphere(arr).density
C_D_0 = 0.07
g = 9.81
t = 1/V
W_init = 1000
W = W_init
TtoW_engine = 180
I_SP = 330
print("total time to climb: ", sum(t))
D = 0.5*rho*V**2*S*C_D_0
integrate_requied_power(D, W_init, I_SP, TtoW_engine)

# T = D+ W*t
# print("Thrust required: ",T)
# I_SP = 330
# m_dot = T/(I_SP*9.81)
# print("total mass for climb: ",sum(m_dot*t))
# print("average thrust required: ", np.average(T))'''