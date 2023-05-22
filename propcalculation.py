import ambiance
import numpy as np
import matplotlib.pyplot as plt

# Inputs
hmax = 33000 #m
climbrate = 2.5 #m/s
mass = 75 #kg
Vcraft = 55 #m/s
cd0 = 0.07
Aspectratio = 8
Wingarea = 15 #m2
oswaldfactor = 0.9 
propefficiency = 0.87
nengines = 6
Plotfunc = True

# Data from Torenbeek for 2 bladed propeller
J = 1.2
CP = 0.06
Mtipmax = 0.8
beta = 25 #deg


def Propeller_sizing(m, h, ROC, V, C_d_0, AR, S, e, prop_efficiency, cp, Mmax, n_engines,PLOT=False):
    '''
    Sizing for a two-blade propeller optimized for a given altitude and velocity. 
    Returns the required thrust, power, propeller diameter and the revolutions per second.
    args:
            m: float,                               system mass
            h: float,                               system altitude
            ROC: float,                             rate of climb that is selected 
            V: float,                               true airspeed
            C_d_0: float,                           zero lift drag
            AR: float,                              aspect ratio of the aircraft
            S: float,                               Wing surface area
            e: float between 0 and 1,               oswald efficiency factor
            prop_efficiency: float,                 propeller efficiency
            cp: float,                              propulsion coefficient
            Mmax: float,                            maximum propeller tip mach number
            n_engines: integer,                     number of engines on the aircraft
            PLOT: optional boolean                  set whether to plot, standardly set to False

    output:
            dictionary{
            P: float,                               required power
            T: float,                               required thrust
            d: float,                               propeller diameter
            n: float,                               rotational velocity
            }
    '''
    ###### Setup
    atmos = ambiance.Atmosphere(h)
    rho = atmos.density
    speedofsound = atmos.speed_of_sound
    theta = np.arcsin(ROC/V)

    # Solve FBD
    W = m*atmos.grav_accel
    L = W/np.cos(theta)
    c_L = L / (0.5*rho*V**2*S)
    c_d = C_d_0 + c_L/(np.pi*AR*e)
    D = c_d*0.5*rho*V**2*S

    Ttot = D + W*np.sin(theta)
    T = Ttot / n_engines

    # Calculate required power
    P = T*V/prop_efficiency

    # Rewritten equations from rotational speed - diameter relations (Torenbeek)
    def rewritten_powereq(d_r,P_r,rho_r,cp_r):
        # _r to distinguish between previously used variables
        return (P_r/(cp_r*rho_r*d_r**5))**(1/3)
    
    def rewritten_tipspeed(d_r,V_r,speedofsound_r,Mmax_r):
        # _r to distinguish between previously used variables
        return ((Mmax_r*speedofsound_r)**2 - V_r**2)**(1/2)/(np.pi*d_r)

    def solve(d_r):
        deltan = np.abs(rewritten_powereq(d_r,P,rho,cp) - rewritten_tipspeed(d_r,V,speedofsound,Mmax))
        return deltan

    # Find intersection (Ugly but it works)
    darray = np.linspace(0.1,20,100000)
    vals = solve(darray)
    id = np.where(vals == np.min(vals))
    d = darray[id]
    n = rewritten_powereq(d,P,rho,cp)

    # Assert similarity
    if np.abs(rewritten_powereq(d,P,rho,cp) - rewritten_tipspeed(d,V,speedofsound,Mmax)) > rewritten_powereq(d,P,rho,cp)*10**-4:
        print(f'WARNING: results not similar (power gives n={rewritten_powereq(d,P,rho,cp)} while tip speed gives n={rewritten_tipspeed(d,V,speedofsound,Mtipmax)})')

    # If selected, plot the resultant power and tip speed requirements along with their delta
    if PLOT:
        n1 = rewritten_powereq(darray,P,rho,cp)
        n2 = rewritten_tipspeed(darray,V,speedofsound,Mmax)
        plt.plot(darray,n1,label='power')
        plt.plot(darray,n2,label='tipspeed')
        plt.plot(darray,vals,label='delta')
        plt.legend()
        plt.tight_layout
        plt.show()


    return {
        'T': T,
        'P': P,
        'n': n,
        'd': d
    }

print(Propeller_sizing(mass,hmax,climbrate,Vcraft,cd0,Aspectratio,Wingarea,oswaldfactor,propefficiency,CP,Mtipmax,nengines,PLOT=Plotfunc))