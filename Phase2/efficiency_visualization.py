import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

planet = {
        'gravity' : 3.72 #acceleration due of gravity in m/s^2
        }

power_subsys = {
        'mass' : 90 #mass in kg
        }

science_payload = {
        'mass' : 75 #mass in kg
        }

chassis = {
        'mass' : 659 #mass in kg
        }

motor = {
        'torque_stall' : 170, #torque in Nm
        'torque_noload' : 0, #torque in Nm
        'speed_noload' : 3.8, #speed in rad/s
        'mass' : 5, #mass in kg
        'effcy_tau' : np.array([0, 10, 20, 40, 75, 165]),
        'effcy' : np.array([0, .6, .75, .73, .55, .05])
        }

speed_reducer = {
        'type' : 'reverted', #only valid type in phase 1
        'diam_pinion' : .04, #diamenter in m
        'diam_gear' : .07, #diamenter in m
        'mass' : 1.5 #mass in kg
        }

wheel = {
        'radius' : .3, #radius in m
        'mass' : 1 #mass in kg
        }

wheel_assembly = {
        'wheel' : wheel,
        'speed_reducer' : speed_reducer,
        'motor' : motor  
        }
    
rover = {
        'wheel_assembly' : wheel_assembly,
        'chassis' : chassis,
        'science_payload' : science_payload,
        'power_subsys' : power_subsys
        }


def plot_effcy(rover):
    effcy = rover['wheel_assembly']['motor']['effcy'] #grabbing efficiency list
    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau'] #grabbing efficiency tau list
    
    alpha_fun = interp1d(effcy_tau, effcy, kind = 'cubic')
    
    x = np.linspace(effcy_tau[0], effcy_tau[-1], 101) #setting up x values for 100 points
    
    plt.plot(effcy_tau, effcy, 'r*', x, alpha_fun(x)) #plotting motor torque
    plt.xlabel('Motor Torque')
    plt.ylabel('Motor Efficiency')
    plt.title('Motor Torque vs Motor Efficiency')
    
plot_effcy(rover)