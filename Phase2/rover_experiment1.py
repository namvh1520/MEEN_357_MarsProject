from define_experiment import experiment1
from subfunctions import simulate_rover
import numpy as np
import matplotlib.pyplot as plt

## DICTIONARY
planet = {
        'g' : 3.72 #acceleration due of gravity in m/s^2
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

def createplots(rover):
    experiment, end_event = experiment1() #grabbing the experiment conditions from define_experiment
    rover_updated = simulate_rover(rover, planet, experiment, end_event) #grabbing updated rover dictionary from simulate_rover
    
    time = rover_updated['telemetry']['Time'] #grabbing lists out of the telemetry
    position = rover_updated['telemetry']['position']
    velocity = rover_updated['telemetry']['velocity']
    power = rover_updated['telemetry']['power']
    
    
    fig, axs = plt.subplots(3) #plotting in a subplot
    plt.subplots_adjust(hspace=1.2) #increasing distance between the three plots
    axs[0].plot(time, position)
    axs[0].set_title('Position vs Time')
    axs[0].set_ylabel('Position (m)')
    
    axs[1].plot(time, velocity)
    axs[1].set_title('Velocity vs Time')
    axs[1].set_ylabel('Velocity (m/s)')
    
    axs[2].plot(time, power)
    axs[2].set_title('Power vs Time')
    axs[2].set_ylabel('Power (W)')

    for ax in axs.flat: ax.set(xlabel='Time (s)') #giving all three subplots the same x label

createplots(rover)