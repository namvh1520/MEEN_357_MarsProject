from define_experiment import experiment1
from subfunctions import simulate_rover
import numpy as np
import matplotlib.pyplot as plt

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
    experiment, end_event = experiment1()
    rover_updated = simulate_rover(rover, planet, experiment, end_event)
    
    time = rover_updated['telemetry']['Time']
    position = rover_updated['telemetry']['position']
    velocity = rover_updated['telemetry']['velocity']
    power = rover_updated['telemetry']['power']
    
    
    plt.plot(time, position)
    plt.ylabel('Positon (m)')
    plt.xlabel('Time (s)')
    plt.title('Position vs Time')
    plt.show()
    
    plt.plot(time, velocity)
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Velocity vs Time')
    plt.show()
    
    plt.plot(time, power)
    plt.xlabel('Time (s)')
    plt.ylabel('Power (m/s)')
    plt.title('Power vs Time')
    plt.show()

createplots(rover)