# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 17:57:10 2022

@author: namqv_7yst908
"""

import numpy as np

def create_dictionary():
    planet = {
        'g' : 3.72 #acceleration due of gravity in m/s^2
        }

    power_subsys = {
        'mass' : 90.0 #mass in kg
        }

    science_payload = {
        'mass' : 75.0 #mass in kg
        }

    chassis = {
        'mass' : 659.0 #mass in kg
        }

    motor = {
        'torque_stall' : 170, #torque in Nm
        'torque_noload' : 0.0, #torque in Nm
        'speed_noload' : 3.8, #speed in rad/s
        'mass' : 5.0, #mass in kg
        'effcy_tau' : np.array([0, 10, 20, 40, 75, 165]), #N- element array of torque datapoints at which efficiency data is gathered [N-m]
        'effcy' : np.array([0, 0.60, 0.75, 0.73, 0.55, 0.05])
        }
    
    telemetry = {
        'Time' : np.array([]),
        'completion_time' : None,
        'velocity' : np.array([]),
        'position' : np.array([]),
        'distance_traveled' : None,
        'max_velocity' : None, 
        'average_velocity' : None,
        'power' : np.array([]),
        'battery_energy' : None, 
        'energy_per_distance' : None
        
        }

    speed_reducer = {
        'type' : 'reverted', #only valid type in phase 1
        'diam_pinion' : .04, #diamenter in m
        'diam_gear' : .07, #diamenter in m
        'mass' : 1.5 #mass in kg
        }

    wheel = {
        'radius' : .3, #radius in m
        'mass' : 1.0 #mass in kg
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
        'power_subsys' : power_subsys,
        'telemetry' : telemetry
        }
    
    return planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover
