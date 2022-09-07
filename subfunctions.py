import numpy as np

#ESTABLISHMENT OF DICTIONARIES



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
    'torque_noload' : 0, #toque in Nm
    'speed_noload' : 3.8, #speed in rad/s
    'mass' : 5 #mass in kg
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
    'pwer_subsys' : power_subsys
    }

#DEFINITION OF FUNCTIONS
def get_mass(rover):
    if (type(rover) is not dict):
        raise Exception("Argument passed must be a dictionary")
        
        
    mass = 0
    mass += rover['wheel_assembly']['wheel']['mass'] * 6 #There are 6 wheels
    mass += rover['wheel_assembly']['speed_reducer']['mass']
    mass += rover['wheel_assembly']['motor']['mass']
    mass += rover['chassis']['mass']
    mass += rover['pwer_subsys']['mass']
    mass += rover['science_payload']['mass']
    
    return mass

def get_gear_ratio(speed_reducer):
    if (type(speed_reducer) is not dict):
        raise Exception("Argument passed must be a dictionary")
        
    if (speed_reducer['type'].lower() != 'reverted'):
        raise Exception("Type of speed reducer should be reverted")
        
    Ng = speed_reducer['diam_pinion'] / speed_reducer['diam_gear']
    
    return Ng

def tau_dcmotor(omega, motor):
    if (type(motor) is not dict):
        raise Exception("Motor argument must be a dictionary")
    if (type(omega) is not np.ndarray and type(omega) is not int and type(omega) is not float):
        raise Exception("Omega argument must be an np.array or a scalar")


def F_drive():
    print('F_drive')

def F_gravity():
    print('F_gravity')
    
def F_rolling():
    print('F_rolling')
    
def F_net():
    print('F_net')
    
