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
def get_mass(r):
    print('get_mass')

def get_gear_ratio():
    print('get_gear_ratio')

def tau_dcmotor():
    print('tau_dcmotor')

def F_drive():
    print('F_drive')

def F_gravity():
    print('F_gravity')
    
def F_rolling():
    print('F_rolling')
    
def F_net():
    print('F_net')