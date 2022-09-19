import numpy as np
from scipy import special 

#ESTABLISHMENT OF DICTIONARIES


def create_dictionary():
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
        'power_subsys' : power_subsys
        }
    
    return planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover
    
#DEFINITION OF FUNCTIONS
def get_mass(rover):
    if (type(rover) is not dict):
        raise Exception("Argument passed must be a dictionary")
        
        
    mass = 0
    mass += rover['wheel_assembly']['wheel']['mass'] * 6 #There are 6 wheels
    mass += rover['wheel_assembly']['speed_reducer']['mass']
    mass += rover['wheel_assembly']['motor']['mass']
    mass += rover['chassis']['mass']
    mass += rover['power_subsys']['mass']
    mass += rover['science_payload']['mass']
    
    return mass #in kg

def get_gear_ratio(speed_reducer): #this function compute the gear ratio of the speed reducer
    
    #Type Checking
    if (type(speed_reducer) is not dict):
        raise Exception("Argument passed must be a dictionary")
        
    if (speed_reducer['type'].lower() != 'reverted'):
        raise Exception("Type of speed reducer should be reverted")
        
    #Calculating gear ratio
    Ng = (speed_reducer['diam_gear'] / speed_reducer['diam_pinion'])**2
    
    return Ng

def tau_dcmotor(omega, motor):
    if (type(motor) is not dict):
        raise Exception("Motor argument must be a dictionary")
    if (type(omega) is not np.ndarray):
        raise Exception("Omega argument must be an np.array")
        
        
    torque_stall = motor['torque_stall']
    torque_noload = motor['torque_noload']
    speed_noload = motor['speed_noload']
    
    
    tau = torque_stall - ((torque_stall - torque_noload)/speed_noload)*omega
    
    return tau

def F_drive(omega, rover):
    if (type(rover) is not dict):
        raise Exception("Motor argument must be a dictionary")
    if (type(omega) is not np.ndarray):
        raise Exception("Omega argument must be an np.array")
        
    radius = rover['wheel_assembly']['wheel']['radius']
    motor = rover['wheel_assembly']['motor']
    taudc = tau_dcmotor(omega, motor)
    tauwheel = taudc * get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    Fd = tauwheel/radius * 6
    
    return Fd

def F_gravity(terrain_angle, rover, planet):
    if (type(terrain_angle) is not np.ndarray):
        raise Exception("Terrain Angle must be an np.array")
    if (type(rover) is not dict):
        raise Exception("Rover must be a dictionary")
    if (type(planet) is not dict):
        raise Exception("Planet must be a dictionary")
        
    if np.any((terrain_angle <-75) | (terrain_angle > 75)):
        raise Exception("Terrain Angle must be betwen -75 and 75 degrees")
        
    Fgt = get_mass(rover) * planet['gravity'] * np.sin(terrain_angle * np.pi / 180)
    
    if terrain_angle < 0:
        return abs(Fgt)
    else:
        return abs(Fgt) * -1
    
    
def F_rolling(omega, terrain_angle, rover, planet, Crr):
    if (not isinstance(omega, np.ndarray)):
        raise Exception("Omega must be an np.array")
    if (not isinstance(terrain_angle, np.ndarray)):
        raise Exception("Terrain Angle must be an np.array")
    if (not isinstance(rover, dict)):
        raise Exception("Rover must be a dictionary")
    if (not isinstance(planet, dict)):
        raise Exception("Planet must be a dictionary")
    if (not isinstance(Crr, int) and not isinstance(Crr, float) or Crr < 0):
        raise Exception("Crr must be a scalar")
    
    if (omega.size != terrain_angle.size):
        raise Exception("Omega and terrain angle must be the same size")
    if np.any((terrain_angle < -75) | (terrain_angle > 75)):
        raise Exception("Terrain angle must be between -75 and 75")
        
        
    Fr1 = Crr * get_mass(rover) * planet['gravity'] * np.cos(terrain_angle * np.pi/180)
    
    Frr = special.erf(40 * rover['wheel_assembly']['wheel']['radius'] *  (omega/ get_gear_ratio(rover['wheel_assembly']['speed_reducer'])) ) * Fr1
    return abs(Frr) * -1 * 6
    
def F_net(omega, terrain_angle, rover, planet, Crr):
    if (not isinstance(omega, np.ndarray)):
        raise Exception("Omega must be an np.array")
    if (not isinstance(terrain_angle, np.ndarray)):
        raise Exception("Terrain Angle must be an np.array")
    if (not isinstance(rover, dict)):
        raise Exception("Rover must be a dictionary")
    if (not isinstance(planet, dict)):
        raise Exception("Planet must be a dictionary")
    if (not isinstance(Crr, int) and not isinstance(Crr, float) or Crr < 0):
        raise Exception("Crr must be a scalar")
    
    if (omega.size != terrain_angle.size):
        raise Exception("Omega and terrain angle must be the same size")
    if np.any((terrain_angle < -75) | (terrain_angle > 75)):
        raise Exception("Terrain angle must be between -75 and 75")
        
       
    #Frr = np.sqrt(F_drive(omega, rover)**2 + F_gravity(terrain_angle, rover, planet)**2 +   F_rolling(omega, terrain_angle, rover, planet, Crr)**2)
    Fr1 = F_drive(omega, rover) + F_gravity(terrain_angle, rover, planet) +   F_rolling(omega, terrain_angle, rover, planet, Crr)
    return Fr1
    