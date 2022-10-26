"""###########################################################################
#   This file contains subfunctions for Phase 2 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Engineering Analysis Team
#   Last Modified: 6 October 2022
###########################################################################"""

import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
#from scipy.integrate import simpson
from statistics import mean

chassis = {'mass': 659}
science_payload = {'mass': 75}
power_subsys = {'mass': 90}
planet = {'g': 3.72}
wheel = {'radius': 0.30, 'mass': 1.0}
telemetry = {'Time': np.array([]), 'completion_time': 4 ,'velocity': np.array([]),'position':np.array([]),'distance_traveled': 4,'max_velocity': 4,'average_velocity': 4,'power': np.array([]),'battery_energy': 4,'energy_per_distance': 4} #NOT RIGHT
speed_reducer = {'type': 'reverted', 'diam_pinion': 0.04, 'diam_gear': 0.07, 'mass': 1.5}
motor = {'torque_stall': 170, 'torque_noload': 0, 'speed_noload': 3.80, 'mass': 5.0, 'effcy_tau': np.array([0,10,20,40,75,165]),'effcy': np.array([0,0.60,0.75,0.73,0.55,0.05])}
wheel_assembly = {'wheel': wheel, 'speed_reducer': speed_reducer, 'motor': motor}
rover = {'wheel_assembly': wheel_assembly, 'chassis': chassis, 'science_payload': science_payload, 'power_subsys': power_subsys,'telemetry':telemetry}
experiment = {'time_range': np.array([]),'initial_conditions': np.array([]),'alpha_dist': np.array([]),'alpha_deg': np.array([]),'Crr':0.1}
end_event = {'max_distance':50,'max_time':5000,'min_velocity':0.01}

omega = np.linspace(0,motor['speed_noload'],25)
terrain_angle = np.linspace(-75,75,25)
Crr = 0.2

v_scalar = 5
v_vector = np.linspace(0,50,25)
v_string = 'HI'


def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    # Check that 1 input has been given.
    #   IS THIS NECESSARY ANYMORE????
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """
    
    # Check that 2 inputs have been given
    #   IS THIS NECESSARY ANYMORE????
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:   dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    # Check that 2 inputs have been given.
    #   IS THIS NECESSARY ANYMORE????
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                     planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance 
                                              coefficient [-]
    
    Outputs:           Fnet:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet


def motorW(v, rover):
    """
    Inputs:               v:  numpy array     Array of velocities [m/s]
                      rover:  dict            Data structure specifying rover 
                                              parameters
    
    Outputs:              w:  numpy array     Array of motor speeds [rad/s]
    """
    
    # Check that the v input is a scalar or a vector
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise Exception('v input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float) # make the scalar a numpy array
    elif len(np.shape(v)) != 1:
        raise Exception('v input must be a scalar or a vector. Matrices are not allowed.')
    
    # Check that the rover input is a dict
    if type(rover) != dict:
        raise Exception('rover input must be a dict')
        
        
    # Main Code
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    w = v*Ng/r
    
    return w


def rover_dynamics(t, y, rover, planet, experiment):
    """
    Inputs:         t:  scalar            Time sample [s]
                    y:  numpy array       Two element array of dependent variables 
                                          (i.e., state vector). First element is 
                                          rover velocity [m/s] and second 
                                          element is rover position [m]
                rover:  dict              Data structure specifying rover 
                                          parameters
               planet:  dict              Data dictionary specifying planetary 
                                          parameters
           experiment:  dict              Data dictionary specifying experiment 
                                          definition
    
    Outputs:     dydt:  numpy array       First derivatives of state vector. 
                                          First element is rover acceleration 
                                          [m/s^2] and second element is rover 
                                          velocity [m/s]
    """
    
    # Check that the t input is a scalar
    if (type(t) != int) and (type(t) != float) and (not isinstance(t, np.ndarray)) and (not isinstance(t,np.float64)):
        raise Exception('t input must be a scalar.')
    elif isinstance(t, np.ndarray):
        if len(t) == 1:
            t = float(t) # make a scalar
        else:
            raise Exception('t input must be a scalar.')
    
    # Check that y input is a 2-entry numpy array
    if (not isinstance(y, np.ndarray)) or (len(y) != 2):
        raise Exception('y must be a 2x1 numpy array.')
    elif isinstance(y[0], np.ndarray):
        y = np.array([float(y[0]),float(y[1])]) # this will turn the column vector into a row
    
    # Check that the rover input is a dict
    if type(rover) != dict:
        raise Exception('rover input must be a dict')
    
    # Check that the planet input is a dict
    if type(planet) != dict:
        raise Exception('planet input must be a dict')
    
    # Check that the experiment input is a dict
    if type(experiment) != dict:
        raise Exception('experiment input must be a dict')
    
    # Main code
    v = float(y[0]) # velocity
    x = float(y[1]) # position
    
    omega = motorW(v, rover)   
    alpha_fun = interp1d(experiment['alpha_dist'].ravel(), experiment['alpha_deg'].ravel(), kind = 'cubic', fill_value="extrapolate")
    terrain_angle = float(alpha_fun(x))
    F = F_net(omega, terrain_angle, rover, planet, experiment['Crr'])
    
    m = get_mass(rover)
    accel = float(F/m)
    dydt = np.array([accel, v], dtype = float)
    
    return dydt


def mechpower(v, rover):
    """
    Inputs:               v:  numpy array     Array of velocities [m/s]
                      rover:  dict            Data structure specifying rover 
                                              parameters
    
    Outputs:              P:  numpy array     Array of instantaneous power 
                                              output of a single motor 
                                              corresponding to each element in 
                                              array v [W]
    """
    # Check that the v input is a scalar or a vector
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise Exception('v input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float) # make the scalar a numpy array
    elif len(np.shape(v)) != 1:
        raise Exception('v input must be a scalar or a vector. Matrices are not allowed.')
    
    # Check that the rover input is a dict
    if type(rover) != dict:
        raise Exception('rover input must be a dict')
    
    omega = motorW(v, rover)  
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor']) 
    
    P = tau*omega
    
    return P


def battenergy(t,v,rover):
    """
    Inputs:               t:  numpy array     Array of time samples from a 
                                              rover simulation [s]
                          v:  numpy array     Array of velocities from a rover 
                                              simulation [m/s]
                      rover:  dict            Data structure specifying rover 
                                              parameters
    
    Outputs:              E:  scalar          Total electrical energy consumed 
                                              from the rover battery pack over
                                              the input simulation profile [J]
    """
    # Check that the t input is a scalar or a vector
    if (not isinstance(t, np.ndarray)):
        raise Exception('t input must be a scalar or a vector. If t is a vector, it should be defined as a numpy array.')
    elif len(np.shape(t)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the v input is a scalar or a vector
    if (not isinstance(v, np.ndarray)):
        raise Exception('v input must be a scalar or a vector. If v is a vector, it should be defined as a numpy array.')
    elif len(np.shape(v)) != 1:
        raise Exception('v input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(t) != len(v):
        raise Exception('First two inputs must be the same size')
    
    # Main code
    P = mechpower(v, rover) # calculate power at each time/velocity
    omega = motorW(v, rover) # calculate motor speed (used in next line)
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor']) # calculate torque (used for efficiency info)
    
    # Determine efficiency for each time/velocity
    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau'].ravel() # change to 1D array
    effcy = rover['wheel_assembly']['motor']['effcy'].ravel()
    effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic') # fit the cubic spline
    effcy_dat = effcy_fun(tau)
    
    # Calculate battery power for each time/velocity
    P_batt = P/effcy_dat
    
    # integrate to calculate energy    
    E_motor = simpson(P_batt, t)
   # E_motor = np.trapz(P_batt,t)
    E = 6*E_motor # 6 wheels, each with a dedicated motor
    
    return E


def simulate_rover(rover,planet,experiment,end_event):
    """
    Inputs:     rover:  dict              Data structure specifying rover 
                                          parameters
               planet:  dict              Data dictionary specifying planetary 
                                          parameters
           experiment:  dict              Data dictionary specifying experiment 
                                          definition
            end_event:  dict              Data dictionary containing the 
                                          conditions necessary and sufficient 
                                          to terminate simulation of rover 
                                          dynamics                 
    
    Outputs:    rover:  dict              Updated rover structure including 
                                          telemetry information
    """
    # Check that the rover input is a dict
    if type(rover) != dict:
        raise Exception('rover input must be a dict')
    
    # Check that the planet input is a dict
    if type(planet) != dict:
        raise Exception('planet input must be a dict')
    
    # Check that the experiment input is a dict
    if type(experiment) != dict:
        raise Exception('experiment input must be a dict')
        
    # Check that the end_event input is a dict
    if type(end_event) != dict:
        raise Exception('end_event input must be a dict')
    
    # Main Code
    fun = lambda t,y: rover_dynamics(t, y, rover, planet, experiment) # differential equation
    t_span = experiment['time_range'] # time span
    y0 = experiment['initial_conditions'].ravel() # initial conditions
    events = end_of_mission_event(end_event) # stopping criteria
    sol = solve_ivp(fun, t_span, y0, method = 'BDF', events=events) #t_eval=(np.linspace(0, 3000, 1000)))  # need a stiff solver like BDF
    
    # extract necessary data
    v_max = max(sol.y[0,:])
    v_avg = mean(sol.y[0,:])
    P = mechpower(sol.y[0,:], rover)
    E = battenergy(sol.t,sol.y[0,:],rover)
    
    # Add telemetry info to rover dict
    telemetry = {'Time' : sol.t,
                 'completion_time' : sol.t[-1],
                 'velocity' : sol.y[0,:],
                 'position' : sol.y[1,:],
                 'distance_traveled' : sol.y[1,-1],  # Matlab version had an integration of velocity over time, but if velocity must be positive (defined by min_velocity), then the final position is the distance traveled
                 'max_velocity' : v_max,
                 'average_velocity' : v_avg,
                 'power' : P,
                 'battery_energy' : E,
                 'energy_per_distance' : E/sol.y[1,-1]}
    
    rover['telemetry'] = telemetry
    return rover

def end_of_mission_event(end_event):
    """
    Defines an event that terminates the mission simulation. Mission is over
    when rover reaches a certain distance, has moved for a maximum simulation 
    time or has reached a minimum velocity.            
    """
    
    mission_distance = end_event['max_distance']
    mission_max_time = end_event['max_time']
    mission_min_velocity = end_event['min_velocity']
    
    # Assume that y[1] is the distance traveled
    distance_left = lambda t,y: mission_distance - y[1]
    distance_left.terminal = True
    
    time_left = lambda t,y: mission_max_time - t
    time_left.terminal = True
    time_left.direction = -1
    
    velocity_threshold = lambda t,y: y[0] - mission_min_velocity;
    velocity_threshold.terminal = True
    
    # terminal indicates whether any of the conditions can lead to the
    # termination of the ODE solver. In this case all conditions can terminate
    # the simulation independently.
    
    # direction indicates whether the direction along which the different
    # conditions is reached matter or does not matter. In this case, only
    # the direction in which the velocity treshold is arrived at matters
    # (negative)
    
    events = [distance_left, time_left, velocity_threshold]
    
    return events




#%% Testing (delete after use!!!)

# from define_rovers import *
# from define_experiment import * 

# rover, planet, end_event = define_rover_4()
# rover, planet = define_rover_3()
# experiment = experiment1()

# t = [0, 10,20,30]
# y = np.array([[0.33,0],[0.30,100],[0.25,500],[0.10,750]])
# for i in range(len(t)):
#     dydt = rover_dynamics(t[i], y[i], rover, planet, experiment)
#     print(dydt)

# v = np.array([0,0.1,0.2,0.3,0.33])
# P = mechpower(v, rover)
# print(P)

# t = np.array(range(7), dtype = float)
# v = np.array([.33,.32,.33,.2,.2,.25,.28], dtype = float)
# E = battenergy(t,v,rover)
# print(E)

# sol = simulate_rover(rover,planet,experiment,end_event)