import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from scipy.integrate import simpson
from statistics import mean

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
                    planet:  dict            Data dictionary specifying planetary 
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
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
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
    Inputs:               v:  np.ndarray     
                      rover:  dictionary            
    
    Outputs:              omega:  np.ndarray    angular velocity of motor
    """
    #Type Checking
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise Exception('v must be a scalar or a vector of type np.ndarray')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float)
    elif len(np.shape(v)) != 1:
        raise Exception('v must be a scalar or a vector of type np.ndarray')
    
    if type(rover) != dict:
        raise Exception('rover input must be a dict')
        
        
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    radius = rover['wheel_assembly']['wheel']['radius']
    
    omega = v*Ng/radius
    
    return omega


def rover_dynamics(t, y, rover, planet, experiment):
    """
    Inputs:         t:  scalar            time[s]
                    y:  np.ndarray        vector and position array
                rover:  dictionary        
               planet:  dictionary              
           experiment:  dictionary              
    
    Outputs:     dydt:  np.ndarray       derivative of y
    """
    #Type checking
    if (type(t) != int) and (type(t) != float) and (not isinstance(t, np.ndarray)) and (not isinstance(t,np.float64)):
        raise Exception('t should be a scalar')
    elif isinstance(t, np.ndarray):
        if len(t) == 1:
            t = float(t)
        else:
            raise Exception('T should be a scalar')
    
    if (not isinstance(y, np.ndarray)) or (len(y) != 2):
        raise Exception('Y should be a 2x1 np.ndarray')
    elif isinstance(y[0], np.ndarray):
        y = np.array([float(y[0]),float(y[1])]) 
        
    if type(rover) != dict:
        raise Exception('Rover should be a dictionary')

    if type(planet) != dict:
        raise Exception('Planet should be a dictionary')
    
    if type(experiment) != dict:
        raise Exception('Experiment should be a dictionary')
    
    
    omega = motorW(float(y[0]), rover)   
    alpha_fun = interp1d(experiment['alpha_dist'].ravel(), experiment['alpha_deg'].ravel(), kind = 'cubic', fill_value="extrapolate")
    
    angle = float(alpha_fun(float(y[1])))
    
    Force = F_net(omega, angle, rover, planet, experiment['Crr'])
    
    mass = get_mass(rover)
    acc = float(Force /mass)
    
    
    dy = np.array([acc, float(y[0])], dtype = float)
    
    return dy


def mechpower(v, rover):
    """
    Inputs:               v:  np.ndarray
                      rover:  dictionary          
    
    Outputs:            mechanicalPower:  np.ndarray
    """
    #Type checking
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise Exception('v must be a scalar or a vector of np.ndarray')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float)
    elif len(np.shape(v)) != 1:
        raise Exception('v must be a scalar or a vector of np.ndarray')
    
    if type(rover) != dict:
        raise Exception('rover should be a dictionary')
    
    omega = motorW(v, rover)  
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor']) 
    
    mechanicalPower = omega * tau
    
    return mechanicalPower


def battenergy(t,v,rover):
    """
    Inputs:               t:  np.ndarray     time
                          v:  np.ndarray     velocity
                      rover:  dictionary            
    
    Outputs:         Energy:  scalar         total energy output by battery
    """
    #Type Checking
    if (not isinstance(t, np.ndarray)):
        raise Exception('t must be a scalar or vector of np.ndarray')
    elif len(np.shape(t)) != 1:
        raise Exception('t must be a scalar or vector of np.ndarray')
    if (not isinstance(v, np.ndarray)):
        raise Exception('v must be a scalar or vector of np.ndarray')
    elif len(np.shape(v)) != 1:
        raise Exception('v must be a scalar or vector of np.ndarray') 
    if len(t) != len(v):
        raise Exception('t and v should be the same size')
    
    #efficiency function
    effcy_fun = interp1d(rover['wheel_assembly']['motor']['effcy_tau'].ravel(), rover['wheel_assembly']['motor']['effcy'].ravel(), kind = 'cubic')
    
    #calculating battery power
    batteryPower = mechpower(v, rover)/effcy_fun(tau_dcmotor(motorW(v, rover), rover['wheel_assembly']['motor']))
    
    #calculaing energy output for each motor
    motorEnergy = simpson(batteryPower, t)
    
    #Energy output of battery
    Energy = motorEnergy * 6
    
    return Energy


def simulate_rover(rover,planet,experiment,end_event):
    """
    Inputs:     rover:  dictionary              
               planet:  dictionary              
           experiment:  dictionary              
            end_event:  dictionary                              
    
    Outputs:    rover:  dictionary              
    """
    #Type Checking
    if type(rover) != dict:
        raise Exception('rover should be a dictionary')
    if type(planet) != dict:
        raise Exception('planet should be a dictionary')
    if type(experiment) != dict:
        raise Exception('experiment should be a dictionary')
    if type(end_event) != dict:
        raise Exception('end_event should be a dictionary')
    
    #creating differential equation
    func = lambda t,y: rover_dynamics(t, y, rover, planet, experiment)
    
    #obtaining necessary variables
    time_span = experiment['time_range']
    y_ini = experiment['initial_conditions'].ravel()
    
    #solving differential equations
    sol = solve_ivp(func, time_span, y_ini, method = 'BDF', events=end_of_mission_event(end_event))
    
    #creating telemetry data
    telemetry = {'Time' : sol.t,
                 'completion_time' : sol.t[-1],
                 'velocity' : sol.y[0,:],
                 'position' : sol.y[1,:],
                 'distance_traveled' : sol.y[1,-1],
                 'max_velocity' : max(sol.y[0,:]),
                 'average_velocity' : mean(sol.y[0,:]),
                 'power' : mechpower(sol.y[0,:], rover),
                 'battery_energy' : battenergy(sol.t,sol.y[0,:],rover),
                 'energy_per_distance' : battenergy(sol.t,sol.y[0,:],rover)/sol.y[1,-1]}
    
    #adding telemetry data to rover
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