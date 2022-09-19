import numpy as np
import matplotlib.pyplot as plt
import subfunctions as sf
import A2
from mpl_toolkits.mplot3d import Axes3D


### VARIABLES FROM SUBFUNCTIONS
planet = {
    'gravity' : -3.72 #acceleration due of gravity in m/s^2
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


### PRELIMINARIES
Crr_array = np.linspace(0.01, .4, 25)
slope_array_deg = np.linspace(-10, 35, 25)

CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)
VMAX = np.zeros(np.shape(CRR), dtype = float)

'''
### FUNCTIONS USED WITHIN THE CODE
def secant(fun, ini_guess, err_max = 1e-5, iter_max = 10000):
    
    isDone = False
    numIter = 0
    err_est = np.nan
    root = np.nan
    
    x0 = ini_guess
    x1 = ini_guess + .001
    
    if fun(x0) * fun(x1) == 0 :
        isDone = True
        err_est = 0.0
        exitFlag = 1
        if abs(fun(x0)) < abs(fun(x1)):
            root = x0
        else:
            root = x1

    if fun(x0) == np.nan:
        isDone = True
        exitFlag = -2
    
    while not isDone:
        
        numIter += 1
        
        x2 = x1 - (fun(x0) * (x1 - x0)) / (fun(x1) - fun(x0))
        
        if fun(x2) == np.nan:
            isDone = True
            exitFlag = -2
            break
        
        if fun(x0) * fun(x2) == 0:
            isDone = True
            exitFlag = 1
            if abs(fun(x0)) <= abs(fun(x2)):
                root = x0
            else:
                root = x2
        
        err_est = 100 * abs((x2 - x0) / x2)
        
        if numIter > iter_max:
            isDone = True
            exitFlag = 0
            root = x2
            break
        elif err_est < err_max: 
            isDone = True
            exitFlag = 1
            root = x2
            break
        else:
            x0 = x1
            x1 = x2
        
    return root, err_est, numIter, exitFlag
'''

### CODE
N = np.shape(CRR)[0]

def fun1(x):
    drive = sf.F_drive(np.array([x]),rover)
    grav = sf.F_gravity(np.array([slope_sample]),rover,planet)
    roll = sf.F_rolling(np.array([x]),np.array([slope_sample]),rover,planet,(Crr_sample))
    
    return float(drive + grav + roll ) 

for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPE[i,j])
        VMAX[i,j] = (A2.secant(fun1,1) * .3)
        
figure = plt.figure() 
ax = Axes3D(figure, elev = 20, azim = 25) 
ax.plot_surface(CRR, SLOPE, VMAX) 