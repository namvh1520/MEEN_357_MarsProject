import numpy as np
import subfunctions as sf
import A2 as A2
import matplotlib.pyplot as plt

Crr = .2

slope_array_deg1 = np.linspace(-10,35,25)


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



def fun1(x):
    drive = sf.F_drive(np.array([x]),rover)
    grav = sf.F_gravity(np.array([slope_array_deg]),rover,planet)
    roll = sf.F_rolling(np.array([x]),np.array([slope_array_deg]),rover,planet,Crr)
    mot = sf.tau_dcmotor(np.array([x]), motor)
    
    return float(drive + grav + roll )

def fun(x):
    return float(sf.F_net(np.array([x]),np.array([slope_array_deg]),rover,planet,Crr))

temp = []
hope =[]

for i in slope_array_deg1:
    slope_array_deg = i
    hope.append(i)
    temp.append(A2.secant(fun1,3) * .3)


v_max = np.array(temp)

plt.plot(slope_array_deg1,v_max)
plt.xlabel('Slope (Degrees)')
plt.ylabel('Velocity (M/s)')
plt.title('Velocity vs Slope')
