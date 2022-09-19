import matplotlib.pyplot as plt
import numpy as np
import subfunctions as sf
import A2 as A2

Crr = np.linspace(.1,.4,25)

slope_array_deg1 = np.array([0])


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
    grav = sf.F_gravity(np.array([0]),rover,planet)
    roll = sf.F_rolling(np.array([x]),np.array([0]),rover,planet,float(np.array([Crr1])))
    
    return float(drive + grav + roll )

temp = []
hope =[]

for i in Crr:
    Crr1 = i
    hope.append(i)
    temp.append(A2.secant(fun1,3) * .3)


v_max = np.array(temp)

plt.plot(Crr,v_max)
plt.xlabel('Rolling Resistance')
plt.ylabel('Velocity (M/s)')
plt.title('Velocity vs Rolling Resistance')
