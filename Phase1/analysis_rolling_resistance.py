import matplotlib.pyplot as plt
import numpy as np
import subfunctions as sf
import A2 as A2

Crr = np.linspace(.1,.4,25)

slope_array_deg1 = np.array([0])

planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover = sf.create_dictionary()


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
