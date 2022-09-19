import matplotlib.pyplot as plt
import numpy as np
import subfunctions as sf
import A2 as A2

Crr = .2

slope_array_deg1 = np.linspace(0.01,0.4,25)


rover, planet = sf.create_dictionary()



def fun1(x):
    drive = sf.F_drive(np.array([x]),rover)
    grav = sf.F_gravity(np.array([slope_array_deg]),rover,planet)
    roll = sf.F_rolling(np.array([x]),np.array([slope_array_deg]),rover,planet,Crr)
    
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
plt.xlabel('Rolling Resistance')
plt.ylabel('Velocity (M/s)')
