import numpy as np
import matplotlib.pyplot as plt
import subfunctions as sf
from mpl_toolkits.mplot3d import Axes3D

planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover = sf.create_dictionary()


### PRELIMINARIES
Crr_array = np.linspace(0.01, .4, 25)
slope_array_deg = np.linspace(-10, 35, 25)

CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)
VMAX = np.zeros(np.shape(CRR), dtype = float)


### CODE
N = np.shape(CRR)[0]

def fun1(x):
    drive = sf.F_drive(np.array([x]),rover)
    grav = sf.F_gravity(np.array([slope_sample]),rover,planet)
    roll = sf.F_rolling(np.array([x]),np.array([slope_sample]),rover,planet,(Crr_sample))
    
    return float(drive + grav + roll ) 

def fun2(x):
    return sf.F_net(np.array([x]), np.array([slope_sample]), rover, planet, Crr_sample)

for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPE[i,j])
        VMAX[i,j] = (sf.secant(fun1,1) * wheel['radius'])
        
figure = plt.figure() 
ax = Axes3D(figure, elev = 15, azim = 30) 
ax.plot_surface(CRR, SLOPE, VMAX) 
plt.title('Speed vs Slope and Rolling Resistance Coefficient')
plt.ylabel('Slope (Degrees)')
ax.set_zlabel('Speed (M/s)')
plt.xlabel('CRR')
ax.set_xticks([.0,.1,.2,.3,.4])