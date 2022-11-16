# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:52:21 2022

@author: namqv_7yst908
"""
import numpy as np
import matplotlib.pyplot as plt
from define_edl_system import *
from subfunctions_EDL import *
from define_planet import *
from define_mission_events import *
from redefine_edl_system import *


def simulate_edl_parachutestudy(edl_system, planet, mission_events, tmax, ITER_INFO):
    events = edl_events(edl_system, mission_events)
        
    tspan = (0, tmax)
    
    y0 = np.array([edl_system['velocity'], 
                   edl_system['altitude'], 
                   edl_system['rocket']['initial_fuel_mass'] * edl_system['num_rockets'],
                   0,
                   0,
                   0,
                   0])
    if ITER_INFO:
        print('Commencing simulation run...\n')
    
    T = np.array([])
    Y = np.array([[], [], [], [], [], [], []])
    TERMINATE_SIM = False
    while not(TERMINATE_SIM):
        
        fun = lambda t, y: edl_dynamics(t, y, edl_system, planet)
        sol = solve_ivp(fun, tspan, y0, method='DOP853', events=events, max_step=0.1)
        t_part = sol.t
        Y_part = sol.y
        TE = sol.t_events
        YE = sol.y_events

        [edl_system, y0, TERMINATE_SIM] = update_edl_state(edl_system, TE, YE, Y_part, ITER_INFO)
        
        tspan = (t_part[-1], tmax)

        if TERMINATE_SIM:
            time = TE[8][0]
            altitude = YE[8][0, 1]
            speed = YE[8][0, 0]
            rover_rel_pos = YE[8][0, 6]
            rover_rel_vel = YE[8][0, 5]
            rover_touchdown_speed = speed + rover_rel_vel
            
            if altitude >= edl_system["sky_crane"]["danger_altitude"] and abs(
                rover_touchdown_speed
            ) <= abs(edl_system["sky_crane"]["danger_speed"]):
                success = 1
            else:
                success = 0
            
        T = np.append(T, t_part)
        Y = np.hstack((Y, Y_part))
        
        if tspan[0] >= tspan[1]:
            TERMINATE_SIM = True
    
    return T, Y, edl_system, time, rover_touchdown_speed, success

edl_system = define_edl_system_1()
edl_system = redefine_edl_system(edl_system)
mars = define_planet()
mission_events = define_mission_events()

parachute_size = np.arange(14, 19.1, 0.5)
tmax = 2000

parachutediameter = []
simulatedtime = []
roverlandingsuccess = []
roverspeed = []

for i in parachute_size:
    parachutediameter.append(i)
    edl_system = define_edl_system_1()
    edl_system = redefine_edl_system(edl_system)
    mars = define_planet()
    mission_events = define_mission_events()
    edl_system['parachute']['diameter'] = i
    [t, Y, edl_system, time, rover_touchdown_speed, success] = simulate_edl_parachutestudy(edl_system, mars, mission_events, tmax, False)
    
    simulatedtime.append(time)
    roverlandingsuccess.append(success)
    roverspeed.append(rover_touchdown_speed)
    
fig, axs = plt.subplots(3) #plotting in a subplot
plt.subplots_adjust(hspace=1.2) #increasing distance between the three plots
axs[0].plot(parachutediameter, simulatedtime)
axs[0].set_title('Simulated Time vs Diameter')
axs[0].set_ylabel('Time (s)')
axs[0].set_xlabel('Diameter (m)')

axs[1].plot(parachutediameter, roverspeed)
axs[1].set_title('Speed vs Diameter')
axs[1].set_ylabel('Speed (m/s)')
axs[1].set_xlabel('Diameter (m)')

axs[2].plot(parachutediameter, roverlandingsuccess)
axs[2].set_title('Success vs Diameter')
axs[2].set_ylabel('Success')
axs[2].set_xlabel('Diameter (m)')
    


