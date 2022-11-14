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

for i in prachute_size:
    parachutediameter.append(i)
    edl_system['parachute']['diameter'] = i
    [t, Y, edl_system] = simulate_edl(edl_system, mars, mission_events, tmax, False)
    


