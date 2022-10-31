# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 09:17:49 2022

@author: namqv_7yst908
"""

from define_experiment import experiment1
from define_rover_planet import create_dictionary
from subfunctions import simulate_rover

planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover = create_dictionary()
experiment, end_event = experiment1()

rover_updated = simulate_rover(rover, planet, experiment, end_event)
print(rover_updated)