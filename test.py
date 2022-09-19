# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 23:51:14 2022

@author: namqv_7yst908
"""

import unittest
import subfunctions as sf

class Test_get_mass():
    
    def test_value_checking_true(self):
        rover, planet = sf.create_dictionary()
        
        mass = 0
        mass += rover['wheel_assembly']['wheel']['mass'] * 6 #There are 6 wheels
        mass += rover['wheel_assembly']['speed_reducer']['mass']
        mass += rover['wheel_assembly']['motor']['mass']
        mass += rover['chassis']['mass']
        mass += rover['power_subsys']['mass']
        mass += rover['science_payload']['mass']
        
        self.assertEqual(sf.get_mass(rover),mass)
    
    def test_value_checking_false(self):
        self.assertRaises(Exception, sf.get_mass(9))
        
        

        

if __name__ == '__main__':
    unittest.main()
    print("All test passed successfully")


