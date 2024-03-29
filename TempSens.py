"""Module with useful classes and functions for the pointing analysis.

Classes:   TempSens
Functions: 
Uses:      numpy,math
Author:    FPS
Date:      May 5,2014
Changes:
"""
import numpy as np
import math

class TempSens():
    """Class to manage the temperature sensor data."""
    def __init__(self,array):
        """__init__ loads array with temperature sensor data, checks it, and calculates special values."""
        self.tempsens = array
        self.valid = np.zeros(len(array))
        self.validate()
        self.calculate()

    def calculate(self):
        """Calculates special values from the array of data."""
        # mean temperature
        self.mean_temp = self.process_list(range(64))
        # RTopUpper,RTopLower,LBotUpper,LBotLower,RBotUpper,RBotLower,LTopUpper,LTopLower
        self.legs = self.process_list((13,14,15,16,17,18,19,20))
        self.legs_left = self.process_list((15,16,19,20))
        self.legs_right = self.process_list((13,14,17,18))
        self.legs_top = self.process_list((13,14,19,20))
        self.legs_bottom = self.process_list((15,16,17,18))
        # Dish Az: 0,60,120,180,240,300 outer,middle,inner top then bottom
        self.dish = self.process_list(
            (10,9,31,30,51,50,
             12,11,33,32,63,62,
             2,21,23,22,43,42,
             4,3,25,24,45,44,
             6,5,27,26,47,46,
             8,7,29,28,49,48)
            )
        self.dish_left = self.process_list(
            (12,11,33,32,63,62,
             2,21,23,22,43,42)
            )
        self.dish_right = self.process_list(
            (6,5,27,26,47,46,
             8,7,29,28,49,48)
            )
        self.dish_top = self.process_list(
            (8,7,29,28,49,48,
             10,9,31,30,51,50,
             12,11,33,32,63,62)
            )
        self.dish_bottom = self.process_list(
            (2,21,23,22,43,42,
             4,3,25,24,45,44,
             6,5,27,26,47,46)
            )
        self.dish_front = self.process_list(
            (10,31,51,
             12,33,63,
             2,23,43,
             4,25,45,
             6,27,47,
             8,29,49)
            )
        self.dish_back = self.process_list(
            (9,30,50,
             11,32,62,
             21,22,42,
             3,24,44,
             5,26,46,
             7,28,48)
            )
        self.dish_inner = self.process_list(
            (51,50,
             63,62,
             43,42,
             45,44,
             47,46,
             49,48)
            )
        self.dish_middle = self.process_list(
            (31,30,
             33,32,
             23,22,
             25,24,
             27,26,
             29,28)
            )
        self.dish_outer = self.process_list(
            (10,9,
             12,11,
             2,21,
             4,3,
             6,5,
             8,7)
            )

        # LBack RBack LFront RFront CFrontBase RBase LBase CFrontX
        self.upper_alidade = self.process_list((34,35,55,54,53,57,61,64))
        # LFront LBack RBack RFront
        self.alidade_base = self.process_list((36,37,38,39))
        # LBot LTop RTop RBot
        self.ballast = self.process_list((40,41,59,60))
        if self.valid[0]:
            self.subref = self.tempsens[0]
        else:
            self.subref = -1
    
                                          
    def validate(self):
        """Checks the values in the array to see if valid."""
        for i in range(len(self.tempsens)):
            if self.tempsens[i] < -100 or self.tempsens[i] > 310:
                self.valid[i] = 0
            else:
                self.valid[i] = 1
    
    def process_list(self,thelist):
        """Calculates the average of a list of sensors you provide."""
        temp = 0.0
        count = 0
        for i in thelist:
            index = i-1
            if self.valid[index] == 1:
                temp = temp + self.tempsens[index]
                count = count + 1
        if count >0:
            temp = temp/count
        else:
            temp = -1
        return temp
