"""Module with photo_meta class which manages the metadata associated with photogrammetry

Classes:  photo_meta
Uses:     numpy, csv, TempSens
Author:   FPS
Date:     November 6, 2016
Changes:  Now maintaining versions for python 2 and python 3
          THIS IS THE PYTHON 3 VERSION
          December 14, 2017 updated default number of actuators to 720
"""

import numpy as np
import csv
from TempSens import TempSens

# this is the python 3 version

class photo_meta():
    """metadata for photogrammetry map"""
    def __init__(self,filename_1,filename_2,nactuators=720):
        f = open(filename_1)
        r = csv.reader(f)
        row1_1 = r.__next__()
        row2_1 = r.__next__()
        f.close()

        f = open(filename_2)
        r = csv.reader(f)
        row1_2 = r.__next__()
        row2_2 = r.__next__()
        f.close()

        # start date time string
        self.date = row2_1[0]

        # read in the actuators
        self.actuators = np.zeros((720,3))
        act_index = 6
        idx = 0
        for i in range(act_index,act_index+nactuators):
            self.actuators[idx,0] = eval(row2_1[i])
            self.actuators[idx,1] = eval(row2_2[i])
            self.actuators[idx,2] = (self.actuators[idx,0]+self.actuators[idx,1])/2.
            idx = idx + 1

        if nactuators == 336:
            # initialize r4 and r5 to zero for now
            for i in range(336,720):
                self.actuators[i,0] = 0.
                self.actuators[i,1] = 0.
                self.actuators[i,2] = 0.

        # read in the temperature sensors
        self.temperatures = np.zeros((64,3))
        temp_index = 6+nactuators
        idx = 0
        for i in range(temp_index,temp_index+64):
            self.temperatures[idx,0] = eval(row2_1[i])
            self.temperatures[idx,1] = eval(row2_2[i])
            self.temperatures[idx,2] = (self.temperatures[idx,0]+self.temperatures[idx,1])/2.
            idx = idx + 1

        self.TempSens = TempSens(self.temperatures[:,2])

        # read in the elevation
        elevation_index = 6+nactuators+64+5
        self.elevation = (eval(row2_1[elevation_index]) + eval(row2_2[elevation_index]))/2.
        print('Elevation %f %f %f'%(self.elevation,eval(row2_1[elevation_index]),eval(row2_2[elevation_index])))
                         
    def get_actuators(self,idx=2):
        return(self.actuators[:,idx])

    def get_temperatures(self,idx=2):
        return(self.temperatures[:,idx])

    def get_date(self):
        return self.date

    def get_elevation(self):
        return self.elevation

    def print_csv_dish_temperatures(self):
        print('%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f'%(self.TempSens.dish,
                                self.TempSens.dish_left,
                                self.TempSens.dish_right,
                                self.TempSens.dish_top,
                                self.TempSens.dish_bottom,
                                self.TempSens.dish_front,
                                self.TempSens.dish_back,
                                self.TempSens.dish_inner,
                                self.TempSens.dish_outer,
                                self.TempSens.upper_alidade,
                                self.TempSens.alidade_base,
                                self.TempSens.ballast))
        print(' ')
