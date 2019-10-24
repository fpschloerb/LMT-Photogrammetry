"""Module with m1_targets class for managing lists of photogrammetry targets

Classes: m1_targets
Uses:    math,numpy,matplotlb.pyplot,csv
Date:    November 6, 2016
Changes: December 2017 - improved handling of target lists

Parameters:
filename     __init__ file name for target data
n            __init__ number of targets
ring[n]      __init__ ring id's
seg[n]       __init__ segment id's
x[n]         __init__ x location array
y[n]         __init__ y location array
flag[n]      __init__ indicates whether valid
point_id[n+1] __init__ data point id number assocated with target
found_target[n+1] __init__ target has been found (ideally 1)
"""

import numpy as np
import math
import matplotlib.pyplot as pl
import csv

class m1_targets():
    """m1_targets manages the photogrammetry target list"""
    def __init__(self,filename):
        # read in the target list from csv file
        self.filename = filename
        f = open(self.filename)
        target_list = csv.reader(f)
        rr = []
        ss = []
        xx = []
        yy = []
        ee = []
        for row in target_list:
            rr.append(eval(row[0]))
            ss.append(eval(row[1]))
            xx.append(eval(row[2]))
            yy.append(eval(row[3]))
            ee.append(eval(row[4]))
        f.close()    

        self.ring = np.array(rr,dtype='int')
        self.seg = np.array(ss,dtype='int')
        self.x = np.array(xx)
        self.y = np.array(yy)
        self.flag = np.array(ee,dtype='int')
        self.n = len(self.ring)
        self.point_id = np.zeros(self.n+1,dtype='int')
        self.found_target = np.zeros(self.n+1,dtype='int')

    def set_flag(self, target_id):
        self.flag[target_id] = 0

    def clear_flag(self, target_id):
        self.flag[target_id] = 1

    def get_flag(self, target_id):
        return(self.flag[target_id])

    def get_targets_not_found(self):
        return(self.found_target[-1])

    def check_targets(self):
        bad = 0
        bad_multiple = 0
        for i in range(self.n):
            if self.found_target[i]>1:
                print('Multiple points found for Target %d: %d %d %d'
                      %(i,self.ring[i],self.seg[i],self.found_target[i]))
                bad_multiple = bad_multiple + 1
            elif self.found_target[i] < 1:
                print('Target %d not found: %d %d %d'
                      %(i,self.ring[i],self.seg[i],self.found_target[i]))
                bad = bad + 1
        if bad_multiple > 0:
            print('%d TARGETS WITH MULTIPLE POINTS'%(bad_multiple))
        if bad == 0:
            print('ALL TARGETS FOUND')
            return(1)
        else:
            # we set the number of targets not found in the extra array element
            self.found_target[-1] = bad
            print('%d MISSING TARGETS'%(bad))
            return(0)

    def set_point_id(self,target_id,point_id):
        self.point_id[target_id] = point_id
        self.found_target[target_id] += 1

    def get_point_id(self,target_id):
        return(self.point_id[target_id])

    def get_found_target(self,target_id):
        return(self.found_target[target_id])

