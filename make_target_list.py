from photo_files import photo_files
from apex import apex
#from apex3 import *
from photo_meta import photo_meta
from m1_photo import m1_photo
from m1_targets import m1_targets
from m1_show import m1_show
import matplotlib.pyplot as pl
import numpy as np


F = photo_files()
F.set_default_file_names('../2019Oct10Data/','M73','101019','SIMPLE',4)
A = apex(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=True)
META = photo_meta(F.META_FILE_1,F.META_FILE_2)
M = m1_photo(F.DATA_FILE_NAME,META,raw_read=True) 
M.rotate_cloud(A)
M.ring_id()
M.find_segments(print_it=True)
M.write_target_file('M73_TARGET_FILE.csv')
