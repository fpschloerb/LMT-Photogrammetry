"""Module with photo_files class for managing file names within system

Clases:  photo_files
Uses:    numpy, m1_photo, m1_targets, photo_meta
Author:  FPS
Date:    November 6, 2016
Changes: December 2017 - new naming conventions and methods
         Filenames may now be loaded without following any naming convention
Parameters:
FOR AUTOMATICALLY NAMING ACCORDING TO A CONVENTION
apex_root            common root for apex measurement files
data_root            common root for data measurement files
path                 path for all photogrammetry files
map                  map identifier (e.g. M04)
session              photogrammetry session identifier (DDMMYY)
model                model identifier (SIMPLE, TILE, ZERNIKE)
npar                 number of model parameters
q                    additional identifier once used in maps
THE ACTUAL FILE NAMES CONSTRUCTED OR INDIVIDUALLY LOADED
APEX_FILE_NAME       the actual apex file name
RED_APEX_FILE_NAME   the file for reduced apex data (appends _REDUCED)
DATA_FILE_NAME       the actual data file name
RED_DATA_FILE_NAME   the file for reduced data (appends MODEL_NPAR_REDUCED)
META_FILE_1          the metadata file for beginning of map
META_FILE_2          the metadata file for end of map
COMMAND_FILE_NAME    the filename for actuator commands
"""

import numpy as np
from m1_photo import m1_photo
from m1_targets import m1_targets
from photo_meta import photo_meta

class photo_files():
    """photo_files manages filenames for photogrammetry"""
    def __init__(self,apex='APEX_FOTO_GLOBAL_',data='DAT_FOTO_GLOBAL_',q=''):
        self.apex_root = apex
        self.data_root = data
        self.q = q
    def set_default_file_names(self,path,map,session,model,npar):
        self.set_path(path)
        self.set_map(map)
        self.set_session(session)
        self.set_model_string(model,npar)
        self.set_apex_file(self.path+self.apex_root+self.map+'_'+self.session+self.q+'.txt')
        self.set_reduced_apex_file()
        self.set_data_file(self.path+self.data_root+self.map+'_'+self.session+self.q+'.txt')
        self.set_reduced_data_file()
        self.set_meta_1_file(self.path+self.map+'_META_1.csv')
        self.set_meta_2_file(self.path+self.map+'_META_2.csv')
        self.set_command_file(self.path+self.map+'_'+self.session+'_'+self.model_string+'_COMMAND_TOTAL.csv')
    def set_path(self,path):
        self.path = path
    def set_map(self,map):
        self.map = map
    def set_session(self,session):
        self.session = session
    def set_model_string(self,model,npar):
        self.model_string = '%s_%02d'%(model,npar)
    def set_apex_file(self,file):
        self.APEX_FILE_NAME = file
    def set_reduced_apex_file(self):
        root = self.APEX_FILE_NAME.split('.')
        n = len(root)
        if n==2:
            self.RED_APEX_FILE_NAME = root[0]+'_REDUCED.'+root[1]
        else:
            ss = ''
            for i in range(n-2):
                ss = ss+'.'
            self.RED_APEX_FILE_NAME = ss+root[-2]+'_REDUCED.'+root[-1]

    def set_data_file(self,file):
        self.DATA_FILE_NAME = file
    def set_reduced_data_file(self):
        root = self.DATA_FILE_NAME.split('.')
        n = len(root)
        if n==2:
            self.RED_DATA_FILE_NAME = root[0]+'_'+self.model_string+'_REDUCED.'+root[1]
        else:
            ss = ''
            for i in range(n-2):
                ss = ss+'.'
            self.RED_DATA_FILE_NAME = ss+root[-2]+'_'+self.model_string+'_REDUCED.'+root[-1]
    def set_meta_1_file(self,file):
        self.META_FILE_1 = file
    def set_meta_2_file(self,file):
        self.META_FILE_2 = file
    def set_command_file(self,file):
        self.COMMAND_FILE_NAME = file
    def print_file_names(self):
        print(self.APEX_FILE_NAME)
        print(self.RED_APEX_FILE_NAME)
        print(self.DATA_FILE_NAME)
        print(self.RED_DATA_FILE_NAME)
        print(self.META_FILE_1)
        print(self.META_FILE_2)
        print(self.COMMAND_FILE_NAME)
    def parse_data_file_name(self):
        # parses the data file name to get map number and session
        # this only works if we follow convention - NO CHECKS!
        r = self.DATA_FILE_NAME.split('_')
        return(self.DATA_FILE_NAME.split('DAT')[0], r[-2], r[-1].split('.')[0])
