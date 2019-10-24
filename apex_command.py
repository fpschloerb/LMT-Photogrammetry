import sys
import numpy as np
import math

from apex3 import apex3
from photo_files import photo_files
def main(argv):
    # get file information from command line arguments
    if len(argv)<2:
        print('usage: path, map, session')
        exit(-1)
    else:
        F = photo_files()
        PATH = argv[0]
        MAP = argv[1]
        SESSION = argv[2]
        F.set_default_file_names(PATH,MAP,SESSION,'DUMMY',0)
        
    print('SOLVING WITH RAW APEX FILE: %s'%(F.APEX_FILE_NAME))
    A = apex3(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=True)

main(sys.argv[1:])
