import sys
import numpy as np
import math
import matplotlib.pyplot as pl

from photo_files import photo_files
from apex import apex
from m1_photo import m1_photo
from m1_show import m1_show
from m1_targets import m1_targets
from photo_meta import photo_meta
from segment_photo import segment_photo
from segment_show import segment_show

def main(argv):
    # initialize parameters ***************************************************

    # get file information from command line arguments
    if len(argv)<7:
        print('segment_command usage: path, map, session, model_name, model_npars, ring, segment')
        sys.exit(-1)
    else:
        # set commonly changed options for the program
        PATH = argv[0]              # path for the data files
        MAP = argv[1]               # name of the map (format is 'MXX', where XX is the map number)
        SESSION = argv[2]           # session name (format is 'ddmmyy')
        MODEL_NAME = argv[3]        # name of model for fit: 'SIMPLE' or 'ZERNIKE'
        MODEL_NPARS = eval(argv[4]) # number of parameters in model 
        # this command sets up a set of nominal filenames for all our files
        F = photo_files()
        F.set_default_file_names(PATH,MAP,SESSION,MODEL_NAME,MODEL_NPARS)

        theRing = eval(argv[5])
        theSegment = eval(argv[6])

    # other program options set in script                            
    # set rings for analysis and number of segments per ring              
    SURFACE_FIT_RINGS = [1,2,3]
    SURFACE_RESIDUALS_RINGS = [1,2,3,4,5]
    ACTUATOR_FIT_RINGS = [1,2,3,4,5]
    N_SEGMENTS = [12,24,48,48,48]

    # plotting parameters                              
    PLOT_IT = True     # plot the maps
    PRINT_IT = False   # make png files of figures
    IMAGE_SCALE_BEFORE = [-3.,3.]
    IMAGE_SCALE = [-0.25, 0.25] # sets vertical limits on images (mm)
    IMAGE_EXTENT = [-26., 26., -26., 26.]# sets x-y limits for images (m) Nominal is 26; use 17 for 3 rings
    CIRCLE_SCALE_BEFORE = 0.005
    CIRCLE_SCALE = 0.001    # scale parameter for circle plots (nominal is 0.001)                                                    
    CIRCLE_SCALE_BAR = 1.0  # set size of scale bar in circle plot (nominal is 0.2 mm)                                                        
    CIRCLE_EXTENT = [-26000., 26000.,-26000., 26000.] # sets x-y limits for circle plot (in mm); use 17000 for 3 rings
    DATA_SIGMA = 0.075      # sets a reference for histrogram drawn on residual distribution (mm)

    # target file
    TARGET_FILE = 'M11_NEW_TARGET_LIST.csv'
#    TARGET_FILE = 'M1_TARGET_LIST.csv'   # this is the file with the targets listed
#    TARGET_FILE = 'M06_TARGET_FILE_TEMP_2.csv'   # this is the file with the targets listed

    # command for actuator command file (nominal is 0) 
    ACTUATOR_COMMAND = 0
    # segment deformation model (0=FEM; 1=HOL); Must use HOL option for 5 rings
    SEGMENT_DEFORMATION_MODEL = 1

    # set up files and strings for titles
    MODEL_STRING = '%s_%02d'%(MODEL_NAME,MODEL_NPARS)

    # create MapName for labels                                                          
    MAP_NAME = MAP+'_'+SESSION

    # start processing *******************************************************

    print('MAP: %s  **************************'%(MAP_NAME))
    print('METAFILES: %s %s'%(F.META_FILE_1,F.META_FILE_2))
    META = photo_meta(F.META_FILE_1,F.META_FILE_2)
    print('METAFILE DATE: %s'%(META.get_date()))

    print('USING REDUCED POINT CLOUD DATA: %s'%(F.RED_DATA_FILE_NAME))
    M = m1_photo(F.RED_DATA_FILE_NAME,META,raw_read=False)

    # load the target list and associate targets with cloud points
    T = m1_targets(TARGET_FILE)
    M.find_target_ids(T)
    M.elim_targets(T)
    M.find_segments()
    T.check_targets()
    
    # set the model to be fit to the cloud from command line args
    M.set_model(ModelType=MODEL_NAME,npar=MODEL_NPARS)

    M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
    # last pass gives residuals to all rings
    M.make_residuals_vector(res_rings=SURFACE_RESIDUALS_RINGS)

    seg = segment_photo(M,theRing,theSegment,META,SEGMENT_DEFORMATION_MODEL)

    pl.ion()
    ss = segment_show(seg,figure=99)
    ss.make_segment_plot()
    seg.print_screw_turns(1270)
        

main(sys.argv[1:])
