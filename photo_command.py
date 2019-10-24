import sys
import numpy as np
import math
import matplotlib.pyplot as pl

from photo_files import photo_files
from apex import *
#from apex3 import *
from m1_photo import m1_photo
from m1_show import m1_show
from m1_targets import m1_targets
from photo_meta import photo_meta
from segment_photo import segment_photo
from segment_show import segment_show

def main(argv):
    # initialize parameters ***************************************************

    # get file information from command line arguments
    if len(argv)<8:
        print('usage: path, map, session, model_name, model_npars,use_reduced_apex,use_reduced_data,multipass')
        exit(-1)
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

        # set switches for reduction options
        USE_REDUCED_APEX_FILE = eval(argv[5])# T: uses a previously derived solution for apex geometry
        USE_REDUCED_DATA_FILE = eval(argv[6])# T: uses a previously derived point cloud
        MULTIPASS_DATA_FILE = eval(argv[7])  # T: makes multiple pass/point elimination of point cloud points
        if eval(argv[8])==0:
            PARABOLA_CORRECTION = False
        else:
            PARABOLA_CORRECTION = True


    # other program options set in script                            
    # set rings for analysis and number of segments per ring              
    SURFACE_FIT_RINGS = [1,2,3,4,5]
    print('FIT ALL FIVE RINGS')
    SURFACE_RESIDUALS_RINGS = [1,2,3,4,5]
    ACTUATOR_FIT_RINGS = [1,2,3,4,5]
    N_SEGMENTS = [12,24,48,48,48]

    # plotting parameters                              
    PLOT_IT = True     # plot the maps
    PRINT_IT = False   # make png files of figures
    VIEW_SEGMENT_FITS = False  # pause after main fit and look at individual segment results
    PRINT_SCREW_TURNS = True # print the results of screw turn calculation
    IMAGE_SCALE_BEFORE = [-0.5,0.5]
    IMAGE_SCALE = [-0.25, 0.25] # sets vertical limits on images (mm)
    IMAGE_EXTENT = [-26., 26., -26., 26.]# sets x-y limits for images (m) Nominal is 26; use 17 for 3 rings
    CIRCLE_SCALE_BEFORE = 0.001
    CIRCLE_SCALE = 0.001    # scale parameter for circle plots (nominal is 0.001)                                                    
    CIRCLE_SCALE_BAR = 0.2  # set size of scale bar in circle plot (nominal is 0.2 mm)                                                        
    CIRCLE_EXTENT = [-26000., 26000.,-26000., 26000.] # sets x-y limits for circle plot (in mm); use 17000 for 3 rings
    DATA_SIGMA = 0.075      # sets a reference for histrogram drawn on residual distribution (mm)

    # target file
#    TARGET_FILE = 'M11_NEW_TARGET_LIST.csv'
#    TARGET_FILE = 'M1_TARGET_LIST_ELIM.csv'   # this is the file with the targets listed

    #TARGET_FILE = PATH+'New_M1_Target_List_'+MAP+'.csv'
#    TARGET_FILE = 'New_M1_Target_List_2.csv'   # this is the file with the targets listed
#    TARGET_FILE = 'M1_TARGET_LIST_AFTER_7.csv'   # this is the file with the targets listed
    TARGET_FILE = 'M1_TARGET_LIST_from_M72.csv'

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

    if USE_REDUCED_DATA_FILE == 0:
        # get the apex reference target file
        if USE_REDUCED_APEX_FILE == 1:
            print('USING REDUCED APEX FILE: %s'%(F.RED_APEX_FILE_NAME))
            A = apex(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=False)
#            A = apex3(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=False)
        else:
            print('SOLVING WITH RAW APEX FILE: %s'%(F.APEX_FILE_NAME))
            A = apex(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=True)
#            A = apex3(F.APEX_FILE_NAME,F.RED_APEX_FILE_NAME,solve=True)
 
        print('USING RAW POINT CLOUD DATA: %s'%(F.DATA_FILE_NAME))
        M = m1_photo(F.DATA_FILE_NAME,META,parabola_correction=PARABOLA_CORRECTION,raw_read=True)
        M.rotate_cloud(A)
        M.ring_id()
        #M.find_segments()
    else:
        print('USING REDUCED POINT CLOUD DATA: %s'%(F.RED_DATA_FILE_NAME))
        M = m1_photo(F.RED_DATA_FILE_NAME,META,parabola_correction=PARABOLA_CORRECTION,raw_read=False)
        #M.find_segments(print_it=True)

    # load the target list and associate targets with cloud points
    T = m1_targets(TARGET_FILE)
    M.find_target_ids(T)
    M.elim_targets(T)
    M.find_segments()
    T.check_targets()
    
    # set the model to be fit to the cloud from command line args
    M.set_model(ModelType=MODEL_NAME,npar=MODEL_NPARS)

    # MULTI PASS BRANCH OF SCRIPT ****************************************
    if MULTIPASS_DATA_FILE == 1:

        # PASS 0
        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        M.make_residuals_vector(res_rings=SURFACE_FIT_RINGS)
        
        # PASS 1 - crop at 4mm error                                         
        M.crop_points(4.)
        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        M.make_residuals_vector(res_rings=SURFACE_FIT_RINGS)

        # fit the segments
        for i in SURFACE_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                seg.update_residuals_array(M)

        # PASS 2 - crop at 1mm error
        M.crop_points(1.)
        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        M.make_residuals_vector(res_rings=SURFACE_FIT_RINGS)

        # fit the segments
        for i in SURFACE_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                seg.update_residuals_array(M)

        # PASS 3 - crop at 0.5mm error
        M.crop_points(0.5)
        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        M.make_residuals_vector(res_rings=SURFACE_FIT_RINGS)

        # fit the segments
        for i in SURFACE_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                seg.update_residuals_array(M)

        # after final pass we redo fit and compute residuals for full dish
        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        M.make_residuals_vector(res_rings=SURFACE_RESIDUALS_RINGS)

        M.red_write(T,'special.csv')

        M.crop_points(10.) # special crop

        print('BEFORE SEGMENT FIT *****************')
        rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
        print('Z RMS     = %6.3f'%(rms))
        M.make_phase_residuals()
        print('PHASE RMS = %6.3f mm'%(M.PhaseRMS))
        print('Weighted  = %6.3f mm'%(M.WeightedPhaseRMS))
        print('Ring 1 RMS= %6.3f mm'%(M.RING_RMS[0]))
        print('Ring 2 RMS= %6.3f mm'%(M.RING_RMS[1]))
        print('Ring 3 RMS= %6.3f mm'%(M.RING_RMS[2]))
        print('Ring 4 RMS= %6.3f mm'%(M.RING_RMS[3]))
        print('Ring 5 RMS= %6.3f mm'%(M.RING_RMS[4]))
        print('Efficiency= %6.3f'%(M.Efficiency))
        M.print_csv_model_line(MAP_NAME)
        M.print_csv_rms_line(MAP_NAME,rms)

        if PLOT_IT:
            S = m1_show(M)
            S.imshow('%s %s ZRMS=%5.3f'%(MAP_NAME,MODEL_STRING,rms),figure=1,extent=IMAGE_EXTENT,rlimits=IMAGE_SCALE_BEFORE,show_points=False)
            S.circle_plot_m1_residuals('%s %s'%(MAP_NAME,MODEL_STRING),CIRCLE_SCALE_BEFORE,figure=2,extent=CIRCLE_EXTENT,scale_bar=CIRCLE_SCALE_BAR)
    
            # make histogram of all residuals                                         
            pl.figure(10)
            res_bins = np.linspace(-1.,1.,81)
            pl.hist(M.RESIDUALS,res_bins)
            rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
            # gaussian for the data_sigma provided as an argument for comparison    
            y = M.fit_nobs/math.sqrt(2.*math.pi*DATA_SIGMA**2)*np.exp(-res_bins**2/(2.*DATA_SIGMA**2))*(res_bins[1]-res_bins[0])
            pl.plot(res_bins,y,'r',lw=2)
            pl.xlabel('RMS')
            pl.ylabel('N')
            pl.grid()
            pl.title('%s Residuals RMS=%6.3f ModelRMS=%6.3f'%(MAP_NAME,rms,DATA_SIGMA))


        if VIEW_SEGMENT_FITS:
            show_segment = 1
            while show_segment == 1:
                the_ring_to_show = eval(input('Ring'))
                the_segment_to_show = eval(input('Segment'))
                seg = segment_photo(M,the_ring_to_show,the_segment_to_show,META,SEGMENT_DEFORMATION_MODEL)
                ss = segment_show(seg,figure=99)
                ss.make_segment_plot()
                show_segment = eval(input('View Another?'))
                pl.close(99)
        

        # do a fit of all segments
        for i in ACTUATOR_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                seg.update_residuals_array(M)

        # crop the outlier points and recompute residuals
        M.crop_points(0.5)
        M.make_residuals_vector(res_rings=SURFACE_RESIDUALS_RINGS)

        # fit segments and write commands on last iteration
        command_file=open(F.COMMAND_FILE_NAME,'w')
        if PRINT_SCREW_TURNS:
            screw_file=open('%sADJUSTERS_%s.dat'%(PATH,MAP_NAME),'w')
        for i in ACTUATOR_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                if PRINT_SCREW_TURNS and (i>3):
                    #seg.print_command_order()
                    #seg.print_screw_turns()
                    seg.write_screw_turns(screw_file)
                seg.write_command_actuators(command_file,ACTUATOR_COMMAND)
                seg.update_residuals_array(M)
        if PRINT_SCREW_TURNS:
            screw_file.close()
        command_file.close()

        print('AFTER SEGMENT FIT *****************')
        rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
        print('Z RMS     = %6.3f'%(rms))
        M.make_phase_residuals()
        print('PHASE RMS = %6.3f mm'%(M.PhaseRMS))
        print('Weighted  = %6.3f mm'%(M.WeightedPhaseRMS))
        print('Ring 1 RMS= %6.3f mm'%(M.RING_RMS[0]))
        print('Ring 2 RMS= %6.3f mm'%(M.RING_RMS[1]))
        print('Ring 3 RMS= %6.3f mm'%(M.RING_RMS[2]))
        print('Ring 4 RMS= %6.3f mm'%(M.RING_RMS[3]))
        print('Ring 5 RMS= %6.3f mm'%(M.RING_RMS[4]))
        print('Efficiency= %6.3f'%(M.Efficiency))
        M.print_csv_model_line(MAP_NAME)
        M.print_csv_rms_line(MAP_NAME,rms)

        if PLOT_IT:
            S = m1_show(M)
            S.imshow('%s SegFit ZRMS=%5.3f'%(MAP_NAME,rms),figure=3,extent=IMAGE_EXTENT,rlimits=IMAGE_SCALE,show_points=False)
            S.circle_plot_m1_residuals(MAP_NAME+' SegFit',CIRCLE_SCALE,figure=4,extent=CIRCLE_EXTENT,scale_bar=CIRCLE_SCALE_BAR)

            # make histogram of all residuals                                         
            pl.figure(11)
            res_bins = np.linspace(-0.5,0.5,41)
            pl.hist(M.RESIDUALS,res_bins)
            rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
            # gaussian for the data_sigma provided as an argument for comparison    
            y = M.fit_nobs/math.sqrt(2.*math.pi*DATA_SIGMA**2)*np.exp(-res_bins**2/(2.*DATA_SIGMA**2))*(res_bins[1]-res_bins[0])
            pl.plot(res_bins,y,'r',lw=2)
            pl.xlabel('RMS')
            pl.ylabel('N')
            pl.grid()
            pl.title('%s Residuals RMS=%6.3f ModelRMS=%6.3f'%(MAP_NAME,rms,DATA_SIGMA))

    # END OF MULTIPASS BRANCH

    else: # THIS IS THE SINGLE PASS BRANCH ON REDUCED POINT CLOUD*********************

        M.fit_cloud(fit_rings=SURFACE_FIT_RINGS)
        # last pass gives residuals to all rings
        M.make_residuals_vector(res_rings=SURFACE_RESIDUALS_RINGS)

        print('BEFORE SEGMENT FIT *****************')
        rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
        print('Z RMS     = %6.3f'%(rms))
        M.make_phase_residuals()
        print('PHASE RMS = %6.3f mm'%(M.PhaseRMS))
        print('Weighted  = %6.3f mm'%(M.WeightedPhaseRMS))
        print('Ring 1 RMS= %6.3f mm'%(M.RING_RMS[0]))
        print('Ring 2 RMS= %6.3f mm'%(M.RING_RMS[1]))
        print('Ring 3 RMS= %6.3f mm'%(M.RING_RMS[2]))
        print('Ring 4 RMS= %6.3f mm'%(M.RING_RMS[3]))
        print('Ring 5 RMS= %6.3f mm'%(M.RING_RMS[4]))
        print('Efficiency= %6.3f'%(M.Efficiency))
        M.print_csv_model_line(MAP_NAME)
        M.print_csv_rms_line(MAP_NAME,rms)
    
        if PLOT_IT:
            S = m1_show(M)
            S.imshow('%s ZRMS=%5.3f'%(MAP_NAME,rms),figure=1,extent=IMAGE_EXTENT,rlimits=IMAGE_SCALE)
            S.circle_plot_m1_residuals(MAP_NAME,CIRCLE_SCALE,figure=2,extent=CIRCLE_EXTENT,scale_bar=CIRCLE_SCALE_BAR)
    
            pl.ion()
            pl.figure(10)
            res_bins = np.linspace(-1.,1.,81)
            pl.hist(M.RESIDUALS,res_bins)
            rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
            # gaussian for the data_sigma provided as an argument for comparison    
            y = M.fit_nobs/math.sqrt(2.*math.pi*DATA_SIGMA**2)*np.exp(-res_bins**2/(2.*DATA_SIGMA**2))*(res_bins[1]-res_bins[0])
            pl.plot(res_bins,y,'r',lw=2)
            pl.xlabel('RMS')
            pl.ylabel('N')
            pl.grid()
            pl.title('%s Residuals RMS=%6.3f ModelRMS=%6.3f'%(MAP_NAME,rms,DATA_SIGMA))


        # fit segments and write commands
        command_file=open(F.COMMAND_FILE_NAME,'w')
        if PRINT_SCREW_TURNS:
            screw_file=open('%sADJUSTERS_%s.dat'%(PATH,MAP_NAME),'w')
        for i in ACTUATOR_FIT_RINGS:
            for j in range(N_SEGMENTS[i-1]):
                seg = segment_photo(M,i,j+1,META,SEGMENT_DEFORMATION_MODEL)
                if PRINT_SCREW_TURNS and (i>3):
                    #seg.print_screw_turns()
                    seg.write_screw_turns(screw_file)
                seg.write_command_actuators(command_file,ACTUATOR_COMMAND)
                seg.update_residuals_array(M)
        if PRINT_SCREW_TURNS:
            screw_file.close()
        command_file.close()

        print('AFTER SEGMENT FIT *****************')
        rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
        print('Z RMS     = %6.3f'%(rms))
        M.make_phase_residuals()
        print('PHASE RMS = %6.3f mm'%(M.PhaseRMS))
        print('Weighted  = %6.3f mm'%(M.WeightedPhaseRMS))
        print('Ring 1 RMS= %6.3f mm'%(M.RING_RMS[0]))
        print('Ring 2 RMS= %6.3f mm'%(M.RING_RMS[1]))
        print('Ring 3 RMS= %6.3f mm'%(M.RING_RMS[2]))
        print('Ring 4 RMS= %6.3f mm'%(M.RING_RMS[3]))
        print('Ring 5 RMS= %6.3f mm'%(M.RING_RMS[4]))
        print('Efficiency= %6.3f'%(M.Efficiency))
        M.print_csv_model_line(MAP_NAME)
        M.print_csv_rms_line(MAP_NAME,rms)

        if PLOT_IT:
            S = m1_show(M)
            S.imshow('%s SegFit ZRMS=%5.3f'%(MAP_NAME,rms),figure=3,extent=IMAGE_EXTENT,rlimits=IMAGE_SCALE)
            S.circle_plot_m1_residuals(MAP_NAME+' SegFit',CIRCLE_SCALE,figure=4,extent=CIRCLE_EXTENT,scale_bar=CIRCLE_SCALE_BAR)

            # make histogram of all residuals                                         
            pl.figure(11)
            res_bins = np.linspace(-0.5,0.5,41)
            pl.hist(M.RESIDUALS,res_bins)
            rms = math.sqrt(np.dot(M.RESIDUALS.transpose(),M.RESIDUALS)/M.fit_nobs)
            # gaussian for the data_sigma provided as an argument for comparison    
            y = M.fit_nobs/math.sqrt(2.*math.pi*DATA_SIGMA**2)*np.exp(-res_bins**2/(2.*DATA_SIGMA**2))*(res_bins[1]-res_bins[0])
            pl.plot(res_bins,y,'r',lw=2)
            pl.xlabel('RMS')
            pl.ylabel('N')
            pl.grid()
            pl.title('%s Residuals RMS=%6.3f ModelRMS=%6.3f'%(MAP_NAME,rms,DATA_SIGMA))

    # END SINGLE PASS BRANCH **********************************************

    # write reduced file
    M.red_write(T,F.RED_DATA_FILE_NAME)



main(sys.argv[1:])
