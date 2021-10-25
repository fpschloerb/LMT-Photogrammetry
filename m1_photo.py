For """Module with m1_photo class for solving photogrammetry cloud for m1

Classes: m1_photo
Uses:    math,numpy
Author:  FPS
Date:    November 6, 2016
Changes: November 2017 - added rings 4 and 5 to analysis
                         added blanking of segments with no targets 
         December 2017 - eliminates points not found within 30mm of target

Parameters:
SETUP AND BASIC DATA ******************************************************
NOTE: values with * may be read in from Reduced file
metadata_date             __init__  string with map date
actuator_postions[720,3]  __init__  beginning(0), end(1), average position(2)
temperatures[64,3]        __init__  beginning,end,average temperatures
area[5]                   __init__  computed area of each ring
nobs            *         raw_read  number of data points
data[nobs,3]              raw_read  x,y,z of data points
id[nobs]        *         raw_read  id number (order in file)
old_id[nobs]    *         raw_read  id number in original raw data file 
valid[nobs]     *         raw_read  flag=1 if point is "valid"
valid_indices[] *         raw_read  array of indices where valid=1
valid_nobs      *         raw_read  number of data points with valid set to 1
select[nobs]    *         raw_read  initialized to values of valid
residuals[nobs] *         red_read  residuals vector with nonfit points set 
                                      to value of 9999
XVec[nobs,3]    *         rotate_cloud point cloud data points rotated to 
                                       dish coordinate system
ring[nobs]      *         ring_id  ring identification for data points
seg[nobs]       *         ring_id  segment identification for data points
found_segment_1[12]       find_segments  flag to show whether there are data 
                                         points for the segment in ring 1
found_segment_2[24]       find_segments  flags for ring 2
found_segment_3[48]       find_segments  flags for ring 3 
found_segment_4[48]       find_segments  flags for ring 4 
found_segment_5[48]       find_segments  flags for ring 5
target_id[nobs]           find_target_ids contains target id number for each 
                                          data point
d_target[nobs]            find_target_ids distance of data point from 
                                          identified target
FITTING THE MODEL **********************************************************
fit_indices[]             set_fit_indices array of indices to valid points 
                                          for fit
fit_nobs                  set_fit_indices number of points in fit
i1,i2,i3,i4,i5            set_fit_indices number of valid fit points in 
                                          rings 1-5 
w1,w2,w3,w4,w5            set_fit_indices weight factor for each ring for 
                                          computation of weighted phase rms
model_type                set_model  string with the type of model 
                                     (SIMPLE,ZERNIKE,TILT)
npar                      set_model  number of parameters in model
P[npar]                   fit_cloud  parameter values of fit
fit_rms                   fit_cloud  rms of the fit
E[npar]                   fit_cloud  error of fit parameters from cov. matrix
Presult[npar]             fit_cloud  P converted to "natural" units
Eresult[npar]             fit_cloud  E converted to "natural" units
MODEL RESIDUALS ************************************************************
X[fit_nobs]               make_residuals_vector array with X values using 
                                                valid points from fit
Y[fit_nobs]               make_residuals_vector array with Y values using 
                                                valid points from fit
Phase[fit_nobs]           make_residuals_vector array of factors to go from z 
                                                errors to half-path errors
D[fit_nobs]               make_residuals_vector array of data for fit
ID[fit_nobs]              make_residuals_vector data point id numbers for 
data points in fit
RING[fit_nobs]            make_residuals_vector ring id for data points in fit
SEG[fit_nobs]             make_residuals_vector seg id for data points in fit
RESIDUALS[fit_nobs]       make_residuals_vector residuals to the fit, using 
                                                vector P
PHASE_RESIDUALS[fit_nobs] make_phase_residuals  half-path residuals vector
Phase_RMS                 make_phase_residuals  rms of half-path residuals 
                                                vector
Efficiency                make_phase_residuals  calculated efficiency of 
                                                antenna - illumination weighted
DENOM                     make_phase_residuals  dummy variable
RING_RMS[5]               make_phase_residuals  computed rms of phase 
                                                residuals in each ring
RING_COUNT[5]             make_phase_residuals  number of points in each ring 
                                                for above computation
WeightedPhaseRMS          make_phase_residuals  half-path rms computed from 
                                                Efficiency
"""
import numpy as np
import math
from m1_targets import m1_targets
from photo_meta import photo_meta

class m1_photo():
    """m1_photo manages point cloud from photometry"""
    def __init__(self,filename,metadata,parabola_correction=False,raw_read=True):

        # get metadata information
        self.metadata_date = metadata.get_date()
        self.actuator_positions = metadata.get_actuators()
        self.temperatures = metadata.get_temperatures()
        self.elevation = metadata.get_elevation()

        if(raw_read):
            self.raw_read(filename)
        else:
            self.red_read(filename)

        self.area = np.zeros(5)
        self.area[0] = np.pi*(6650.**2-1625.**2)
        self.area[1] = np.pi*(11550.**2-6650.**2)
        self.area[2] = np.pi*(16250.**2-11550.**2)
        self.area[3] = np.pi*(20750.**2-16250.**2)
        self.area[4] = np.pi*(25000.**2-20750.**2)

        # set the model values for correction
        self.parabola_correction = parabola_correction
        self.model_focal_length = 17485.454 + 13.383*math.sin(self.elevation)
        self.model_y_tilt = (-12.60 + 69.78*math.cos(self.elevation))/206265.
        print('PARABOLA CORRECTION: ',self.parabola_correction)
        print('MODEL: %f %f'%(self.model_focal_length,self.model_y_tilt*206265.))

    def red_read(self,filename):
        input_data = np.loadtxt(filename)
        if(len(input_data[0,:]) != 9):
            print('WARNING: %s is not a reduced file'%(filename))
        self.nobs = len(input_data[:,0])
        self.target_id = np.array(input_data[:,0],dtype='int')
        self.ring = np.array(input_data[:,1],dtype='int')
        self.seg = np.array(input_data[:,2],dtype='int')
        self.XVec = np.zeros((self.nobs,3))
        self.XVec[:,0] = input_data[:,3]
        self.XVec[:,1] = input_data[:,4]
        self.XVec[:,2] = input_data[:,5]
        self.id = np.array(input_data[:,0],dtype='int')
        self.old_id = np.array(input_data[:,6],dtype='int')
        self.valid = np.array(input_data[:,7],dtype='int')
        self.residuals = np.array(input_data[:,8])
        self.valid_indices = np.where(self.valid == 1)[0]
        self.valid_nobs = len(self.valid_indices)
        self.select = self.valid

    def red_write(self,Targets,filename):
        F = np.zeros(11)
        f = open(filename,'w')
        for i in range(Targets.n):
            if Targets.found_target[i] == 1:   # only write targets with data
                j = Targets.point_id[i]
                index = np.where(self.ID==j)[0]
                if len(index) == 1:
                    resid = self.RESIDUALS[index]
                else:
                    resid = 9999.
                f.write('%d %d %d %f %f %f %d %d %f \n'
                        %
                        (i,
                         self.ring[j],
                         self.seg[j],
                         self.XVec[j,0],
                         self.XVec[j,1],
                         self.XVec[j,2],
                         self.id[j],
                         self.valid[j],
                         resid)
                        )
            else:
                #print('red_write: target %d found %d points'
                #      %(i,Targets.found_target[i]))
                f.write('%d %d %d 0.0 0.0 0.0 %d 0 9999. \n'
                        %
                        (i,
                         self.ring[j],
                         self.seg[j],
                         self.id[j])
                        )
        f.close()

    def write_target_file(self,filename):
        # we assume that file has been rotated etc. first
        target_list = []
        for i in range(self.nobs):
            target_list.append( 
                (self.ring[i],
                 self.seg[i],
                 self.XVec[i,0],
                 self.XVec[i,1]) 
              )
        target_type = [
            ('ring','int'), 
            ('seg','int'), 
            ('x','float'), 
            ('y','float')
            ]
        b = np.array(target_list, dtype = target_type)
        print(b)
        a = np.sort(b, order=['ring','seg','x'])
        f = open(filename,'w')
        for i in range(self.nobs):
            f.write('%d,%d,%f,%f,1\n'%(a[i][0],a[i][1],a[i][2],a[i][3]))
            print('%d %d %d %f %f'%(i,a[i][0],a[i][1],a[i][2],a[i][3]))
        f.close()

    def elim_points(self,id_list):
        for i in id_list:
            self.valid[i] = 0
            self.select[i] = 0
            print('Eliminate %d'%(i))
        self.valid_indices = np.where(self.valid == 1)[0]
        self.valid_nobs = len(self.valid_indices)

    def elim_targets(self,Targets):
        for i in range(Targets.n):
            if Targets.flag[i] == 0:
                self.valid[Targets.point_id[i]] = 0
                self.select[Targets.point_id[i]] = 0
                print('Eliminate Target %d Point %d'%(i,Targets.point_id[i]))
        self.valid_indices = np.where(self.valid == 1)[0]
        self.valid_nobs = len(self.valid_indices)

    def raw_read(self,filename):
        """ reads original data file """
        input_data = np.loadtxt(filename)
        if(len(input_data[0,:]) > 3):
           print('WARNING: %s is not an original data file'%(filename))
        self.nobs = len(input_data[:,0])
        self.data = np.zeros((self.nobs,3))
        self.data[:,0] = input_data[:,1] # x is second element 
        self.data[:,1] = input_data[:,2] # y is third element
        self.data[:,2] = input_data[:,0] # z is first element
        self.id = np.array(np.zeros(self.nobs),dtype='int')
        self.old_id = np.array(np.zeros(self.nobs),dtype='int')
        for i in range(self.nobs):
            self.id[i] = i
            self.old_id[i] = i
        self.valid = np.array(np.ones(self.nobs),dtype='int')
        self.valid_indices = np.where(self.valid == 1)[0]
        self.valid_nobs = len(self.valid_indices)
        self.select = self.valid

    def rotate_cloud(self,apex_reference):
        self.XVec = np.zeros((self.nobs,3))
        for i in range(self.nobs):
            self.XVec[i,:] = apex_reference.rotate(
                apex_reference.model,
                self.data[i,:]
               )
            
    def ring_id(self,print_it=False):
        self.ring = np.zeros(self.nobs,dtype='int')
        self.seg = np.zeros(self.nobs,dtype='int')
        for i in range(self.nobs):
            radius = math.sqrt(self.XVec[i,0]**2+self.XVec[i,1]**2)
            theta = math.atan2(self.XVec[i,0],self.XVec[i,1])/math.pi*180.
            if(radius < 6650.):
                self.ring[i] = 1
                if(theta<-15.):
                    theta = theta + 360.
                self.seg[i] = np.floor((theta+15)/30.) + 1
            elif (radius < 11550.):
                self.ring[i] = 2
                if(theta<0.):
                    theta = theta + 360.
                self.seg[i] = np.floor(theta/15.) + 1
            elif (radius < 16250.):
                self.ring[i] = 3
                if(theta<0.):
                    theta = theta + 360.
                self.seg[i] = np.floor(theta/7.5) + 1
            elif (radius < 20750.):
                self.ring[i] = 4
                if(theta<0.):
                    theta = theta + 360.
                self.seg[i] = np.floor(theta/7.5) + 1
            else:
                self.ring[i] = 5
                if(theta<0.):
                    theta = theta + 360.
                self.seg[i] = np.floor(theta/7.5) + 1
            if print_it:
                print('%d %f %d %f %d %f %f %f'
                      %
                      (self.id[i],
                       radius,
                       self.ring[i],
                       theta,
                       self.seg[i],
                       self.XVec[i,0],
                       self.XVec[i,1],
                       self.XVec[i,2])
                     )

    def find_segments(self,print_it=False):
        self.found_segment_1 = np.zeros(12,dtype='int')
        self.found_segment_2 = np.zeros(24,dtype='int')
        self.found_segment_3 = np.zeros(48,dtype='int')
        self.found_segment_4 = np.zeros(48,dtype='int')
        self.found_segment_5 = np.zeros(48,dtype='int')

        for i in range(self.nobs):
            if self.valid[i] == 1:
                radius = math.sqrt(self.XVec[i,0]**2+self.XVec[i,1]**2)
                theta = math.atan2(self.XVec[i,0],self.XVec[i,1])/math.pi*180.
                if(radius < 6650.):
                    if(theta<-15.):
                        theta = theta + 360.
                    sss = int(np.floor((theta+15)/30.))
                    self.found_segment_1[sss] += 1
                elif (radius < 11550.):
                    if(theta<0.):
                        theta = theta + 360.
                    sss = int(np.floor(theta/15.))
                    self.found_segment_2[sss] += 1
                elif (radius < 16250.):
                    if(theta<0.):
                        theta = theta + 360.
                    sss = int(np.floor(theta/7.5))
                    self.found_segment_3[sss] += 1
                elif (radius < 20750.):
                    if(theta<0.):
                        theta = theta + 360.
                    sss = int(np.floor(theta/7.5))
                    self.found_segment_4[sss] += 1
                else:
                    if(theta<0.):
                        theta = theta + 360.
                    sss = int(np.floor(theta/7.5))
                    self.found_segment_5[sss] += 1
        
        if print_it == True:
            for i in range(12):
                if self.found_segment_1[i] == 0:
                    print('No data for Ring 1 Segment %d'%(i+1))
            for i in range(24):
                if self.found_segment_2[i] == 0:
                    print('No data for Ring 2 Segment %d'%(i+1))
            for i in range(48):
                if self.found_segment_3[i] == 0:
                    print('No data for Ring 3 Segment %d'%(i+1))
            for i in range(48):
                if self.found_segment_4[i] == 0:
                    print('No data for Ring 4 Segment %d'%(i+1))
            for i in range(48):
                if self.found_segment_5[i] == 0:
                    print('No data for Ring 5 Segment %d'%(i+1))

    def find_target_ids(self,Targets):
        self.target_id = -1*np.ones(self.nobs,dtype='int')
        self.d_target = np.zeros(self.nobs)
        for i in range(self.nobs):
            which_target = -1
            dmin = 99999.
            for j in range(Targets.n):
                if Targets.ring[j] == self.ring[i]:
                    if Targets.seg[j] == self.seg[i]:
                        d = np.sqrt(
                              (self.XVec[i,0]-Targets.x[j])**2 
                            + (self.XVec[i,1]-Targets.y[j])**2
                            )
                        if d < dmin:
                            dmin = d
                            which_target = j
            if dmin > 30:
                self.target_id[i] = -1
                self.d_target[i] = dmin
                self.valid[i] = 0
                Targets.set_point_id(-1,i)
                if math.sqrt(self.XVec[i,0]**2+self.XVec[i,1]**2+self.XVec[i,2]**2) > 10.:
                    print('Target %d: closest Point %d is outside tolerance: %f'
                      %
                      (which_target,i,dmin)
                     )
#            if which_target == -1:
#                self.target_id[i] = which_target
#                self.d_target[i] = dmin
#                Targets.set_point_id(which_target,i)
            else:
                if Targets.found_target[which_target] == 0:
                    self.target_id[i] = which_target
                    self.d_target[i] = dmin
                    Targets.set_point_id(which_target,i)
                else:
                    if self.d_target[Targets.point_id[which_target]] > dmin:
                        print('Target %d: Point %d overrides %d'
                              %
                              (which_target,
                               i,
                               Targets.point_id[which_target])
                             )
                        self.target_id[i] = which_target
                        self.d_target[i] = dmin
                        Targets.set_point_id(which_target,i)
                    else:
                        print('Target %d: Point %d duplicates target with error of %f'
                              %
                              (which_target,
                               i,
                               dmin)
                             )
                        self.valid[i] = 0
#            print('point %d  target %d   dmin=%f valid=%d'
#                  %
#                  (i,
#                   which_target,
#                   dmin,
#                   self.valid[i])
#                 )
        self.valid_indices = np.where(self.valid == 1)[0]
        self.valid_nobs = len(self.valid_indices)

    def select_segment_indices(self,ring,seg):
        list = []
        for i in range(self.nobs):
            if((self.ring[i] == ring) and (self.seg[i] == seg)):
                if(self.valid[i] == 1):
                    list.append(self.id[i])
        return(np.array(list,dtype='int'))
        
    def select_residual_segment_indices(self,ring,seg):
        list = []
        for i in range(self.fit_nobs):
            if((self.RING[i] == ring) and (self.SEG[i] == seg)):
                if(self.valid[self.ID[i]] == 1):
                    list.append(i)
        return(np.array(list,dtype='int'))
        
    def set_fit_indices(self,fit_rings):
        count = 0
        index = []
        for i in range(self.nobs):
            if self.valid[i] == 1:
                for j in range(len(fit_rings)):
                    if self.ring[i] == fit_rings[j]:
                        index.append(i)
                        count = count + 1
        
        self.fit_indices = np.array(index)
        self.fit_nobs = len(self.fit_indices)

        self.i5 = 0.
        self.w5 = 0.
        self.i4 = 0.
        self.w4 = 0.
        self.i3 = 0.
        self.w3 = 0.
        self.i2 = 0.
        self.w2 = 0.
        self.i1 = 0.
        self.w1 = 0.

        for theRing in fit_rings:
            if theRing == 5:
                iii = np.where(self.ring[self.fit_indices]==5)
                self.i5 = len(iii[0])
                self.w5 = self.area[4]/self.i5
                #print('ring 5 defined %d %f'%(self.i5,self.w5))
            if theRing == 4:
                iii = np.where(self.ring[self.fit_indices]==4)
                self.i4 = len(iii[0])
                self.w4 = self.area[3]/self.i4
                #print('ring 4 defined %d %f'%(self.i4,self.w4))
            if theRing == 3:
                iii = np.where(self.ring[self.fit_indices]==3)
                self.i3 = len(iii[0])
                self.w3 = self.area[2]/self.i3
                #print('ring 3 defined %d %f'%(self.i3,self.w3))
            if theRing == 2:
                iii = np.where(self.ring[self.fit_indices]==2)
                self.i2 = len(iii[0])
                self.w2 = self.area[1]/self.i2
                #print('ring 2 defined %d %f'%(self.i2,self.w2))
            if theRing == 1:
                iii = np.where(self.ring[self.fit_indices]==1)
                self.i1 = len(iii[0])
                self.w1 = self.area[0]/self.i1
                #print('ring 1 defined %d %f'%(self.i1,self.w1))
    
    def set_model(self,ModelType='SIMPLE',npar=4):
        self.model_type = ModelType
        self.npar = npar

    def set_basis_data(self,i):
        if self.parabola_correction == True:
            result = self.XVec[i,2] - (self.XVec[i,0]**2+self.XVec[i,1]**2)/4./self.model_focal_length - self.model_y_tilt*self.XVec[i,1]
            print('parabola correction %f %f'%(result,self.XVec[i,2]))
        else:
            if self.model_type == 'SIMPLE':
                result = self.XVec[i,2]
            elif self.model_type == 'ZERNIKE':
                result = self.XVec[i,2]
            elif self.model_type == 'TILT':
                result = self.XVec[i,2] - (self.XVec[i,0]**2+self.XVec[i,1]**2)/4./self.model_focal_length
            else:
                result = 0.
        return result
            
    def set_basis_functions(self,i,F):
        if self.model_type == 'SIMPLE':
            F[0] = 1.
            F[1] = self.XVec[i,0]**2 + self.XVec[i,1]**2
            F[2] = self.XVec[i,0]
            F[3] = self.XVec[i,1]
        elif self.model_type == 'ZERNIKE':
            rho = math.sqrt(self.XVec[i,0]**2 + self.XVec[i,1]**2)/25000.
            theta = math.atan2(self.XVec[i,1],self.XVec[i,0])
            F[0] = 1.                                                # Piston               [0,  0]
            F[1] = 2.*rho*math.cos(theta)                            # Tip (X)              [1, +1]
            F[2] = 2.*rho*math.sin(theta)                            # Tilt (Y)             [1, -1]
            F[3] = math.sqrt(3.)*(2.*rho**2-1.)                      # Defocus              [2,  0]
            F[4] = math.sqrt(6.)*rho**2*math.cos(2.*theta)           # Vertical Astigmatism [2, +2]
            F[5] = math.sqrt(6.)*rho**2*math.sin(2.*theta)           # Oblique Astigmatism  [2, -2]
            F[6] = math.sqrt(8.)*(3.*rho**3-2.*rho)*math.cos(theta)  # Horizontal Coma      [3, +1]
            F[7] = math.sqrt(8.)*(3.*rho**3-2.*rho)*math.sin(theta)  # Vertical Coma        [3, -1]
            F[8] = math.sqrt(5.)*(6.*rho**4-6.*rho**2+1.)            # Primary Spherical    [4,  0]
            F[9] = math.sqrt(8.)*rho**3*math.cos(3.*theta)           # Oblique Trefoil     [3, +3]
            F[10] = math.sqrt(8.)*rho**3*math.sin(3.*theta)          # Vertical Trefoil    [3, -3]
            F[11] = math.sqrt(10.)*(4.*rho**4-3.*rho**2)*math.cos(2.*theta) # Vertical secondary astigmatism [4, +2]
            F[12] = math.sqrt(10.)*(4.*rho**4-3.*rho**2)*math.sin(2.*theta) # Oblique secondary  astigmatism [4, -2]
            F[13] = math.sqrt(10.)*rho**4*math.cos(4.*theta)         # Vertical Quadrafoil  [4, +4]
            F[14] = math.sqrt(10.)*rho**4*math.sin(4.*theta)         # Oblique Quadrafoil   [4, -4]
            F[15] = math.sqrt(12.)*(10.*rho**5-12.*rho**3+3.*rho)*math.cos(theta) #         [5, +1]
            F[16] = math.sqrt(12.)*(10.*rho**5-12.*rho**3+3.*rho)*math.sin(theta) #         [5, -1]
            F[17] = math.sqrt(12.)*(5.*rho**5-4.*rho**3)*math.cos(3.*theta)       #         [5, +3]
            F[18] = math.sqrt(12.)*(5.*rho**5-4.*rho**3)*math.sin(3.*theta)       #         [5, -3]
            F[19] = math.sqrt(12.)*rho**5*math.cos(5.*theta)                      #         [5, +5]
            F[20] = math.sqrt(12.)*rho**5*math.sin(5.*theta)                      #         [5, -5]
            F[21] = math.sqrt(7.)*(20.*rho**6-30.*rho**4+12.*rho**2-1.0)          #         [6,  0]
            F[22] = math.sqrt(14.)*(15.*rho**6-20.*rho**4+6.*rho**2)*math.cos(2.*theta) #   [6, +2]
            F[23] = math.sqrt(14.)*(15.*rho**6-20.*rho**4+6.*rho**2)*math.sin(2.*theta) #   [6, -2]
            F[24] = math.sqrt(14.)*(6.*rho**6-5.*rho**4)*math.cos(4.*theta)       #         [6, +4]
            F[25] = math.sqrt(14.)*(6.*rho**6-5.*rho**4)*math.sin(4.*theta)       #         [6, -4]
            F[26] = math.sqrt(14.)*rho**6*math.cos(6.*theta)                      #         [6, +6]
            F[27] = math.sqrt(14.)*rho**6*math.sin(6.*theta)                      #         [6, +6]

        elif self.model_type == 'TILT':
            F[0] = 1.
            F[1] = self.XVec[i,0]
            F[2] = self.XVec[i,1]
        else:
            print('MODEL TYPE ERROR')
        
    def fit_cloud(self,fit_rings=[1,2,3]):
        F = np.zeros(28)  # largest number of parameters
        
        # select observations according to ring number 
        # (specify outermost ring above)
        self.set_fit_indices(fit_rings)

        H = np.zeros((self.fit_nobs,self.npar))
        D = np.zeros(self.fit_nobs)

        ipt = 0
        for i in self.fit_indices:
            D[ipt] = self.set_basis_data(i)
            self.set_basis_functions(i,F)
            for j in range(self.npar):
                H[ipt,j] = F[j]
            ipt = ipt + 1

        HTH = np.dot(np.transpose(H),H)
        HTD = np.dot(np.transpose(H),D)

        HTHINV = np.linalg.inv(HTH)

        self.P = np.dot(HTHINV,HTD)

        RESIDUALS = D - np.dot(H,self.P)
        CHISQ = np.dot(np.transpose(RESIDUALS),RESIDUALS)
        DATA_VARIANCE = CHISQ/self.fit_nobs 
        COV_MATRIX = DATA_VARIANCE*HTHINV
        self.fit_rms = math.sqrt(DATA_VARIANCE)

        # save errors and prepare array for printed results
        self.E = np.zeros(self.npar)
        self.Eresult = np.zeros(self.npar)
        self.Presult = np.zeros(self.npar)
        for i in range(self.npar):
            self.E[i] = math.sqrt(COV_MATRIX[i,i])

        # print the results in the "natural units" for each type of fit
        print('Number of points in fit: %d of %d'
              %(self.fit_nobs,
                self.valid_nobs)
             )
        if self.model_type == 'SIMPLE':
            print('Simple fit with %d parameters'%(self.npar))
            if self.npar >= 2:
                self.Presult[0] = self.P[0]
                self.Eresult[0] = self.E[0] # added Dec 6, 2017
                print('Delta Z = %f (%f) millimeters'
                      %(self.Presult[0],self.Eresult[0]))
                self.Presult[1] = 1./4./self.P[1]
                self.Eresult[1] = 1./4./self.P[1]*(self.E[1]/self.P[1])
                print('Focal Length = %f (%f) millimeters'
                      %(self.Presult[1],self.Eresult[1]))
            if self.npar == 4:
                self.Presult[2] = self.P[2]*206265.
                self.Presult[3] = self.P[3]*206265.
                self.Eresult[2] = self.E[2]*206265.
                self.Eresult[3] = self.E[3]*206265.
                print('Delta TX = %f (%f) arcsec'
                      %(self.Presult[2],self.Eresult[2]))
                print('Delta TY = %f (%f) arcsec'
                      %(self.Presult[3],self.Eresult[3]))
        elif self.model_type == 'ZERNIKE':
            for i in range(self.npar):
                self.Presult[i] = self.P[i]
                self.Eresult[i] = self.E[i]
                print('Z[%d] = %e (%e)'
                      %(i,self.Presult[i],self.Eresult[i]))
        elif self.model_type == 'TILT':
            self.Presult[0] = self.P[0]
            self.Presult[1] = self.P[1]*206265.
            self.Presult[2] = self.P[2]*206265.
            self.Eresult[0] = self.E[0]
            self.Eresult[1] = self.E[1]*206265.
            self.Eresult[2] = self.E[2]*206265.
            print('Delta Z  = %f (%f) millimeters'
                  %(self.Presult[0],self.Eresult[0]))
            print('Delta TX = %f (%f) arcsec'
                  %(self.Presult[1],self.Eresult[1]))
            print('Delta TY = %f (%f) arcsec'
                  %(self.Presult[2],self.Eresult[2]))
        print('RMS = %f millimeters'%(self.fit_rms))


    def crop_points(self,limit):
        for i in range(self.fit_nobs):
            if(abs(self.RESIDUALS[i]) > limit):
                self.valid[self.ID[i]] = 0
                print('Ring=%d Segment=%d Target=%d Point=%d Value=%6.3f eliminated'
                      %(self.ring[self.ID[i]],
                        self.seg[self.ID[i]],
                        self.target_id[self.ID[i]],
                        self.id[self.ID[i]],
                        self.RESIDUALS[i])
                     )

    def make_residuals_vector(self,res_rings=[1,2,3]):
        """ create residuals vectors after model fit """
        F = np.zeros(28)  # largest number of parameters

        # select observations according to ring number 
        # (specify outermost ring above)
        self.set_fit_indices(res_rings)

        H = np.zeros((self.fit_nobs,self.npar))

        self.X = np.zeros(self.fit_nobs)
        self.Y = np.zeros(self.fit_nobs)
        self.Phase = np.zeros(self.fit_nobs)
        self.D = np.zeros(self.fit_nobs)
        self.ID = np.zeros(self.fit_nobs,dtype='int')
        self.RING = np.zeros(self.fit_nobs,dtype='int')
        self.SEG = np.zeros(self.fit_nobs,dtype='int')

        ipt = 0
        for i in self.fit_indices:
            self.X[ipt] = self.XVec[i,0]
            self.Y[ipt] = self.XVec[i,1]
            self.Phase[ipt] = 1./(1.+(self.X[ipt]**2+self.Y[ipt]**2)/4./17500.**2)
            self.D[ipt] = self.set_basis_data(i)
            self.ID[ipt] = self.id[i]
            self.RING[ipt] = self.ring[i]
            self.SEG[ipt] = self.seg[i]

            self.set_basis_functions(i,F)
            for j in range(self.npar):
                H[ipt,j] = F[j]
            ipt = ipt + 1

        self.RESIDUALS = self.D - np.dot(H,self.P)
        
    def make_phase_residuals(self):
        self.PHASE_RESIDUALS = self.RESIDUALS*self.Phase
        self.PhaseRMS = np.sqrt(
            np.dot(
                self.PHASE_RESIDUALS.transpose(),
                self.PHASE_RESIDUALS
                )/self.fit_nobs
            )
        
        self.Efficiency = 0.
        self.DENOM = 0.
        self.RING_RMS = np.zeros(5)
        self.RING_COUNT = np.zeros(5)
        WAVELENGTH = 3.333 # 90 GHz
        # for three rings with RSR scale was 10700. 
        for i in range(len(self.fit_indices)):
            #factor = np.exp(np.sqrt(self.X[i]**2+self.Y[i]**2)/10700.) # 20k is a guess
            factor = 1. # uniform illumination
            phi = np.cos(4.*np.pi/WAVELENGTH * self.PHASE_RESIDUALS[i])
            if self.RING[i] == 3:
                self.Efficiency = self.Efficiency + phi*self.w3*factor
                self.DENOM = self.DENOM + self.w3*factor
                self.RING_RMS[2] = self.RING_RMS[2] + self.PHASE_RESIDUALS[i]**2
                self.RING_COUNT[2] = self.RING_COUNT[2] + 1
            elif self.RING[i] == 2:
                self.Efficiency = self.Efficiency + phi*self.w2*factor
                self.DENOM = self.DENOM + self.w2*factor
                self.RING_RMS[1] = self.RING_RMS[1] + self.PHASE_RESIDUALS[i]**2
                self.RING_COUNT[1] = self.RING_COUNT[1] + 1
            elif self.RING[i] == 1:
                self.Efficiency = self.Efficiency + phi*self.w1*factor
                self.DENOM = self.DENOM + self.w1*factor
                self.RING_RMS[0] = self.RING_RMS[0] + self.PHASE_RESIDUALS[i]**2
                self.RING_COUNT[0] = self.RING_COUNT[0] + 1
            elif self.RING[i] == 4:
                self.Efficiency = self.Efficiency + phi*self.w4*factor
                self.DENOM = self.DENOM + self.w4*factor
                self.RING_RMS[3] = self.RING_RMS[3] + self.PHASE_RESIDUALS[i]**2
                self.RING_COUNT[3] = self.RING_COUNT[3] + 1        
            elif self.RING[i] == 5:
                self.Efficiency = self.Efficiency + phi*self.w5*factor
                self.DENOM = self.DENOM + self.w5*factor
                self.RING_RMS[4] = self.RING_RMS[4] + self.PHASE_RESIDUALS[i]**2
                self.RING_COUNT[4] = self.RING_COUNT[4] + 1        
        self.Efficiency = (self.Efficiency/self.DENOM)**2
        self.WeightedPhaseRMS = WAVELENGTH/(4.*np.pi) * np.sqrt(-np.log(self.Efficiency))
        self.RING_RMS[0] = np.sqrt(self.RING_RMS[0]/self.RING_COUNT[0])
        self.RING_RMS[1] = np.sqrt(self.RING_RMS[1]/self.RING_COUNT[1])
        self.RING_RMS[2] = np.sqrt(self.RING_RMS[2]/self.RING_COUNT[2])
        self.RING_RMS[3] = np.sqrt(self.RING_RMS[3]/self.RING_COUNT[3])
        self.RING_RMS[4] = np.sqrt(self.RING_RMS[4]/self.RING_COUNT[4])

    def print_csv_rms_line(self,id,final_rms):
        print('%s,%s,%f,%f,%f,%f,%f,%f,%f'
              %
              (id,
               self.model_type,
               final_rms,
               self.WeightedPhaseRMS,
               self.RING_RMS[0],
               self.RING_RMS[1],
               self.RING_RMS[2],
               self.RING_RMS[3],
               self.RING_RMS[4])
             )
        print(' ')
    
    def print_csv_model_line(self,id):
        # obviously this could be done better.....
        if self.npar == 2:
            print('%s,%s,%f,%f'%(id,
                                 self.model_type,
                                 self.Presult[0],
                                 self.Presult[1]))
            print(' ')
        elif self.npar == 4:
            print('%s,%s,%f,%f,%f,%f'%(id,
                                       self.model_type,
                                       self.Presult[0],
                                       self.Presult[1],
                                       self.Presult[2],
                                       self.Presult[3]))
            print(' ')
        elif self.npar == 3:
            print('%s,%s,%f,%f,%f'%(id,
                                       self.model_type,
                                       self.Presult[0],
                                       self.Presult[1],
                                       self.Presult[2]))
            print(' ')
        elif self.npar == 6:
            print('%s,%s,%f,%f,%f,%f,%f,%f'%(id,
                                       self.model_type,
                                       self.Presult[0],
                                       self.Presult[1],
                                       self.Presult[2],
                                       self.Presult[3],
                                       self.Presult[4],
                                       self.Presult[5]))
            print(' ')
        elif self.npar == 8:
            print('%s,%s,%f,%f,%f,%f,%f,%f,%f,%f'%(id,
                                       self.model_type,
                                       self.Presult[0],
                                       self.Presult[1],
                                       self.Presult[2],
                                       self.Presult[3],
                                       self.Presult[4],
                                       self.Presult[5],
                                       self.Presult[6],
                                       self.Presult[7]))
            print(' ')
        else:
            print('MODEL OPTION NOT VALID')

