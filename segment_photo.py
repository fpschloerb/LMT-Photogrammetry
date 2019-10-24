"""Module with segment_photo class for solving shape for individual segment

Classes: segment_photo
Uses:    math,numpy
Author:  FPS
Date:    November 6, 2016
Changes: December 2017 - added manual adjuster screw turns
         December 15, 2017 - fixed error in HOLOG model of deformation

Parameters:
ring                    __init__ ring of the segment under consideration
seg                     __init__ segment id of segment in the ring
deformation_model       __init__ 0=FEM 1=HOL id's the model used for fit.
deformation_model_label __init__ label for which model is being fit: HOL or FEM
nseg                    __init__ number of points for segment fit
id                      __init__ indices of the points on the segment
xr[nseg]                __init__ x position in segment coordinates
yr[nseg]                __init__ y position in segment coordinates
zr[nseg]                __init__ z position in segment coordinates
prefit_rms              __init__ rms of segment points before fit
actuator_positions[4]   __init__ array of actuator position solutions 
segment_residuals[nseg] __init__ array of residuals at segment points  
rms                     __init__ rms of segment points after fit
"""
import numpy as np
import math

class segment_photo():
    def __init__(self,M1,ring,seg,metadata,model=0):
        # set ring and segment id
        self.ring = ring
        self.seg = seg
        if model == 0:
            self.deformation_model = 0
            self.deformation_model_label = 'FEM'
        else:
            self.deformation_model = 1
            self.deformation_model_label = 'HOL'
        
        # initialize the actuator values from metadata
        self.initialize_actuator_values(metadata)

        # get the indices of the segment data from the residuals list
        segment_indices = M1.select_residual_segment_indices(ring,seg)
        self.nseg = len(segment_indices)
        # get indices of original data file
        self.id = M1.select_segment_indices(ring,seg)

        # get the model
        ac,rotation,center = self.look_up_model()

        # correct data for segment coordinates
        ct = math.cos(rotation/180.*math.pi)
        st = math.sin(rotation/180.*math.pi)
        self.xr = np.zeros(self.nseg)
        self.yr = np.zeros(self.nseg)
        self.zr = np.zeros(self.nseg)
        if self.nseg >= 4:
            self.prefit_rms = 0.
            for i in range(self.nseg):
                self.xr[i] = ct*M1.X[segment_indices[i]] + st*M1.Y[segment_indices[i]] - center
                self.yr[i] =-st*M1.X[segment_indices[i]] + ct*M1.Y[segment_indices[i]]
                self.zr[i] = M1.RESIDUALS[segment_indices[i]]
                self.prefit_rms = self.prefit_rms + M1.RESIDUALS[segment_indices[i]]**2
            self.prefit_rms = math.sqrt(self.prefit_rms/self.nseg)

            # solve for actuator positions
            self.solve(ac)
        else:
            # not enough points for solution
            print('WARNING: Not enough points for solution of Ring %d Seg %d'%(self.ring,self.seg))
            self.actuator_positions = np.zeros(4)
            self.segment_residuals = np.zeros(self.nseg)
            self.rms = 0.0
            
    def initialize_actuator_values(self,META):
        if self.ring == 1:
            index = (self.seg-1)*4
        elif self.ring == 2:
            index = 48 + (self.seg-1)*4
        elif self.ring == 3:
            index = 144 + (self.seg-1)*4
        elif self.ring == 4:
            index = 336 + (self.seg-1)*4
        elif self.ring == 5:
            index = 528 + (self.seg-1)*4
        else:
            print('WARNING: RING %d OUT OF RANGE'%(self.ring))
            index = 0
        self.actuator_values = META.get_actuators()[index:index+4]
        
    def resolve(self,limit):
        # resolve
        idx = np.where(abs(self.segment_residuals)<limit)[0]
        nseg_2 = len(idx)
        # get the model
        ac,rotation,center = self.look_up_model()

        # correct data for segment coordinates
        ct = math.cos(rotation/180.*math.pi)
        st = math.sin(rotation/180.*math.pi)
        xr_2 = self.xr[idx]
        yr_2 = self.yr[idx]
        zr_2 = self.zr[idx]
        if nseg_2 >= 4:
            self.prefit_rms = 0.
            for i in range(nseg_2):
                self.prefit_rms = self.prefit_rms + zr_2[i]**2
            self.prefit_rms = math.sqrt(self.prefit_rms/nseg_2)

            # solve for actuator positions
            npar = 4
            H = np.zeros((nseg_2,npar))
            for i in range(nseg_2):
                for j in range(npar):
                    H[i,j] = (ac[j,0] + ac[j,1]*xr_2[i] + ac[j,2]*yr_2[i] + ac[j,3]*xr_2[i]*yr_2[i])/100.

            HTH = np.dot(H.transpose(),H)
            HTD = np.dot(H.transpose(),zr_2)
            HTHINV = np.linalg.inv(HTH)
            self.actuator_positions = np.dot(HTHINV,HTD)
            segment_model = np.dot(H,self.actuator_positions)
        
            segment_residuals_2 = zr_2 - segment_model
            self.rms = math.sqrt(np.dot(segment_residuals_2.transpose(),segment_residuals_2)/nseg_2)

            # now write over the points
            self.nseg = nseg_2
            self.xr = np.zeros(self.nseg)
            self.yr = np.zeros(self.nseg)
            self.zr = np.zeros(self.nseg)
            self.segment_residuals = np.zeros(self.nseg)
            for i in range(self.nseg):
                self.xr[i] = xr_2[i]
                self.yr[i] = yr_2[i]
                self.zr[i] = zr_2[i]
                self.segment_residuals[i] = segment_residuals_2[i]
        else:
            # not enough points for solution
            print('WARNING: Not enough points for solution of Ring %d Seg %d'%(self.ring,self.seg))
            self.nseg = nseg_2
            self.xr = np.zeros(self.nseg)
            self.yr = np.zeros(self.nseg)
            self.zr = np.zeros(self.nseg)
            self.segment_residuals = np.zeros(self.nseg)
            self.actuator_positions = np.zeros(4)
            self.rms = 0.0


    def update_residuals_array(self,M1):
        # get the indices of the segment data from the residuals list
        segment_indices = M1.select_residual_segment_indices(self.ring,self.seg)
        for i in range(self.nseg):
            M1.RESIDUALS[segment_indices[i]] = self.segment_residuals[i]
        self.solve(ac)
        

    def set_deformation_model(self,model):
        if model == 0:
            self.deformation_model = 0
            self.deformation_model_label = 'FEM'
        elif model == 1:
            self.deformation_model = 1
            self.deformation_model_label = 'HOL'

    def look_up_model(self):
        if(self.ring == 1):
            if self.deformation_model == 0:
                ac = np.array([[2.005e-02, -1.163e-05, +1.659e-05, -1.105e-08],
                               [2.005e-02, -1.163e-05, -1.659e-05, +1.105e-08],
                               [3.018e-02, +1.172e-05, +2.711e-05, +3.058e-09],
                               [3.018e-02, +1.172e-05, -2.711e-05, -3.058e-09]])
            else:
#                ac = np.array([[  1.97418777e-02,  -1.17450938e-05,   4.77202748e-05, -2.83903644e-08],
#                               [  1.97418777e-02,  -1.17450938e-05,  -4.77202748e-05,  2.83903644e-08],
#                               [  2.97285915e-02,   1.15713938e-05,   1.89619795e-05,  7.38065682e-09],
#                               [  2.97285915e-02,   1.15713938e-05,  -1.89619795e-05, -7.38065682e-09]])
# NEW 12/15/17
                ac = np.array([[  1.98077988e-02,  -1.17843124e-05,   4.78796200e-05,  -2.84851641e-08],
                               [  1.98077988e-02,  -1.17843124e-05,  -4.78796200e-05,   2.84851641e-08],
                               [  3.07300823e-02,   1.19612085e-05,   1.96007669e-05,   7.62929490e-09],
                               [  3.07300823e-02,   1.19612085e-05,  -1.96007669e-05,  -7.62929490e-09]])

            rotation = 90. - (self.seg-1)*30.;
            center = 4548.87;
        elif(self.ring == 2):
            if self.deformation_model == 0:
                ac = np.array([[2.233E-002, -1.054E-005, 1.960E-005, -9.892E-009],
                               [2.232E-002, -1.054E-005, -1.961E-005, 9.866E-009],
                               [2.558E-002, 1.058E-005, 2.342E-005, 5.509E-009],
                               [2.558E-002, 1.058E-005, -2.342E-005, -5.509E-009]])
            else:
#                ac = np.array([[  2.23083886e-02,  -1.09188037e-05,   2.93242045e-05, -1.43526831e-08],
#                               [  2.23083886e-02,  -1.09188037e-05,  -2.93242045e-05,  1.43526831e-08],
#                               [  2.59824841e-02,   1.05850273e-05,   1.89653169e-05,  7.72629730e-09],
#                               [  2.59824841e-02,   1.05850273e-05,  -1.89653169e-05, -7.72629730e-09]])
# NEw 12/15/17
                ac = np.array([[  2.31241215e-02,  -1.13180628e-05,   3.03964792e-05,  -1.48775062e-08],
                               [  2.31241215e-02,  -1.13180628e-05,  -3.03964792e-05,   1.48775062e-08],
                               [  2.86578669e-02,   1.16749539e-05,   2.09181511e-05,   8.52186415e-09],
                               [  2.86578669e-02,   1.16749539e-05,  -2.09181511e-05,  -8.52186415e-09]])

            rotation = 90. - 7.5 - (self.seg-1)*15.;
            center = 9104.07;
        elif(self.ring == 3):
            if self.deformation_model == 0:
                ac = np.array([[2.283E-002, -1.036E-005, 2.413E-005, -1.166E-008],
                               [2.283E-002, -1.036E-005, -2.413E-005, 1.166E-008],
                               [2.419E-002, 1.033E-005, 3.046E-005, 8.063E-009],
                               [2.419E-002, 1.033E-005, -3.046E-005, -8.063E-009]])
            else:
#                ac = np.array([[  2.16862040e-02,  -1.09776053e-05,   3.39642976e-05,  -1.71928039e-08],
#                               [  2.16862040e-02,  -1.09776053e-05,  -3.39642976e-05,   1.71928039e-08],
#                               [  2.47250313e-02,   1.05259729e-05,   2.63452651e-05,   1.12157410e-08],
#                               [  2.47250313e-02,   1.05259729e-05,  -2.63452651e-05,  -1.12157410e-08]])
# NEW 12/15/17
                ac = np.array([[  2.40572885e-02,  -1.21778536e-05,   3.76778207e-05,  -1.90725976e-08],
                               [  2.40572885e-02,  -1.21778536e-05,  -3.76778207e-05,   1.90725976e-08],
                               [  2.98325729e-02,   1.27003623e-05,   3.17875044e-05,   1.35326183e-08],
                               [  2.98325729e-02,   1.27003623e-05,  -3.17875044e-05,  -1.35326183e-08]])

            rotation = 90 - 3.75 - (self.seg-1)*7.5;
            center = 13904.47;
        elif(self.ring == 4):
            # using HOL coefficients
            if self.deformation_model == 0:
                ac = np.array([[  2.53606594e-02,  -1.36841732e-05,   2.66113950e-05,  -1.43590484e-08],
                               [  2.53606594e-02,  -1.36841732e-05,  -2.66113950e-05,   1.43590484e-08],
                               [  3.12737732e-02,   1.43489595e-05,   2.53434143e-05,   1.16280061e-08],
                               [  3.12737732e-02,   1.43489595e-05,  -2.53434143e-05,  -1.16280061e-08]])
            else:
#                ac = np.array([[  2.08185395e-02,  -1.12333238e-05,   2.18452671e-05,  -1.17873283e-08],
#                               [  2.08185395e-02,  -1.12333238e-05,  -2.18452671e-05,   1.17873283e-08],
#                               [  2.33488934e-02,   1.07128847e-05,   1.89213074e-05,   8.68143004e-09],
#                               [  2.33488934e-02,   1.07128847e-05,  -1.89213074e-05,  -8.68143004e-09]])
# NEW 12/15/17
                ac = np.array([[  2.53606594e-02,  -1.36841732e-05,   2.66113950e-05,  -1.43590484e-08],
                               [  2.53606594e-02,  -1.36841732e-05,  -2.66113950e-05,   1.43590484e-08],
                               [  3.12737732e-02,   1.43489595e-05,   2.53434143e-05,   1.16280061e-08],
                               [  3.12737732e-02,   1.43489595e-05,  -2.53434143e-05,  -1.16280061e-08]])

            rotation = 90. - 3.75 - (self.seg-1)*7.5;
            center = 18500;
        elif(self.ring == 5):
            # using HOL coefficients
            if self.deformation_model == 0:
                ac = np.array([[  2.75888455e-02,  -1.46828677e-05,   2.21064467e-05,  -1.17651184e-08],
                               [  2.75888455e-02,  -1.46828677e-05,  -2.21064467e-05,   1.17651184e-08],
                               [  3.22232928e-02,   1.54642315e-05,   2.11300281e-05,   1.01404797e-08],
                               [  3.22232928e-02,   1.54642315e-05,  -2.11300281e-05,  -1.01404797e-08]])
            else:
#                ac = np.array([[  2.03735443e-02,  -1.08428624e-05,   1.63249553e-05,    -8.68819100e-09],
#                               [  2.03735443e-02,  -1.08428624e-05,  -1.63249553e-05,     8.68819100e-09],
#                               [  2.14520135e-02,   1.02950032e-05,   1.40668941e-05,     6.75082179e-09],
#                               [  2.14520135e-02,   1.02950032e-05,  -1.40668941e-05,    -6.75082179e-09]])
# new 12/15/17
                ac = np.array([[  2.75888455e-02,  -1.46828677e-05,   2.21064467e-05,  -1.17651184e-08],
                               [  2.75888455e-02,  -1.46828677e-05,  -2.21064467e-05,   1.17651184e-08],
                               [  3.22232928e-02,   1.54642315e-05,   2.11300281e-05,   1.01404797e-08],
                               [  3.22232928e-02,   1.54642315e-05,  -2.11300281e-05,  -1.01404797e-08]])

            rotation = 90. - 3.75 - (self.seg-1)*7.5;
            center = 22875.;
        else:
            print('ILLEGAL RING %d'%(self.ring))
        return(ac,rotation,center)


    def solve(self,ac):
        npar = 4
        H = np.zeros((self.nseg,npar))
        for i in range(self.nseg):
            for j in range(npar):
                H[i,j] = (ac[j,0] + ac[j,1]*self.xr[i] + ac[j,2]*self.yr[i] + ac[j,3]*self.xr[i]*self.yr[i])/100.

        HTH = np.dot(H.transpose(),H)
        HTD = np.dot(H.transpose(),self.zr)
        HTHINV = np.linalg.inv(HTH)
        self.actuator_positions = np.dot(HTHINV,HTD)
        segment_model = np.dot(H,self.actuator_positions)
        self.segment_residuals = self.zr - segment_model
        self.rms = math.sqrt(np.dot(self.segment_residuals.transpose(),self.segment_residuals)/self.nseg)

    def update_residuals_array(self,M1):
        # get the indices of the segment data from the residuals list
        segment_indices = M1.select_residual_segment_indices(self.ring,self.seg)
        for i in range(self.nseg):
            M1.RESIDUALS[segment_indices[i]] = self.segment_residuals[i]

    def dir_turns_frac(self, value, one_turn=1270):
        turns = np.floor( abs(value)/one_turn )
        flats = np.floor( (abs(value)/one_turn - turns)*6. )
        quarters = round( (abs(value)/one_turn - turns - flats/6.)*24. )
        if value < 0:
            dir = ' UP '
        else:
            dir = 'DOWN'
        return([dir,turns,flats,quarters])
    
    def print_screw_turns(self,one_turn=1270.):
        print('-----------------------------------')
        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[2],one_turn)
        print('%d %d EL %6.0f'%(self.ring,self.seg,self.actuator_positions[2]))
        print('%d %2d EL %s %2.0f turns %2.0f flats %2.0f quarters'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[3],one_turn)
        print('%d %d ER %6.0f'%(self.ring,self.seg,self.actuator_positions[3]))
        print('%d %2d ER %s %2.0f turns %2.0f flats %2.0f quarters'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[1],one_turn)
        print('%d %d IR %6.0f'%(self.ring,self.seg,self.actuator_positions[1]))
        print('%d %2d IR %s %2.0f turns %2.0f flats %2.0f quarters'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[0],one_turn)
        print('%d %d IL %6.0f'%(self.ring,self.seg,self.actuator_positions[0]))
        print('%d %2d IL %s %2.0f turns %2.0f flats %2.0f quarters'%(self.ring,self.seg,dir,turns,flats,quarters))
        print('-----------------------------------')

    def write_screw_turns(self,F,one_turn=1270.):
        F.write('----------------------------------\n')

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[2],one_turn)
        F.write('%d %2d EL %s %2.0f turns %2.0f flats %2.0f quarters\n'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[3],one_turn)
        F.write('%d %2d ER %s %2.0f turns %2.0f flats %2.0f quarters\n'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[1],one_turn)
        F.write('%d %2d IR %s %2.0f turns %2.0f flats %2.0f quarters\n'%(self.ring,self.seg,dir,turns,flats,quarters))

        [dir,turns,flats,quarters] = self.dir_turns_frac(self.actuator_positions[0],one_turn)
        F.write('%d %2d IL %s %2.0f turns %2.0f flats %2.0f quarters\n'%(self.ring,self.seg,dir,turns,flats,quarters))

        F.write('----------------------------------\n')

    def print_command_order(self):
        # NOTE: desired actuator position is negative of actuator_positions
        print('%d %d EL %f %f %f'%(self.ring,self.seg,-self.actuator_positions[2],self.actuator_values[0],-self.actuator_positions[2]+self.actuator_values[0]))
        print('%d %d ER %f %f %f'%(self.ring,self.seg,-self.actuator_positions[3],self.actuator_values[1],-self.actuator_positions[3]+self.actuator_values[1]))
        print('%d %d IR %f %f %f'%(self.ring,self.seg,-self.actuator_positions[1],self.actuator_values[2],-self.actuator_positions[1]+self.actuator_values[2]))
        print('%d %d IL %f %f %f'%(self.ring,self.seg,-self.actuator_positions[0],self.actuator_values[3],-self.actuator_positions[0]+self.actuator_values[3]))
        
    def print_fem_order(self):
        # NOTE: desired actuator position is negative of actuator_positions
        print('%d %d EL %f'%(self.ring,self.seg,-self.actuator_positions[2]))
        print('%d %d ER %f'%(self.ring,self.seg,-self.actuator_positions[3]))
        print('%d %d IL %f'%(self.ring,self.seg,-self.actuator_positions[0]))
        print('%d %d IR %f'%(self.ring,self.seg,-self.actuator_positions[1]))
                
    def write_command_actuators(self,F,cmd):
        # NOTE: desired actuator position is negative of actuator_positions
        F.write('%d,%d,EL,%d,%f\n'%(self.ring,self.seg,cmd,-self.actuator_positions[2]+self.actuator_values[0]))
        F.write('%d,%d,ER,%d,%f\n'%(self.ring,self.seg,cmd,-self.actuator_positions[3]+self.actuator_values[1]))
        F.write('%d,%d,IR,%d,%f\n'%(self.ring,self.seg,cmd,-self.actuator_positions[1]+self.actuator_values[2]))
        F.write('%d,%d,IL,%d,%f\n'%(self.ring,self.seg,cmd,-self.actuator_positions[0]+self.actuator_values[3]))

    def write_residuals(self,F):
        for i in range(self.nseg):
            F.write('%d, %d, %f, %f, %d, %f\n'%(self.ring,self.seg,self.xr[i],self.yr[i],self.id[i],self.segment_residuals[i]))
        
