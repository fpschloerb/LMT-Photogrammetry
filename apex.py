"""Module with apex class for solving photogrammetry geometry from apex targets

Classes: apex
Uses:    math,numpy,scipy.optimize
Author:  FPS
Date:    November 6, 2016
Changes: December 2017 - added speedup to solver using 
                         scipy.optimize.leastsq
Parameters:
nref         __init__ number of reference targets [4]
apex_pt[4,3] __init__ true locations of apex reference points [4]
reference_dimensions[nref,nref] 
             __init__ distances between reference point pairs
data[4,3]    __init__ photogrammetry points [4]
data_dimensions[nref,nref] 
             __init__ distances between all point pairs
dref[nref,3] __init__ working vector with rotated reference locations
ddif[nref,3] __init__ working vector with difference between dref and apex_pt
model[6]     search   array with geometry model parameters
error[6]     search   array with geometry model parameter errors
rms          search   rms of fit
"""
import math
import numpy as np
from scipy.optimize import leastsq


def compute_the_residuals(v,data,ref):
    residuals = np.zeros(12)

    phi = v[3]*np.pi/180.
    theta = v[4]*np.pi/180.
    psi = v[5]*np.pi/180.
    # rotate function
    D = np.array([[math.cos(phi),math.sin(phi),0.],
                  [-math.sin(phi),math.cos(phi),0.],
                  [0.,0.,1.]])
    C = np.array([[math.cos(theta),0.,-math.sin(theta)],
                  [0.,1.,0.],
                  [math.sin(theta),0.,math.cos(theta)]])
    B = np.array([[1.,0.,0.],
                  [0.,math.cos(psi),math.sin(psi)],
                  [0.,-math.sin(psi),math.cos(psi)]])
    R = np.dot(B,np.dot(C,D))
    
    k = 0
    for i in range(4):
        d = np.dot(R,data[i,:]) + v[0:3]
        for j in range(3):
            residuals[k] = d[j] - ref[i,j]
            k = k+1
    return(residuals)

class apex():
    """apex solves for photogrammetry geometry"""
    def __init__(self,
                 filename,redfile,solve=True,
                 log_model_file='APEX_MODELS.csv',
                 log_dim_file='APEX_DIMENSIONS.csv'):

        # this is the location of the reference targers                                 
        self.nref = 4
        self.apex_pt = np.array([
                [1097.78475,1019.3272,245.364655],
                [1009.0703,-1108.9863,233.06653],
                [-1139.65305,-983.49458,224.26834],
                [-1024.4449,1109.39855,223.30923]])
        self.reference_dimensions = self.check_dimensions(self.apex_pt)

        # set up the data array
        input_data = np.loadtxt(filename)
        self.data = np.zeros((self.nref,3))
        self.data[:,0] = input_data[:,1] # x is second element
        self.data[:,1] = input_data[:,2] # y is third element
        self.data[:,2] = input_data[:,0] # z is first element
        self.data_dimensions = self.check_dimensions(self.data)

        self.dref = np.zeros((self.nref,3))
        self.ddif = np.zeros((self.nref,3))

        # print comparison of the distances between points
        print('Reference and Data Dimension Difference (mm)')
        print(self.reference_dimensions-self.data_dimensions)

        # solve for model if needed, otherwise read from reduced file
        if solve:
            # we now (Dec 2017) use the scipy function to search for minimum 
            # this procedure speeds things up considerably!            
            v0 = np.array([111.60, 8.637, 1321.0, -89.3435, 0.5172, 0.866])
            theData = self.data
            theApexRef = self.apex_pt
            fit = leastsq(compute_the_residuals, v0, args=(theData,theApexRef))
            # then once minimum is found, we iterate on a grid to find 
            # minimum and errors
            #print(fit)
            self.search(np.array(fit[0]),
                        np.array([0.005, 0.005, 0.005, 0.0001, 0.0001, 0.0001])
                        ,4)
            self.print_solution()
            self.write_model(redfile)
        else:
            self.read_model(redfile)
            self.print_solution()
            
        f = open(log_model_file,'a')
        self.write_csv_model_line(f,filename)
        f.close()

        f = open(log_dim_file,'a')
        self.write_csv_dimension_line(f,filename)
        f.close()


    def check_dimensions(self,a):
        al = np.zeros((self.nref,self.nref))
        for i in range(self.nref):
            for j in range(self.nref):
                al[i,j] = np.sqrt((a[i,0]-a[j,0])**2
                                  +(a[i,1]-a[j,1])**2
                                  +(a[i,2]-a[j,2])**2
                                  )
        return al

    def rotate(self,v,data):
        phi = v[3]*np.pi/180.
        theta = v[4]*np.pi/180.
        psi = v[5]*np.pi/180.

        D = np.array([[math.cos(phi),math.sin(phi),0.],
                      [-math.sin(phi),math.cos(phi),0.],
                      [0.,0.,1.]])
        C = np.array([[math.cos(theta),0.,-math.sin(theta)],
                      [0.,1.,0.],
                      [math.sin(theta),0.,math.cos(theta)]])
        B = np.array([[1.,0.,0.],
                      [0.,math.cos(psi),math.sin(psi)],
                      [0.,-math.sin(psi),math.cos(psi)]])
        R = np.dot(B,np.dot(C,D))

        return(np.dot(R,data) + v[0:3])

    def compute_residuals(self,v):
        residuals = np.zeros(self.nref)
        for i in range(self.nref):
            self.dref[i,:] = self.rotate(v,self.data[i,:])
            self.ddif[i,:] = self.dref[i,:] - self.apex_pt[i,:];
            residuals[i] = math.sqrt(self.ddif[i,0]**2 + self.ddif[i,1]**2 + self.ddif[i,2]**2)
        return(residuals)

    def compute_chisq(self,v):
        chisq = 0.
        for i in range(self.nref):
            self.dref[i,:] = self.rotate(v,self.data[i,:])
            self.ddif[i,:] = self.dref[i,:] - self.apex_pt[i,:];
            chisq = chisq + self.ddif[i,0]**2 + self.ddif[i,1]**2 + self.ddif[i,2]**2
        return(chisq)

    def search(self,initial_model,delta,niterations):
        tmodel = initial_model
        P_variance = np.zeros(6)
        P_error = np.zeros(6)
        for iteration in range(niterations):
            #print('iteration: ',iteration)
            for ax in range(6):
                r_0 = self.compute_chisq(tmodel)
                x_0 = tmodel[ax]
                tmodel[ax] = tmodel[ax] + delta[ax]
                r_up = self.compute_chisq(tmodel)
                x_up = tmodel[ax]
                if(r_up>r_0):
                    tmodel[ax] = tmodel[ax] - 2*delta[ax]
                    r_down =self.compute_chisq(tmodel)
                    x_down = tmodel[ax]
                    if(r_down<r_0):
                        while(r_down<r_0):
                            tmodel[ax] = tmodel[ax]-delta[ax]
                            r_up = r_0
                            x_up = x_0
                            r_0 = r_down
                            x_0 = x_down
                            r_down = self.compute_chisq(tmodel)
                            x_down = tmodel[ax]
                        tmodel[ax] = tmodel[ax] + delta[ax]
                        P = np.polyfit(np.array([x_down-x_0, 0., x_up-x_0]),np.array([r_down, r_0, r_up]),2)
                        tmodel[ax] = x_0-P[1]/(2*P[0])
                        P_variance[ax] = np.polyval(P,-P[1]/(2.*P[0]))/(self.nref*3)/P[0]
                        #print('A %d %f %f %f %f %f %f %f\n'%(ax,x_down,x_0,x_up,r_down,r_0,r_up,x_0-P[1]/(2*P[0])))
                    else:
                        tmodel[ax] = tmodel[ax] + delta[ax]
                        P = np.polyfit([x_down-x_0, 0., x_up-x_0],[r_down, r_0, r_up],2)
                        #print('B %d %f %f %f %f %f %f %f\n'%(ax,x_down,x_0,x_up,r_down,r_0,r_up,x_0-P[1]/(2*P[0])))
                        tmodel[ax] = x_0-P[1]/(2*P[0])
                        P_variance[ax] = np.polyval(P,-P[1]/(2.*P[0]))/(self.nref*3)/P[0]
                else:
                    while(r_up<r_0):
                        tmodel[ax] = tmodel[ax] + delta[ax]
                        r_down = r_0
                        x_down = x_0
                        r_0 = r_up
                        x_0 = x_up
                        r_up = self.compute_chisq(tmodel)
                        x_up = tmodel[ax]
                    tmodel[ax] = tmodel[ax] - delta[ax]
                    P = np.polyfit([x_down-x_0, 0., x_up-x_0],[r_down, r_0, r_up],2)
                    #print('C %d %f %f %f %f %f %f %f\n'%(ax,x_down,x_0,x_up,r_down,r_0,r_up,x_0-P[1]/(2*P[0])))
                    tmodel[ax] = x_0-P[1]/(2*P[0]);
                    P_variance[ax] = np.polyval(P,-P[1]/(2.*P[0]))/(self.nref*3)/P[0]
                #print('%d %d %f %f\n'%(iteration,ax,tmodel[ax],np.polyval(P,-P[1]/(2.*P[0]))))
            delta = delta/2.;
            
        self.model = tmodel
        self.error = np.sqrt(P_variance)
        self.rms  = math.sqrt(self.compute_chisq(tmodel)/(3.*self.nref));

    def print_solution(self):
        # print results
        print('TRANSLATION:')
        print('X0    = %8.3f (%8.3f) millimeters'%(self.model[0],self.error[0]))
        print('Y0    = %8.3f (%8.3f) millimeters'%(self.model[1],self.error[1]))
        print('Z0    = %8.3f (%8.3f) millimeters'%(self.model[2],self.error[2]))
        print('ROTATION:')
        print('Z Axis= %10.6f (%10.6f) degrees'%(self.model[3]%360.,self.error[3])) # mod 360
        print('Y Axis= %10.6f (%10.6f) degrees'%(self.model[4]%360.,self.error[4]))
        print('X Axis= %10.6f (%10.6f) degrees'%(self.model[5]%360.,self.error[5]))
        print('RMS   = %f millimeters'%(self.rms))

    def write_model(self, filename):
        f = open(filename,'w')
        for i in range(6):
            f.write('%d %f %f\n'%(i,self.model[i],self.error[i]))
        f.write('%f\n'%(self.rms))
        z = self.reference_dimensions-self.data_dimensions
        f.write('%f %f %f %f %f %f\n'%(z[0,1],z[0,2],z[0,3],z[1,2],z[1,3],z[2,3]))
        f.close()

    def read_model(self,filename):
        f = open(filename,'r')
        s = f.readlines()
        self.model = np.zeros(6)
        self.error = np.zeros(6)
        self.rms = 0.
        for i in range(6):
            fields = s[i].split()
            self.model[i] = eval(fields[1])
            self.error[i] = eval(fields[2])
        fields = s[6].split()
        self.rms = eval(fields[0])

    def print_csv_dimension_line(self,id):
        z = self.reference_dimensions-self.data_dimensions
        print('%s,%f,%f,%f,%f,%f,%f'%(id,z[0,1],z[0,2],z[0,3],z[1,2],z[1,3],z[2,3]))
        print(' ')

    def write_csv_dimension_line(self,f,id):
        z = self.reference_dimensions-self.data_dimensions
        f.write('%s,%f,%f,%f,%f,%f,%f\n'%(id,z[0,1],z[0,2],z[0,3],z[1,2],z[1,3],z[2,3]))

    def print_csv_model_line(self,id):
        print('%s,%8.3f,%8.3f,%8.3f,%10.6f,%10.6f,%10.6f,%f'%(id,
                                                              self.model[0],
                                                              self.model[1],
                                                              self.model[2],
                                                              self.model[3],
                                                              self.model[4],
                                                              self.model[5],
                                                              self.rms))
        print(' ')

    def write_csv_model_line(self,f,id):
        f.write('%s,%f,%f,%f,%f,%f,%f,%f\n'%(id,
                                                              self.model[0],
                                                              self.model[1],
                                                              self.model[2],
                                                              self.model[3],
                                                              self.model[4],
                                                              self.model[5],
                                                              self.rms))
