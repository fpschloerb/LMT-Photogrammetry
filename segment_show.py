"""Module with segment_show class for viewing results of fits to segments

Classes: segment_show
Uses:    math,numpy,matplotlib.pyplot
Author:  FPS
Date:    November 6, 2016
Changes:
"""

import numpy as np
import math
import matplotlib.pyplot as pl

class segment_show():
    def __init__(self,segment,scale=0.002,figure=2):
        self.segment = segment
        self.scale = scale
        self.figure = figure
    def set_segment(self,segment):
        self.segment = segment
    def set_scale(self,scale):
        self.scale = scale
    def set_figure(self,figure):
        self.figure = figure
    def make_segment_plot(self):
        pl.figure(num=self.figure,figsize=(6,8),dpi=100)
        pl.rc('xtick', labelsize=10) 
        pl.rc('ytick', labelsize=10) 
        # original data
        pl.subplot(211)
        self.draw_circle_plot(self.segment.zr)
        self.subpanel()
        pl.title('Before fit: RMS=%4.0f microns'%(self.segment.prefit_rms*1000.),fontsize=10)
        pl.axis('equal')
        # residuals after segment fit
        pl.subplot(212)
        self.draw_circle_plot(self.segment.segment_residuals)
        self.subpanel()
        pl.title('After %s fit: [EL=%6.0f ER=%6.0f IR=%6.0f IL=%6.0f] RMS=%4.0f microns'
                 %(self.segment.deformation_model_label,
                   self.segment.actuator_positions[2],
                   self.segment.actuator_positions[3],
                   self.segment.actuator_positions[1],
                   self.segment.actuator_positions[0],
                   self.segment.rms*1000.),fontsize=10)
        pl.axis('equal')
        pl.suptitle('R %d S %d'%(self.segment.ring,self.segment.seg))
    
    def make_single_segment_plot(self):
        pl.figure(num=self.figure,figsize=(6,6),dpi=100)
        pl.rc('xtick', labelsize=10) 
        pl.rc('ytick', labelsize=10) 
        # original data
        self.draw_text_plot(self.segment.zr)
        self.subpanel()
        pl.title('%d %d: RMS=%4.0f microns'%(self.segment.ring,self.segment.seg,self.segment.prefit_rms*1000.),fontsize=10)
        pl.axis('equal')

    def draw_circle_plot(self,zz):
        tt = np.linspace(0.,2.*math.pi,21)
        for i in range(self.segment.nseg):
            if(zz[i]>0.):
                l = np.floor(abs(zz[i])/self.scale)+1
                xp = self.segment.xr[i] + l*np.cos(tt)
                yp = self.segment.yr[i] + l*np.sin(tt)
                pl.plot(xp,yp,'b',lw=1)
            else:
                l = np.floor(abs(zz[i])/self.scale)+1
                xp = self.segment.xr[i] + l*np.cos(tt)
                yp = self.segment.yr[i] + l*np.sin(tt)
                pl.plot(xp,yp,'r',lw=1)
        pl.plot(self.segment.xr,self.segment.yr,'k+',ms=4,lw=1)
        # plot scale circle
        l = np.floor(0.5/self.scale)+1
        xp = -2000 + l*np.cos(tt)
        yp = -1250 + l*np.sin(tt)
        pl.plot(xp,yp,'k',lw=2)

    def draw_text_plot(self,zz):
        for i in range(self.segment.nseg):
            if(zz[i]>0):
                pl.text(self.segment.xr[i],
                        self.segment.yr[i],
                        '%3.1f'%(zz[i]),
                        fontsize=12,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='b')
            else:
                pl.text(self.segment.xr[i],
                        self.segment.yr[i],
                        '%3.1f'%(abs(zz[i])),
                        fontsize=12,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='r')

    def subpanel(self):
        if(self.segment.ring == 1):
            self.subpanel1()
        elif(self.segment.ring == 2):
            self.subpanel2()
        elif(self.segment.ring == 3):
            self.subpanel3()
        elif(self.segment.ring == 4):
            self.subpanel4()
        elif(self.segment.ring == 5):
            self.subpanel5()

    def subpanel1(self):
        # ring 1
        ri = np.array([1625.000, 2867.112, 4136.256, 5389.170])
        ro = np.array([2867.112, 4136.256, 5389.170, 6650.000])
        phi = 30./180.*math.pi
        xoff = 4548.87

        # 1.1
        ang = np.linspace(-phi/2,phi/2,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

        # 1.2
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

        # 1.3
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

        # 1.4
        ang = np.linspace(-phi/2.,-phi/6.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(-phi/6.,phi/6.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(phi/6.,phi/2.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)


    def subpanel2(self):
        ri = np.array([6650.000, 7894.300, 9126.400, 10345.300])
        ro = np.array([7894.300, 9126.400, 10345.300, 11550.000])
        phi = 15./180.*math.pi;
        xoff = 9104.07;
        # 2.1
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 2.2
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 2.3
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 2.4
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

    def subpanel3(self):
        ri = np.array([11550.000, 12746.060, 13928.440, 15096.580])
        ro = np.array([12746.060, 13928.440, 15096.580, 16250.000])
        phi = 30./4./180.*math.pi;
        xoff = 13904.47;

        # 3.1
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 3.2
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 3.3
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 3.4
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

    def subpanel4(self):
        # BOGUS DISTANCES
        ri = np.array([16250., 16250.+1125., 16250.+2*1125., 16250.+3*1125.])
        ro = np.array([16250.+1125., 16250.+2*1125., 16250.+3*1125., 20750.000])
        phi = 30./4./180.*math.pi;
        xoff = 16250.+2*1125.;

        # 4.1
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 4.2
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 4.3
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 4.4
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

    def subpanel5(self):
        # BOGUS DISTANCES
        ri = np.array([20750., 20750.+1062.5, 20750.+2*1062.5, 20750.+3*1062.5])
        ro = np.array([20750.+1062.5, 20750.+2*1062.5, 20750.+3*1062.5, 25000.000])
        phi = 30./4./180.*math.pi;
        xoff = 20750.+2*1062.5;

        # 5.1
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[0]*np.cos(ang)-xoff;
        yi = ri[0]*np.sin(ang);
        xo = ro[0]*np.cos(ang)-xoff;
        yo = ro[0]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 5.2
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[1]*np.cos(ang)-xoff;
        yi = ri[1]*np.sin(ang);
        xo = ro[1]*np.cos(ang)-xoff;
        yo = ro[1]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 5.3
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[2]*np.cos(ang)-xoff;
        yi = ri[2]*np.sin(ang);
        xo = ro[2]*np.cos(ang)-xoff;
        yo = ro[2]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        # 5.4
        ang = np.linspace(-phi/2.,0.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)
        ang = np.linspace(0.,phi/2.,101)
        xi = ri[3]*np.cos(ang)-xoff;
        yi = ri[3]*np.sin(ang);
        xo = ro[3]*np.cos(ang)-xoff;
        yo = ro[3]*np.sin(ang);
        pl.plot(xi,yi,'k',lw=1)
        pl.plot(xo,yo,'k',lw=1)
        pl.plot([xi[0], xo[0]],[yi[0], yo[0]],'k',lw=1)
        pl.plot([xi[len(ang)-1], xo[len(ang)-1]],[yi[len(ang)-1], yo[len(ang)-1]],'k',lw=1)

