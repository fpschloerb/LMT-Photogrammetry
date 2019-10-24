"""Module with m1_show class which plots m1_photo data

Classes: m1_show
Uses:    math,numpy,m1_photo
Author:  FPS
Date:    November 6, 2016
Changes: 11/21/17: added segment dimensions for rings 4 and 5
         11/21/17: added blackout for segments with no stickers
         12/10/17: reversed jet for image so that blue is UP
         12/17/17: blackout subpanels at tetrapod
         10/25/18: mlab.griddata deprecated - use scipy.interpolate.griddata
"""
import numpy as np
import math
import matplotlib.pyplot as pl
#from scipy.interpolate import griddata
import matplotlib.mlab as mlab

class m1_show():
    """m1_show makes plots of products of m1_photo"""
    def __init__(self,m1):
        self.M1 = m1
        pl.ion()
    def data_xy(self,figure=1):
        pl.figure(figure)
        pl.plot(self.M1.X,self.M1.Y,'b.')
    def circle_plot_m1_residuals(self,mytitle,scale,figure=1,extent=[-26000.,26000.,-26000.,26000.],scale_bar=1.):
        pl.ioff()
        pl.figure(num=figure, figsize=(8,8), dpi=100)
        my_width = 1
        my_marker_size = 4
        tt = np.linspace(0.,2.*math.pi,21)
        if extent[1]>25000.:
            scale_location = 23000.
        else:
            scale_location = 15000.

        # compute RMS for title
        rms = np.sqrt(np.dot(self.M1.RESIDUALS.transpose(),self.M1.RESIDUALS)/self.M1.fit_nobs)
        for i in range(self.M1.fit_nobs):
            if(self.M1.RESIDUALS[i] > 0.):
                # positive residuals plotted with BLUE circle
                l = np.floor(abs(self.M1.RESIDUALS[i]/scale)+1);
                xp = self.M1.X[i] + l*np.cos(tt);
                yp = self.M1.Y[i] + l*np.sin(tt);
                pl.plot(xp,yp,'b',lw=my_width);
            else:
                # negative residuals plotted with RED circle
                l = np.floor(abs(self.M1.RESIDUALS[i]/scale)+1);
                xp = self.M1.X[i] + l*np.cos(tt);
                yp = self.M1.Y[i] + l*np.sin(tt);
                pl.plot(xp,yp,'r',lw=my_width);

        # point locations are plotted with BLACK dots
        pl.plot(self.M1.X,self.M1.Y,'k+',ms=my_marker_size,lw=my_width);
        # plot a scale circle
        l = np.floor(abs(scale_bar/scale)+1);
        xp = -scale_location + l*np.cos(tt);
        yp = -scale_location + l*np.sin(tt);
        pl.plot(xp,yp,'k');
        pl.text(-scale_location,-scale_location,'%4.0f'%(scale_bar*1000))
        # draw segment boundaries
        self.template(figure=figure,my_color='k')
        # make axes have equal dimensions for proper aspect ratio
        pl.axis('square')
        pl.axis(extent)
        # label plot
        pl.xlabel('X Panel (mm)')
        pl.ylabel('Y Panel (mm)')
        pl.title('%s RMS=%6.3f'%(mytitle,rms))
        pl.ion()
        pl.show()
        
    def template(self,figure=1,my_color='k'):
        ring = np.zeros((5,2))
        nring = np.zeros(5,dtype='int')
        ring[0,0] = 1650.
        ring[0,1] = 6650
        nring[0] = 12
        ring[1,0] = 6650.
        ring[1,1] = 11550.
        nring[1] = 24
        ring[2,0] = 11550.
        ring[2,1] = 16250.
        nring[2] = 48
        ring[3,0] = 16250.
        ring[3,1] = 20750.
        nring[3] = 48
        ring[4,0] = 20750.
        ring[4,1] = 25000.
        nring[4] = 48

        t = np.linspace(0.,2.*math.pi,1001)
        my_width = 1
        x1 = ring[0,0]*np.cos(t);
        y1 = ring[0,0]*np.sin(t);
        pl.plot(x1,y1,my_color,lw=my_width)
        x2 = ring[1,0]*np.cos(t);
        y2 = ring[1,0]*np.sin(t);
        pl.plot(x2,y2,my_color,lw=my_width)
        x3 = ring[2,0]*np.cos(t);
        y3 = ring[2,0]*np.sin(t);
        pl.plot(x3,y3,my_color,lw=my_width)
        x4 = ring[3,0]*np.cos(t);
        y4 = ring[3,0]*np.sin(t);
        pl.plot(x4,y4,my_color,lw=my_width)
        x5 = ring[4,0]*np.cos(t);
        y5 = ring[4,0]*np.sin(t);
        pl.plot(x5,y5,my_color,lw=my_width)
        x6 = ring[4,1]*np.cos(t);
        y6 = ring[4,1]*np.sin(t);
        pl.plot(x6,y6,my_color,lw=my_width)

        for ir in range(5):
            for i in range(nring[ir]):
                tt = (i-1)*2.*math.pi/nring[ir];
                if(ir == 0):
                    tt = tt - 2.*math.pi/nring[ir]/2.;
                pl.plot([ring[ir,0]*math.cos(tt), ring[ir,1]*math.cos(tt)], [ring[ir,0]*math.sin(tt), ring[ir,1]*math.sin(tt)],my_color,lw=my_width)

    def imshow(self, mytitle, figure, rlimits=[-.5,.5], extent=[-26.,26.,-26.,26.], resolution= 0.2, show_points=True, show_segments=True, show_subpanels=True):
        # create the grid
        xi,yi = np.mgrid[extent[0]:extent[1]:resolution,
                         extent[2]:extent[3]:resolution]

        # open figure
        pl.figure(figure, figsize=(8,8),dpi=100)
        # fixed facecolor since axisbg was deprecated 10/2/2018
        pl.subplot('111',facecolor='black')

        # templates
        if show_segments:
            self.show_segments()
        if show_points:
            self.show_points()
        if show_subpanels:
            self.show_subpanels()

        # grid the data to rectangular defined in meters 
        myGrid = mlab.griddata(self.M1.X/1000.,self.M1.Y/1000.,self.M1.RESIDUALS,xi,yi,interp='linear')

        # mask center hole in the image
        for i in range(len(xi)):
            for j in range(len(yi)):
                rrr = np.sqrt(xi[i,j]**2+yi[i,j]**2)
                ttt = math.atan2(xi[i,j],yi[i,j])/math.pi*180.
                if(rrr<1.625):
                    myGrid.mask[i,j] = True
                elif (rrr < 6.650): # ring 1
                    if(ttt<-15.):
                        ttt = ttt + 360.
                    sss = int(np.floor((ttt+15)/30.))
                    if(self.M1.found_segment_1[sss] == 0):
                        #print('Ring 1 Seg %d'%(sss+1))
                        myGrid.mask[i,j] = True
                elif (rrr < 11.550): # ring 2
                    if (ttt<0.):
                        ttt = ttt + 360.
                    sss = int(np.floor(ttt/15.))
                    if(self.M1.found_segment_2[sss] == 0):
                        #print('Ring 2 Seg %d'%(sss+1))
                        myGrid.mask[i,j] = True
                elif (rrr < 16.250): # ring 3
                    if (ttt<0.):
                        ttt = ttt + 360.
                    sss = int(np.floor(ttt/7.5))
                    if(self.M1.found_segment_3[sss] == 0):
                        #print('Ring 3 Seg %d'%(sss+1))
                        myGrid.mask[i,j] = True
                elif (rrr < 20.750): # ring 4
                    if(ttt<0.):
                        ttt = ttt + 360.
                    sss = int(np.floor(ttt/7.5))
                    if(self.M1.found_segment_4[sss] == 0):
                        #print('Ring 4 Seg %d'%(sss+1))
                        myGrid.mask[i,j] = True
                    # the following masks the cover plates at tetrapod
                    if sss == 5: # segment number - 1
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt>41.25): 
                            myGrid.mask[i,j] = True
                    elif sss == 6:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt<48.75): 
                            myGrid.mask[i,j] = True
                    elif sss == 17:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt>131.25): 
                            myGrid.mask[i,j] = True
                    elif sss == 18:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt<138.75): 
                            myGrid.mask[i,j] = True
                    elif sss == 29:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt>221.25): 
                            myGrid.mask[i,j] = True
                    elif sss == 30:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt<228.75): 
                            myGrid.mask[i,j] = True
                    elif sss == 41:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt>311.25): 
                            myGrid.mask[i,j] = True
                    elif sss == 42:
                        if (rrr > 17.325) and (rrr < 19.625) and (ttt<318.75): 
                            myGrid.mask[i,j] = True
                elif (rrr < 25.000): # ring 5
                    if(ttt<0.):
                        ttt = ttt + 360.
                    sss = int(np.floor(ttt/7.5))
                    if(self.M1.found_segment_5[sss] == 0):
                        #print('Ring 5 Seg %d'%(sss+1))
                        myGrid.mask[i,j] = True
                else:
                    myGrid.mask[i,j] = True

        # display the image
        pl.imshow(myGrid.T,cmap=pl.cm.jet_r,clim=rlimits,origin='lower',extent=extent)
        # add labels and colorbar
        pl.colorbar()
        pl.axis(extent)
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')
        pl.title(mytitle)


    def show_points(self):
        # draw the points on the image
        pl.plot(self.M1.X/1000.,self.M1.Y/1000.,'k+',markersize=1)

    def show_segments(self, show_rings=5):
        # draw the segment template in the image
        lw = 1
        zzz = np.linspace(0.,2.*np.pi,1001)
        ring=[[1.625,6.650],[6.550,11.550],[11.550,16.250],[16.250, 20.750],[20.750,25.000]]
        nring = [12,24,48,48,48]
        pl.plot(ring[0][0]*np.cos(zzz),ring[0][0]*np.sin(zzz),'k',linewidth=lw)
        for i in range(show_rings):
            pl.plot(ring[i][1]*np.cos(zzz),ring[i][1]*np.sin(zzz),'k',linewidth=lw)
        for ir in range(show_rings):
            for i in range(nring[ir]):
                tt = (i-1)*2.*np.pi/nring[ir];
                if(ir == 0):
                    tt = tt - 2.*np.pi/nring[ir]/2.;
                pl.plot([ring[ir][0]*np.cos(tt), ring[ir][1]*np.cos(tt)], [ring[ir][0]*np.sin(tt), ring[ir][1]*np.sin(tt)],'k',linewidth=lw)

    def show_subpanels(self):
        zzz = np.linspace(0.,2.*np.pi,1001)
        lw = 1
        ri = np.array([1.625, 2.867112, 4.136256, 5.389170])
        ro = np.array([2.867112, 4.136256, 5.389170, 6.650000])
        phi = 30./180.*np.pi
        for i in range(4):
            pl.plot(ro[i]*np.cos(zzz),ro[i]*np.sin(zzz),'k',linewidth=0.5*lw)
            if i>0:
                for j in range(12):
                    if i==3:
                        ang = -phi/2.+phi/3. + j*phi
                        pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
                        ang = -phi/2.+2*phi/3. + j*phi
                        pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
                    else:
                        ang = j*phi
                        pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
        # ring 2 subpanels
        ri = np.array([6.650000, 7.894300, 9.126400, 10.345300])
        ro = np.array([7.894300, 9.126400, 10.345300, 11.550000])
        phi = 15./180.*math.pi;
        for i in range(4):
            pl.plot(ro[i]*np.cos(zzz),ro[i]*np.sin(zzz),'k',linewidth=0.5*lw)
            for j in range(24):
                ang = phi/2. + j*phi
                pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
        # ring 3
        ri = np.array([11.550000, 12.746060, 13.928440, 15.096580])
        ro = np.array([12.746060, 13.928440, 15.096580, 16.250000])
        phi = 30./4./180.*math.pi;
        for i in range(4):
            pl.plot(ro[i]*np.cos(zzz),ro[i]*np.sin(zzz),'k',linewidth=0.5*lw)
            for j in range(48):
                ang = phi/2. + j*phi
                pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
        # ring 4 - bogus distances
        ri = np.array([16.250000, 17.375000, 18.500000, 19.625000])
        ro = np.array([17.375000, 18.500000, 19.625000, 20.750000])
        phi = 30./4./180.*math.pi;
        for i in range(4):
            pl.plot(ro[i]*np.cos(zzz),ro[i]*np.sin(zzz),'k',linewidth=0.5*lw)
            for j in range(48):
                ang = phi/2. + j*phi
                pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)
        # ring 5 - bogus distances
        ri = np.array([20.750000, 21.821500, 22.875000, 23.937500])
        ro = np.array([21.821500, 22.875000, 23.937500, 25.000000])
        phi = 30./4./180.*math.pi;
        for i in range(4):
            pl.plot(ro[i]*np.cos(zzz),ro[i]*np.sin(zzz),'k',linewidth=0.5*lw)
            for j in range(48):
                ang = phi/2. + j*phi
                pl.plot([ri[i]*np.cos(ang),ro[i]*np.cos(ang)],[ri[i]*np.sin(ang),ro[i]*np.sin(ang)],'k',linewidth=0.5*lw)

