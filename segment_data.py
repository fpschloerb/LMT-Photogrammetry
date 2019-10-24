import numpy as np
import math
import matplotlib.pyplot as pl
import matplotlib.mlab as mlab
import math
from numpy.random import RandomState
from scipy.special import jn
from template import template_meters

class dish_data():
    def __init__(self,nrings=5):
        self.DISH = []
        self.DISH.append(segment_data(1))
        self.DISH.append(segment_data(2))
        self.DISH.append(segment_data(3))
        self.DISH.append(segment_data(4))
        self.DISH.append(segment_data(5))
        self.nrings = nrings
        self.nseg = np.array([12,24,48,48,48],dtype=int)
        self.area = math.pi*self.DISH[self.nrings-1].ro**2 - math.pi*self.DISH[0].ri**2
        self.ndish = 0
        for i in range(self.nrings):
            self.ndish = self.ndish + self.nseg[i]*self.DISH[i].npoints

    def compute_surface_points(self,act):
        self.act_rms = np.std(act)
        self.xdish = np.zeros(self.ndish)
        self.ydish = np.zeros(self.ndish)
        self.zdish = np.zeros(self.ndish)
        sss = 0
        count = 0
        for i in range(self.nrings):
            ring = i+1
            for j in range(self.nseg[i]):
                seg = j+1
                X,Y = self.DISH[i].get_XY(seg)
                Z = self.DISH[i].get_deformation_model([act[sss,3],act[sss,2],act[sss,0],act[sss,1]])
                for k in range(self.DISH[i].npoints):
                    self.xdish[count] = X[k]
                    self.ydish[count] = Y[k]
                    self.zdish[count] = Z[k]
                    count = count + 1
                sss = sss + 1

    def get_segment_data(self,ring):
        index = ring-1
        return(self.DISH[index])

    def show_surface(self, grid_dimension, resolution, limits):        
        xi,yi = np.mgrid[-grid_dimension:grid_dimension:resolution,-grid_dimension:grid_dimension:resolution]
        myGrid = mlab.griddata(self.xdish/1000.,self.ydish/1000,self.zdish,xi,yi,interp='linear')
        # mask central hole
        for i in range(len(xi)):
            for j in range(len(yi)):
                rrr = np.sqrt(xi[i,j]*xi[i,j]+yi[i,j]*yi[i,j])
                if (rrr < 1.625):
                    myGrid.mask[i,j] = True
        pl.imshow(myGrid,cmap=pl.cm.jet,clim=limits,origin='lower',extent=[-grid_dimension,grid_dimension,-grid_dimension,grid_dimension])
        template_meters()
        pl.axis([-grid_dimension,grid_dimension,-grid_dimension,grid_dimension])
        pl.colorbar()
        pl.title('ACT - RMS=%5.0f microns'%(self.act_rms))
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')
        

class segment_data():
    """segment_data manages the information needed for segment deformation and integrations"""
    def __init__(self,ring):
        self.ring = ring
        self.init_data()
    def init_data(self):
        # note that ac gives Z error for actuator movement
        if self.ring == 1:
            ac = np.array([[  1.98077988e-02,  -1.17843124e-05,   4.78796200e-05,  -2.84851641e-08],
                           [  1.98077988e-02,  -1.17843124e-05,  -4.78796200e-05,   2.84851641e-08],
                           [  3.07300823e-02,   1.19612085e-05,   1.96007669e-05,   7.62929490e-09],
                           [  3.07300823e-02,   1.19612085e-05,  -1.96007669e-05,  -7.62929490e-09]])
            
            self.center = 4548.87
            self.ri = 1650.
            self.ro = 6650.
            self.nradius = 20
            self.dr = (self.ro-self.ri)/self.nradius
            self.ntheta = 40
            self.dtd = 30./self.ntheta
            self.dt = self.dtd/180.*math.pi
            self.rvalues = np.linspace(self.ri+self.dr/2.,self.ro-self.dr/2.,self.nradius)
            self.tvalues = np.linspace(-15.+self.dtd/2.,15.-self.dtd/2.,self.ntheta)/180.*math.pi

            

        elif self.ring == 2:
            ac = np.array([[  2.31241215e-02,  -1.13180628e-05,   3.03964792e-05,  -1.48775062e-08],
                           [  2.31241215e-02,  -1.13180628e-05,  -3.03964792e-05,   1.48775062e-08],
                           [  2.86578669e-02,   1.16749539e-05,   2.09181511e-05,   8.52186415e-09],
                           [  2.86578669e-02,   1.16749539e-05,  -2.09181511e-05,  -8.52186415e-09]])
            
            self.center = 9104.07
            self.ri = 6650.
            self.ro = 11550.
            self.nradius = 20
            self.dr = (self.ro-self.ri)/self.nradius
            self.ntheta = 20
            self.dtd = 15./self.ntheta
            self.dt = self.dtd/180.*math.pi
            self.rvalues = np.linspace(self.ri+self.dr/2.,self.ro-self.dr/2.,self.nradius)
            self.tvalues = np.linspace(-7.5+self.dtd/2.,7.5-self.dtd/2.,self.ntheta)/180.*math.pi

        elif self.ring == 3:
            ac = np.array([[  2.40572885e-02,  -1.21778536e-05,   3.76778207e-05,  -1.90725976e-08],
                           [  2.40572885e-02,  -1.21778536e-05,  -3.76778207e-05,   1.90725976e-08],
                           [  2.98325729e-02,   1.27003623e-05,   3.17875044e-05,   1.35326183e-08],
                           [  2.98325729e-02,   1.27003623e-05,  -3.17875044e-05,  -1.35326183e-08]])
        
            self.center = 13904.47
            self.ri = 11550.
            self.ro = 16250.
            self.nradius = 20
            self.dr = (self.ro-self.ri)/self.nradius
            self.ntheta = 10
            self.dtd = 7.5/self.ntheta
            self.dt = self.dtd/180.*math.pi
            self.rvalues = np.linspace(self.ri+self.dr/2.,self.ro-self.dr/2.,self.nradius)
            self.tvalues = np.linspace(-3.75+self.dtd/2.,3.75-self.dtd/2.,self.ntheta)/180.*math.pi
        elif self.ring == 4:
            ac = np.array([[  2.53606594e-02,  -1.36841732e-05,   2.66113950e-05,  -1.43590484e-08],
                           [  2.53606594e-02,  -1.36841732e-05,  -2.66113950e-05,   1.43590484e-08],
                           [  3.12737732e-02,   1.43489595e-05,   2.53434143e-05,   1.16280061e-08],
                           [  3.12737732e-02,   1.43489595e-05,  -2.53434143e-05,  -1.16280061e-08]])
            
            self.center = 18500
            self.ri = 16250.
            self.ro = 20750.
            self.nradius = 20
            self.dr = (self.ro-self.ri)/self.nradius
            self.ntheta = 10
            self.dtd = 7.5/self.ntheta
            self.dt = self.dtd/180.*math.pi
            self.rvalues = np.linspace(self.ri+self.dr/2.,self.ro-self.dr/2.,self.nradius)
            self.tvalues = np.linspace(-3.75+self.dtd/2.,3.75-self.dtd/2.,self.ntheta)/180.*math.pi
        elif self.ring == 5:
            ac = np.array([[  2.75888455e-02,  -1.46828677e-05,   2.21064467e-05,  -1.17651184e-08],
                           [  2.75888455e-02,  -1.46828677e-05,  -2.21064467e-05,   1.17651184e-08],
                           [  3.22232928e-02,   1.54642315e-05,   2.11300281e-05,   1.01404797e-08],
                           [  3.22232928e-02,   1.54642315e-05,  -2.11300281e-05,  -1.01404797e-08]])
            
            self.center = 22875.
            self.ri = 20750.
            self.ro = 25000.
            self.nradius = 20
            self.dr = (self.ro-self.ri)/self.nradius
            self.ntheta = 10
            self.dtd = 7.5/self.ntheta
            self.dt = self.dtd/180.*math.pi
            self.rvalues = np.linspace(self.ri+self.dr/2.,self.ro-self.dr/2.,self.nradius)
            self.tvalues = np.linspace(-3.75+self.dtd/2.,3.75-self.dtd/2.,self.ntheta)/180.*math.pi
        else:
            print('init_data: illegal ring ')

        self.npoints = self.nradius*self.ntheta
        self.xr = np.zeros(self.npoints)
        self.yr = np.zeros(self.npoints)
        self.rr = np.zeros(self.npoints)
        self.tr = np.zeros(self.npoints)
        k = 0;
        for r in self.rvalues:
            for t in self.tvalues:
                self.xr[k] = r*math.cos(t) - self.center
                self.yr[k] = r*math.sin(t)
                self.rr[k] = r
                self.tr[k] = t
                k = k+1
        self.phase = 1./(1.+self.rr**2/4./17500.**2)
        npar = 4
        self.H = np.zeros((self.npoints,npar))
        for i in range(self.npoints):
            for j in range(npar):
                self.H[i,j] = (ac[j,0]+ac[j,1]*self.xr[i]+ac[j,2]*self.yr[i]+ac[j,3]*self.xr[i]*self.yr[i])/100.


    def get_rotation(self,seg):
        if self.ring == 1:
            rotation = (90. - (seg-1)*30.)/180.*math.pi
        elif self.ring == 2:
            rotation = (90. - 7.5 - (seg-1)*15.)/180.*math.pi
        elif self.ring == 3:
            rotation = (90 - 3.75 - (seg-1)*7.5)/180.*math.pi
        elif self.ring == 4:
            rotation = (90. - 3.75 - (seg-1)*7.5)/180.*math.pi
        elif self.ring == 5:
            rotation = (90. - 3.75 - (seg-1)*7.5)/180.*math.pi
        else:
            print('rotation: illegal ring')
            rotation = 0.
        return(rotation)

    def get_XY(self,seg):
        rotation = self.get_rotation(seg)
        ct = math.cos(-rotation)
        st = math.sin(-rotation)

        X =  ct*(self.xr+self.center) + st*self.yr
        Y = -st*(self.xr+self.center) + ct*self.yr
        return(X,Y)

    def get_deformation_model(self,actuator_positions):
        model = np.dot(self.H,actuator_positions)
        return(model)

    def integrate_segment(self,bu,bv,wavelength,seg,actuator_positions):
        tpw = 2.*np.pi/wavelength
        rotation = self.get_rotation(seg)
        model = self.get_deformation_model(actuator_positions)
        phase = tpw*model*self.phase #make this path error not z error
        dreal = np.cos(phase)
        dimag = np.sin(phase)
        cpd = np.cos(tpw*self.rr*(bu*np.cos(self.tr+rotation)+bv*np.sin(self.tr+rotation)))
        spd = np.sin(tpw*self.rr*(bu*np.cos(self.tr+rotation)+bv*np.sin(self.tr+rotation)))
        rpart = np.sum((dreal*cpd-dimag*spd)*self.rr)*self.dr*self.dt
        ipart = np.sum((dreal*spd+dimag*cpd)*self.rr)*self.dr*self.dt
        return(rpart,ipart)

