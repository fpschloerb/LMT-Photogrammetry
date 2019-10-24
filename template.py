import numpy as np
import math
import matplotlib.pyplot as pl

def template_meters(my_color='k'):
    ring = np.zeros((5,2))
    nring = np.zeros(5,dtype='int')
    ring[0,0] = 1.650
    ring[0,1] = 6.650
    nring[0] = 12
    ring[1,0] = 6.650
    ring[1,1] = 11.550
    nring[1] = 24
    ring[2,0] = 11.550
    ring[2,1] = 16.250
    nring[2] = 48
    ring[3,0] = 16.250
    ring[3,1] = 20.750
    nring[3] = 48
    ring[4,0] = 20.750
    ring[4,1] = 25.000
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

