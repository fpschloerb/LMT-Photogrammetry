import numpy as np
import csv
import matplotlib.pyplot as pl

def read_target_file(filename):
    f = open(filename,'r')
    r = csv.reader(f)

    rlist = []
    slist = []
    xlist = []
    ylist = []
    alist = []

    for row in r:
        rlist.append(eval(row[0]))
        slist.append(eval(row[1]))
        xlist.append(eval(row[2]))
        ylist.append(eval(row[3]))
        alist.append(eval(row[4]))

    f.close()

    ring = np.array(rlist)
    segment = np.array(slist)
    x = np.array(xlist)
    y = np.array(ylist)
    
    return ring,segment,x,y

def select_segment_points(r,s,x,y,iring,iseg):
    idxr = np.where(r==iring)[0]
    idx = np.where(s[idxr]==iseg)[0]
    return x[idxr[idx]],y[idxr[idx]],len(idxr[idx])



#r57,s57,x57,y57 = read_target_file('M57_TARGET_FILE.csv')
#r58,s58,x58,y58 = read_target_file('M58_TARGET_FILE.csv')
#r59,s59,x59,y59 = read_target_file('M59_TARGET_FILE.csv')
#r60,s60,x60,y60 = read_target_file('M60_TARGET_FILE.csv')

old_r,old_s,old_x,old_y = read_target_file('M1_TARGET_LIST.csv')
new_r,new_s,new_x,new_y = read_target_file('New_M1_Target_List.csv')

nseg = [12,24,48,48,48]

for iring in range(5):
    for iseg in range(nseg[iring]):
        X_OLD,Y_OLD,N_OLD = select_segment_points(old_r,old_s,old_x,old_y,iring+1,iseg+1)
        X_NEW,Y_NEW,N_NEW = select_segment_points(new_r,new_s,new_x,new_y,iring+1,iseg+1)
        if not (N_OLD == N_NEW):
            print('%d %d : %d %d'%(iring+1,iseg+1,N_OLD,N_NEW))
            pl.figure()
            pl.plot(X_NEW,Y_NEW,'b+')
            pl.plot(X_OLD,Y_OLD,'rx')
            pl.axis('equal')
            pl.title('Ring: %d Segment: %d'%(iring+1,iseg+1))
            pl.show()
            pl.close()
        X_NEW = []
        Y_NEW = []
        X_OLD = []
        Y_OLD = []
#        X_57,Y_57,n_57 = select_segment_points(r57,s57,x57,y57,iring+1,iseg+1)
#        X_58,Y_57,n_58 = select_segment_points(r58,s58,x58,y58,iring+1,iseg+1)
#        X_59,Y_57,n_59 = select_segment_points(r59,s59,x59,y59,iring+1,iseg+1)
#        X_60,Y_60,n_60 = select_segment_points(r60,s60,x60,y60,iring+1,iseg+1)
#        if not (n_57==n_58 and n_58==n_59 and n_59==n_60 and n_60==n_57):
#            print('%d %d : %d %d %d %d'%(iring+1,iseg+1,n_57,n_58,n_59,n_60))
