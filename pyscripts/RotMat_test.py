# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:19:13 2016

@author: noel
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def MulMat(a, b, c, nn, n):
    for i in range(0,nn):
        for j in range(0,nn):
            a[0][i + nn * j] = 0
            for k in range(0,nn):
                a[0][i + (nn * j)] = a[0][i + (nn * j)] + b[0][i + (nn * k)] * c[0][(n*9)+(k + (nn * j))];

def BuildStepRmat(mp, ax, ay, az):
    print("mp", mp)
    m1 = np.zeros((1,9))         # double[] m1 = new double[9]        
    m2 = np.zeros((1,9))         # double[] m2 = new double[9]
    c = np.zeros((1,3))          # double[] c = new double[3]         
    s = np.zeros((1,3))          # double[] s = new double[3]
    # double ak, c0c2, c0s2, s0c2, s0s2, t    
    ak = 0
    c0c2 = 0
    c0s2 = 0
    s0c2 = 0
    s0s2 = 0
    t = 0

    for k in range(0,3):   # for(int k = 0; k < 3; k++):
        if k == 0:         #    if(k == 0){  ak = ax
            ak = ax        #        ak = ax
        elif k == 1:       #    }else if(k == 1){   
            ak = ay        #        ak = ay
        elif k == 2:       #    }else if(k == 2){   
            ak = az        #        ak = az}
                           
        t = 0.25*ak*ak     #    t = 0.25 * ak*ak        
        c[0][k] = (1-t)/(1+t) #    c[k] = (1.0-t)/(1.0+t)
        s[0][k] = ak/(1+t)    #    s[k] = ak/(1.0+t)
        #print("c[",k,"] =",c[0][k]," s[",k,"]=",s[0][k])
        
    # The rest of this function is unchanged.
    c0c2 = c[0][0] * c[0][2]
    c0s2 = c[0][0] * s[0][2]
    s0c2 = s[0][0] * c[0][2]
    s0s2 = s[0][0] * s[0][2]
    
    #print("c0c2=",c0c2," c0s2=",c0s2," s0c2=",s0c2," s0s2=",s0s2)
    m1[0][0] = c[0][1] * c[0][2]
    m1[0][1] = s0c2 * s[0][1] + c0s2
    m1[0][2] = -c0c2 * s[0][1] + s0s2
    m1[0][3] = -c[0][1] * s[0][2]
    m1[0][4] = -s0s2 * s[0][1] + c0c2
    m1[0][5] = c0s2 * s[0][1] + s0c2
    m1[0][6] = s[0][1]
    m1[0][7] = -s[0][0] * c[0][1]
    m1[0][8] = c[0][0] * c[0][1]
    m2[0][0] = m1[0][0]
    m2[0][1] = -m1[0][3]
    m2[0][2] = -m1[0][6]
    m2[0][3] = s0c2 * s[0][1] - c0s2
    m2[0][4] = s0s2 * s[0][1] + c0c2
    m2[0][5] = -m1[0][7]
    m2[0][6] = c0c2 * s[0][1] + s0s2
    m2[0][7] = c0s2 * s[0][1] - s0c2
    m2[0][8] = m1[0][8]
    MulMat(mp, m1, m2, 3, 0)

###############################################################################
def genMatrix(TMrMatT, eAngx, eAngy, eAngz):
    p = np.zeros((1,10))
    tq = np.zeros((1,4))
    a1 = 0.5 * eAngy
    a2 = 0.5 * (eAngx - eAngz)
    a3 = 0.5 * (eAngx + eAngz)
   
    q1 = np.sin(a1) * np.cos(a2)
    q2 = np.sin(a1) * np.sin(a2)
    q3 = np.cos(a1) * np.sin(a3)
    q4 = np.cos(a1) * np.cos(a3)
 
    tq[0][0] = q1
    tq[0][1] = q2
    tq[0][2] = q3
    tq[0][3] = q4
    
    k = 0
    k2 = 0
    for i in range(0,4):       
        k1 = k2
        for j in range(i,4):
            p[0][k] = 2*tq[0][k1]*tq[0][k2]
            k1 = k1 + 1
            k = k + 1
        k2 = k2 + 1
        
    ''' The following four lines were slightly modified from java file (removed . from TM.rMat)'''
    TMrMatT[0][0] = p[0][0] + p[0][9] - 1;  TMrMatT[0][4] = p[0][4] + p[0][9] - 1;   TMrMatT[0][8] = p[0][7] + p[0][9] - 1;
    s = 1.0;    #Transpose = 1
    TMrMatT[0][1] = p[0][1] + s * p[0][8];  TMrMatT[0][3] = p[0][1] - s * p[0][8];   TMrMatT[0][2] = p[0][2] - s * p[0][6];
    TMrMatT[0][6] = p[0][2] + s * p[0][6];  TMrMatT[0][5] = p[0][5] + s * p[0][3];   TMrMatT[0][7] = p[0][5] - s * p[0][3];
###############################################################################
coords = [1,0,0]
print(0,coords)
TMrMatT = np.zeros((1,9))

itr = 100

x = np.zeros((1,itr))
y = np.zeros((1,itr))
z = np.zeros((1,itr))
mc = np.zeros((1,9))
mt = np.zeros((1,9))
for i in range(0,itr):
    BuildStepRmat(mc, 0.0, 0.0, 0.1)
    MulMat(mt, mc, TMrMatT, 3, 0)
    TMrMatT[0][0] = mt[0][0]
    TMrMatT[0][1] = mt[0][1]
    TMrMatT[0][2] = mt[0][2]
    TMrMatT[0][3] = mt[0][3]
    TMrMatT[0][4] = mt[0][4]
    TMrMatT[0][5] = mt[0][5]
    TMrMatT[0][6] = mt[0][6]
    TMrMatT[0][7] = mt[0][7]
    TMrMatT[0][8] = mt[0][8]
    
    tx = TMrMatT[0][0]*coords[0] + TMrMatT[0][3]*coords[1] + TMrMatT[0][6]*coords[2]
    ty = TMrMatT[0][1]*coords[0] + TMrMatT[0][4]*coords[1] + TMrMatT[0][7]*coords[2]
    tz = TMrMatT[0][2]*coords[0] + TMrMatT[0][5]*coords[1] + TMrMatT[0][8]*coords[2]
    coords[0] = tx
    coords[1] = ty
    coords[2] = tz
    x[0][i] = coords[0]
    y[0][i] = coords[1]
    z[0][i] = coords[2]

plt.scatter(x,y) 
plt.show()
