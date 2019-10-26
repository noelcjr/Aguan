# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 16:47:04 2016

@author: noel
"""
import numpy as np
import pandas as pd
# Two Waters
W1 = np.array([[7.82315,8.56378,2.99730],[7.56563,7.87379,2.38589],[7.89423,8.12090,3.84289]])
W2 = np.array([[6.10071,7.84564,4.29736],[5.17240,7.63135,4.38982],[6.55182,7.00766,4.39998]])
qW1 = np.array([-0.834,0.417,0.417])
qW2 = np.array([-0.834,0.417,0.417])
eW1 = np.array([-0.152100,-0.046000,-0.046000])
eW2 = np.array([-0.152100,-0.046000,-0.046000])
rW1 = np.array([1.768200,0.224500,0.224500])
rW2 = np.array([1.768200,0.224500,0.224500])
qO = -0.834
O_Eps = -0.152100     
O_RminH = 1.768200
#####   Hydrogen Parameters ############
qH = 0.417
H_Eps = -0.046000
H_RminH = 0.224500
##################################################
# http://localscf.com/localscf.com/LJPotential.aspx.html
k = 332.0716
VWE = 0.0
EEE = 0.0
num_molec = 3
fW1 = np.zeros((num_molec,3))
fW2 = np.zeros((num_molec,3))
###############   Calculations  ##################
for i in range(num_molec):
    for j in range(num_molec):
        d = np.linalg.norm(W1[i]-W2[j])
        dx = W1[i][0]-W2[j][0]
        dy = W1[i][1]-W2[j][1]
        dz = W1[i][2]-W2[j][2]
        dd = d*d
        Rmin = rW1[i] + rW2[j]
        EPS = np.sqrt(eW1[i]*eW2[j])
        A = Rmin/d
        A2 = A*A
        A6 = A2*A2*A2
        A12 = A6*A6
        VWE += EPS*(A12-(2*A6))
        EEE += k*qW1[i]*qW2[j]/d
        eef = -k*qW1[i]*qW2[j]/d
        vwf = 12*EPS*(A6-A12)
        fW1[i][0] += eef*dx/dd + vwf*dx/dd
        fW1[i][1] += eef*dy/dd + vwf*dy/dd
        fW1[i][2] += eef*dz/dd + vwf*dz/dd
        fW2[j][0] -= eef*dx/dd + vwf*dx/dd
        fW2[j][1] -= eef*dy/dd + vwf*dy/dd
        fW2[j][2] -= eef*dz/dd + vwf*dz/dd
####################################################################################
##############   Results     #######################################################
# This script matches results of a CHARMM simulation on the following directory    #
# /home/noel/Projects/FixMD/Aguan/tests/validation/charmm/gen_2_waters.out         #
##############    FORCES     #######################################################
#    1    1 TIP3 OH2 -126.65564 -45.37471  94.38930 WAT  2172   2.20000 #
#    2    1 TIP3 H1     3.50931  -3.76788  -3.84990 WAT  2172   0.00000 #
#    3    1 TIP3 H2    11.65112  -6.92971  -1.28434 WAT  2172   0.00000 #
#    4    2 TIP3 OH2  103.14963  56.09669 -86.85922 WAT  2701   2.20000 #
#    5    2 TIP3 H1     1.57826  -1.77955   0.04244 WAT  2701   0.00000 #
#    6    2 TIP3 H2     6.76732   1.75516  -2.43827 WAT  2701   0.00000 #
#########################################################################
import matplotlib.pyplot as plt

ro = 0.315061
ep = 0.6
#= (2^(1/6))*TM.ro;
#wallCutOff = 1.122462*ro
regionX = 1.8662
x = [float(i)/100 for i in range(200,240,1)]
fcValx = []
uSum = []
Rmin = rW1[0] + rW2[0]
EPS = np.sqrt(eW1[0]*eW2[0])
for i in x:
    d = i
    dd = d*d
    A = Rmin/d
    A2 = A*A
    A6 = A2*A2*A2
    A12 = A6*A6
    uSum.append(EPS*(A12-(2*A6)))
    fcValx.append(12*EPS*(A6-A12)*d/d)
    
DF = pd.DataFrame(columns=['d','E','E2','F'])
DF['d'] = pd.Series(x)
DF['E'] = pd.Series(uSum)
DF['F'] = pd.Series(fcValx)
plt.figure()
plt.plot(DF['d'],DF['E'])
#plt.plot(DF['d'],DF['F'],color='red')
plt.axvline(x=np.linalg.norm(W1[0]-W2[0]))
plt.axhline(y=0.0)
plt.title("Lennard-Jones potential for two oxygens")
plt.xlabel('distance')
plt.ylabel('energy')
plt.show()
##################################################################################################
line = '{:6}{:>5} {:4}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
line = '{:6}{:>5} {:4}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
lines = []
lines.append(line.format('ATOM',str(1),'OH2','','TP3','A',1,'',W1[0][0],W1[0][1],W1[0][2],1.00,0.00,'O',''))
lines.append(line.format('ATOM',str(2),'H1','','TP3','A',1,'',W1[1][0],W1[1][1],W1[1][2],1.00,0.00,'H',''))
lines.append(line.format('ATOM',str(3),'H2','','TP3','A',1,'',W1[2][0],W1[2][1],W1[2][2],1.00,0.00,'H',''))
lines.append(line.format('ATOM',str(4),'OH2','','TP3','A',2,'',W2[0][0],W2[0][1],W2[0][2],1.00,0.00,'O',''))
lines.append(line.format('ATOM',str(5),'H1','','TP3','A',2,'',W2[1][0],W2[1][1],W2[1][2],1.00,0.00,'H',''))
lines.append(line.format('ATOM',str(6),'H2','','TP3','A',2,'',W2[2][0],W2[2][1],W2[2][2],1.00,0.00,'H',''))  

outFile = open('/home/noel/Projects/FixMD/Aguan/tests/validation/charmm/two_waters.pdb', 'w')
for i in lines:
    outFile.write(i+'\n')
outFile.close()

   