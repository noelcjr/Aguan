# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 00:23:53 2016

@author: noel
"""
from Bio.PDB import *
from Bio.PDB.StructureBuilder import *
import copy
import numpy as np

def _calculate_center_of_mass(structure):
    total_mass = 0
    mx_total = 0
    my_total = 0
    mz_total = 0
    for atom in structure.get_atoms():
        coords = atom.coord.tolist()
        mass = atom.mass
        total_mass += mass
        mx_total += coords[0] * mass
        my_total += coords[1] * mass
        mz_total += coords[2] * mass
    return [mx_total/total_mass, my_total/total_mass, mz_total/total_mass]

def translate_molecule(structure, direction):
    for atoms in structure.get_atoms():
        atoms.set_coord(atoms.get_coord() + direction)

def center_molecule(center):
    return [x * -1 for x in center]

def genMatrix(TMrMatT, eAngx, eAngy, eAngz):
    '''p = new double[10]; tq = new double[4];'''
    p = np.zeros((1,10))
    tq = np.zeros((1,4))
    '''double a1, a2, a3, f;'''
    '''double eAngx, eAngy, eAngz;'''
    '''TM.rx[0] = 0.0;   TM.ry[0] = 0.0;    TM.rz[0] = 0.0;'''
    '''TM.rvx[0] = 0.0;  TM.rvy[0] = 0.0;   TM.rvz[0] = 0.0;'''
    '''TM.rax[0] = 0.0;  TM.ray[0] = 0.0;   TM.raz[0] = 0.0;'''
    '''eAngx = 0.0;      eAngy = 0.0;       eAngz = 0.0;'''
    '''a1 = 0.5 * eAngy;
       a2 = 0.5 * (eAngx - eAngz);
       a3 = 0.5 * (eAngx + eAngz);
    '''
    a1 = 0.5 * eAngy
    a2 = 0.5 * (eAngx - eAngz)
    a3 = 0.5 * (eAngx + eAngz)
    '''TM.q_u1[0] = Math.sin(a1) * Math.cos(a2);'''
    TMq_u1 = np.sin(a1) * np.cos(a2)
    '''TM.q_u2[0] = Math.sin(a1) * Math.sin(a2);'''
    TMq_u2 = np.sin(a1) * np.sin(a2)
    '''TM.q_u3[0] = Math.cos(a1) * Math.sin(a3);'''
    TMq_u3 = np.cos(a1) * np.sin(a3)
    '''TM.q_u4[0] = Math.cos(a1) * Math.cos(a3);'''
    TMq_u4 = np.cos(a1) * np.cos(a3)
    '''tq[0] = TM.q_u1[0];'''
    tq[0][0] = TMq_u1
    '''tq[1] = TM.q_u2[0];'''
    tq[0][1] = TMq_u2
    '''tq[2] = TM.q_u3[0];'''
    tq[0][2] = TMq_u3
    '''tq[3] = TM.q_u4[0];'''
    tq[0][3] = TMq_u4
    
    '''
    for(k = 0, k2 = 0; k2 < 4; k2++){  
           for(k1 = k2; k1 < 4; k1++, k++){
           k  ,  k1,  k2
           0     0     0
           1     1     0
           2     2     0
           3     3     0
           -------------
           4     1     1
           5     2     1
           6     3     1
           -------------
           7     2     2
           8     3     2
           -------------
           9     3     3
               p[k] = 2.0*tq[k1]*tq[k2];
    }}'''
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

def MulMat(a, b, c, nn, n):                    #MulMat(double[] a, double[] b, double[] c, int nn, int n){
    for i in range(0,nn):                     #   for(int i = 0; i < nn; i++){
        for j in range(0,nn):                 #      for(int j = 0; j < nn; j++){
            a[0][i + nn * j] = 0                 #          UNCHANGED
            for k in range(0,nn):             #          for(int k = 0; k < nn; k++){
                a[0][i + (nn * j)] = a[0][i + (nn * j)] + b[0][i + (nn * k)] * c[0][(n*9)+(k + (nn * j))];
                                              #}}}

def BuildStepRmat(mp, ax, ay, az):
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
            ak = az       #        ak = az}
                           
        t = 0.25*ak*ak     #    t = 0.25 * ak*ak
        c[0][k] = (1-t)/(1+t) #    c[k] = (1.0-t)/(1.0+t)
        s[0][k] = ak/(1+t)    #    s[k] = ak/(1.0+t)
    # The rest of this function is unchanged.
    c0c2 = c[0][0] * c[0][2]
    c0s2 = c[0][0] * s[0][2]
    s0c2 = s[0][0] * c[0][2]
    s0s2 = s[0][0] * s[0][2]

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

parser = PDBParser()

s = parser.get_structure('Insulin', '/home/noel/Projects/Protein_design/Insulin/OIPD/2hiu_1H.pdb')

center_of_mass = _calculate_center_of_mass(s)
translate_molecule(s,center_molecule(center_of_mass))
model_number = 1
coords = [5,0,0]
translate_molecule(s,coords)
TMrMatT = np.zeros((1,9))
genMatrix(TMrMatT, 0.0, 0.0, 0.0)
tx = TMrMatT[0][0]*coords[0] + TMrMatT[0][3]*coords[1] + TMrMatT[0][6]*coords[2]
ty = TMrMatT[0][1]*coords[0] + TMrMatT[0][4]*coords[1] + TMrMatT[0][7]*coords[2]
tz = TMrMatT[0][2]*coords[0] + TMrMatT[0][5]*coords[1] + TMrMatT[0][8]*coords[2]
print(tx,ty,tz)
#mc = new double[9];
coords = [1,0,0]
genMatrix(TMrMatT, 0.0, 0.0, 0.0)
mc = np.zeros((1,9))
#mt = new double[9];
mt = np.zeros((1,9))
print(0,coords)
for i in range(1,8):
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
    print(i,coords)
    #a = a + 0.1
      
''' This would be for rotation of the molecule around is center of mass
    s2 = copy.deepcopy(s)
    for atoms in s2.get_atoms():
        coords = atom.coord.tolist()
        tx = TMrMatT[9+0]*coords[0] + TMrMatT[9+3]*coords[1] + TMrMatT[9+6]*coords[2];
        ty = TMrMatT[9+1]*coords[0] + TMrMatT[9+4]*coords[1] + TMrMatT[9+7]*coords[2];
        tz = TMrMatT[9+2]*coords[0] + TMrMatT[9+5]*coords[1] + TMrMatT[9+8]*coords[2];

        atoms.set_coord([atoms.get_coord()[0] + tx, atoms.get_coord()[1] + ty, atoms.get_coord()[2] + tz])

    s2[0].id = model_number
    s.add(s2[0])
'''
        
io = PDBIO()
io.set_structure(s)
io.save('/home/noel/Projects/Protein_design/Insulin/out.pdb')
