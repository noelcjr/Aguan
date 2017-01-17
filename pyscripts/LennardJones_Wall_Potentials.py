# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:27:42 2016

@author: noel Carrascal
"""
import pandas as pd
import matplotlib.pyplot as plt

ro = 0.315061
ep = 0.6
#= (2^(1/6))*TM.ro;
#wallCutOff = 1.122462*ro
regionX = 1.8662
x = [float(i)/100 for i in range(95,150,1)]
fcValx = []
uSum = []
uSum2 = []
for i in x:
    dr = i
    rri = 1/(dr*dr)
    rri3 = rri*rri*rri
    fcValx.append(48 * rri3 * (rri3 - 0.5) *rri *dr)
    uSum.append(4 * rri3 * (rri3 - 1))
    uSum2.append(4 * rri3 * (rri3 - 1) + 1)

DF = pd.DataFrame(columns=['d','E','E2','F'])
DF['d'] = pd.Series(x)
DF['E'] = pd.Series(uSum)
DF['E2'] = pd.Series(uSum2)
DF['F'] = pd.Series(fcValx)
plt.figure()
plt.plot(DF['d'],DF['E'])
plt.plot(DF['d'],DF['E2'],color='red')
plt.axvline(x=1.12)
plt.axhline(y=0.0)
plt.title("Lennard-Jones and Wall Potentials")
plt.xlabel('ro (ro=0.315061 nm)')
plt.ylabel('e (e=0.6364 kj/mol)')
plt.show()
DF.tail(55)

ms1 = m1 * TM.sitesMol;
for(int j1 = 0; j1 < TM.sitesMol ; j1++){
      TM.fxs[ms1+j1] = TM.fxs[ms1+j1] + fcValx
