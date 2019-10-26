# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
#import glob
import matplotlib.pyplot as plt

dirpath = "/home/noel/Projects/FixMD/Phase_Diag/"
systems = ["t3_216_wt0_108_150EE_T2/","t3_2016_wt0_108_0EE_T2"]
temperature = ["298/","288/"]
density_runs = ["20X_20Y_30Z/","20X_20Y_29Z/","20X_20Y_28Z/","20X_20Y_27Z/",\
                "20X_20Y_26Z/","20X_20Y_25Z/","20X_20Y_24Z/","20X_20Y_23Z/",\
                "20X_20Y_22Z/","20X_20Y_21Z/"]
simulations = ['A','B','C','D','E','F']
cols = ['t','ENR','PE','VDW','EE','RF1','RF2','EF','KE','RKE','TKE','T',"TT","RT","P","Sim"]
col = {'P':'Pres','PE':'Pot. Eng.','KE':'Kit. Eng.','T':'Temp','VDW':'VDW Eng','EE':'EE Eng'}
ycolU = {'P':'KPa','PE':'Kj/mol','KE':'Kj/mol','T':'K','VDW':'Kj/mol','EE':'Kj/mol'}
###############################################################################
''' Inputs '''
syst = 0
temp = 0
simu = 0
d_run = 1
file_range = range(1,5+1)
value = 'P'
###############################################################################
###############################################################################
### on a given density, I want to compare results with simulations with alternative
### starting points. WARNINING THis will only work for six plots arrange as (2x3)
for d_r in range(len(density_runs)):
    print('-------------  Runing '+density_runs[d_r]+'-------------------')
    all_data = []
    for sim in range(len(simulations)):
        file_2_plot = dirpath+systems[syst]+temperature[temp]+density_runs[d_r]+simulations[sim]
        for i in file_range:
            fl = file_2_plot+"_"+str(i)+".out"
            #df = pd.DataFrame(columns=cols)
            df = pd.read_csv(fl,delimiter=' ', header=-1)
            df[15] = [simulations[sim] for i in range(min(df.index),max(df.index)+1)]
            all_data.append(df)
    
    DFout = pd.DataFrame()
    DFout = pd.concat(all_data)
    DFout.columns = cols
    DFout.index = DFout.t
    del DFout['t']
    DFout = DFout.drop(df.index[[0.0]])
    DFout.index = DFout.index.astype(int)
    
    df = pd.DataFrame()
    for i in simulations:
        df[i] = DFout[value][DFout.Sim == i].describe()
    
    mean_mean = np.mean(df.ix['mean'])
    
    f, axarr = plt.subplots(2, 3)
    axarr[0, 0].plot(DFout[value][DFout.Sim == simulations[0]])
    axarr[0, 0].axhline(y=df[simulations[0]]['mean'], color='red')
    axarr[0, 0].axhline(mean_mean, color='blue')
    axarr[0, 0].set_title(col[value]+" "+ycolU[value]+" "+simulations[0])
    
    
    axarr[0, 1].plot(DFout[value][DFout.Sim == simulations[1]])
    axarr[0, 1].axhline(mean_mean, color='blue')
    axarr[0, 1].axhline(y=df[simulations[1]]['mean'], color='red')
    axarr[0, 1].set_title(col[value]+" "+ycolU[value]+" "+simulations[1])
    axarr[0, 2].plot(DFout[value][DFout.Sim == simulations[2]])
    axarr[0, 2].axhline(mean_mean, color='blue')
    axarr[0, 2].axhline(y=df[simulations[2]]['mean'], color='red')
    axarr[0, 2].set_title(col[value]+" "+ycolU[value]+" "+simulations[2])
    axarr[1, 0].plot(DFout[value][DFout.Sim == simulations[3]])
    axarr[1, 0].axhline(mean_mean, color='blue')
    axarr[1, 0].axhline(y=df[simulations[3]]['mean'], color='red')
    axarr[1, 0].set_title(col[value]+" "+ycolU[value]+" "+simulations[3])
    axarr[1, 1].plot(DFout[value][DFout.Sim == simulations[4]])
    axarr[1, 1].axhline(mean_mean, color='blue')
    axarr[1, 1].axhline(y=df[simulations[4]]['mean'], color='red')
    axarr[1, 1].set_title(col[value]+" "+ycolU[value]+" "+simulations[4])
    axarr[1, 2].plot(DFout[value][DFout.Sim == simulations[5]])
    axarr[1, 2].axhline(mean_mean, color='blue')
    axarr[1, 2].axhline(y=df[simulations[5]]['mean'], color='red')
    axarr[1, 2].set_title(col[value]+" "+ycolU[value]+" "+simulations[5])
    # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
    plt.show()
    print('The red horizontal line is the mean for the simulation.')
    print('In blue horizontall line is the mean for all simulations = '+str(mean_mean))
    print df
    print('\n')
###############################################################################
###############################################################################
## This compares properties across different density on files that 
## have identical initial configurations
for den_run in density_runs:
    file_2_plot = dirpath+systems[syst]+temperature[temp]+den_run+simulations[simu]
    #DF = pd.read_csv(file_2_plot,delimiter=' ')
    
    #allFiles = glob.glob(dirpath + temperature + density_runs[run] + "/*A*.out")
    DFout = pd.DataFrame()
    all_data = []
    for i in file_range:
        fl = file_2_plot+"_"+str(i)+".out"
        #df = pd.DataFrame(columns=cols)
        df = pd.read_csv(fl,delimiter=' ', header=-1)
        all_data.append(df)
    DFout = pd.concat(all_data)
    DFout.columns = cols
    DFout.index = DFout.t
    del DFout['t']
    del DFout['NA']
    DFout = DFout.drop(df.index[[0.0]])
    #print DFout[value].describe()
    
    plt.figure()
    plt.plot(DFout[value])
    #plt.plot(DF['d'],DF['E2'],color='red')
    #plt.axvline(x=1.12)
    #plt.axhline(y=0.0)
    plt.title(col[value]+" "+den_run)
    plt.xlabel('t(ps)')
    plt.ylabel(ycolU[value])
    plt.show()

