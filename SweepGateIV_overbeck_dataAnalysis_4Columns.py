# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:40:09 2017

@author: ovj
"""

import numpy as np
import re
import os
import codecs
import copy
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
import numpy as np
import matplotlib.pyplot as plt
import csv
import os


#==============================================================================
# Define Function for searching
#==============================================================================

""" Input Path here!!!"""

#dir = os.getcwd()

Wdir = r"C:\Users\ovj\Desktop"
dir = Wdir + '\input'


listFiles = os.listdir(dir)
#==============================================================================
# removing directories from listfiles
#==============================================================================
newfiles = []
for file in listFiles:
    if file.find(".") != -1:
        newfiles.append(file)
listFiles=newfiles
del newfiles


#==============================================================================
# make figure for overall v-gate plot for all files in input
#==============================================================================

fig3 = plt.figure(fname+"R_v_VgateYoko_all") # a new figure window
ax3 = fig3.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
ax3.set_title(fname)
ax3.set_xlabel("Vgate[V]")
ax3.set_ylabel("Resistance [Ohm]")

  
#==============================================================================
# Cycle over all input files                
#==============================================================================
for x in listFiles:
    print("Opening: " + x)
#    filepath=os.path.join(dir, file) 
    filepath=os.path.join(dir, x) 
    
    fname = filepath.split("\\")[-1] #Get the last part of the filepath as filename.
#    header = [["Wavelength","Intensity", "Raman shift", "Intensity - Bgrd.", "Bgrd.", "Lorentz Fit"],["nm","CCD counts","rel. cm-1","CCD counts","CCD counts","CCD counts"]]
#    outdata = [] # will get the above structure
    
#    with open(filepath, 'r') as csvfile:  # because I've used with, I don't need to close it afterwards...
#        dialect = csv.Sniffer().sniff(csvfile.read(1024))
#        delim = dialect.delimiter
    
    a = np.loadtxt(filepath ,skiprows=12)#, delimiter=delim ) # No need to skip rows without header 
    Vgate = a[:,0]
    Rvgate = a[:,1]
    rawdata = a[:,2:] # has structure |Vbias|Current|GateVoltageDAQ|GateVoltageYoko|Vbias|Current|GateVoltageDAQ|GateVoltageYoko|...
    rawVbias = rawdata[:,0::4]
    rawCurrent = rawdata[:,1::4]
    rawGateDAQ = rawdata[:,2::4]
    rawGateYoko = rawdata[:,3::4]
    listResistances = np.ndarray(len(rawVbias[0]))
    cycleIndex = np.arange(len(rawVbias[0]))
    

    for i in cycleIndex:
        
        x_bias = rawVbias[:,i]# applied bias in V
        y_curr = rawCurrent[:,i] # measured current in Amps
        
        # fitting the background to a line #
        m, c = np.polyfit(x_bias, y_curr, 1)
        Rfit = 1/m
        listResistances[i]= Rfit
                       
#==============================================================================
#        Plot R vs. cycle number                
#==============================================================================
       
    fig1 = plt.figure(fname+"R_v_CycleNo") # a new figure window
    ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
    ax1.set_title(fname)
    ax1.set_xlabel("Cycle number")
    ax1.set_ylabel("Resistance [Ohm]")
    ax1.plot(cycleIndex[1:-1],listResistances[1:-1])
    
    try:
        DateSampleName= fname.split("_J")[0]
    except:
        DateSampleName=""
    
    savepath = filepath.split(fname)[0]+DateSampleName+"_Plots\\" # create subdirectory
    if not os.path.exists(savepath):
        os.makedirs(savepath)
            
    fig1.savefig(savepath+fname+"_R_v_CycleNo"+".png")
    

                       
#==============================================================================
#        Plot R vs. Vgate Yoko                
#==============================================================================
       
    fig2 = plt.figure(fname+"R_v_VgateYoko") # a new figure window
    ax2 = fig2.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
    ax2.set_title(fname)
    ax2.set_xlabel("Vgate[V]")
    ax2.set_ylabel("Resistance [Ohm]")
    ax2.plot(rawGateYoko[0][1:-1],listResistances[1:-1])
    fig2.savefig(savepath+fname+"_R_v_VgateYoko"+".png")
    
    
#==============================================================================
#     Add R vs. Vgate Yoko to overall V-gate-Plot
#==============================================================================
    ax3.plot(rawGateYoko[0][1:-1],listResistances[1:-1])


#==============================================================================
# Save overall V-gate-Plot
#==============================================================================
fig3.savefig(savepath+fname+"_R_v_VgateYoko"+".png")
#        
#        
##==============================================================================
##         Plot all IV-curves
##==============================================================================
#                       
#              
#       
#        
#    for i in range(450,500):
#        fig1 = plt.figure(fname+"R_v_Bias") # a new figure window
#        ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#        ax1.set_title(fname)
#        ax1.set_xlabel("Gate Volatage [mV]")
#        ax1.set_ylabel("Current [nA]")
#        ax1.plot(rawVbias[:,i]*10**3,rawCurrent[:,i]*10**9,'x')
#
#
#    fig1 = plt.figure(fname+"R_v_Gate") # a new figure window
#    ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#    ax1.set_title(fname)
#    ax1.set_xlabel("Gate Volatage [V]")
#    ax1.set_ylabel("Resistance [kOhm]")
#    ax1.plot(Vgate,Rvgate)
#    #ax1.legend(['Centre: %d rel. cm-1'%best_parameters[1], 'Bgrd. substr.'])
#    #fig1        
#    
#    
#    fig2 = plt.figure(fname+"R_v_Bias") # a new figure window
#    ax1 = fig2.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#    ax1.set_title(fname)
#    ax1.set_xlabel("Bias Volatage [mV]")
#    ax1.set_ylabel("Current [nA]")
#    ax1.plot(rawdata[:,0]*10**3,rawdata[:,1]*10**9)
#    #ax1.legend(['Centre: %d rel. cm-1'%best_parameters[1], 'Bgrd. substr.'])
#    #fig2        
    
    