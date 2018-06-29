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
    rawdata = a[:,2:] # has structure |Vbias|Current|GateVoltage|Vbias|Current|GateVoltage|...
    rawVbias = rawdata[:,0::3]
    rawCurrent = rawdata[:,1::3]
    rawGateDAQ = rawdata[:,2::3]
    rawGateYoko = rawdata[:,3::3]
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
    

                       
#==============================================================================
#        Plot R vs. Vgate Yoko                
#==============================================================================
       
    fig1 = plt.figure(fname+"R_v_VgateYoko") # a new figure window
    ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
    ax1.set_title(fname)
    ax1.set_xlabel("Vgate[V]")
    ax1.set_ylabel("Resistance [Ohm]")
    ax1.plot(rawGateYoko[0],listResistances[1:])
#        
#        
##==============================================================================
##         Plot IV-curves
##==============================================================================
#                       
#              
    l1 = [450,470,500]   # select indices to plot
#   l1 = np.arange[len(rawVbias[0])]
        
    for i in l1:
        fig1 = plt.figure(fname+"R_v_Bias") # a new figure window
        ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
        ax1.set_title(fname)
        ax1.set_xlabel("Gate Volatage [mV]")
        ax1.set_ylabel("Current [nA]")
        ax1.plot(rawVbias[:,i]*10**3,rawCurrent[:,i]*10**9,'x')
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
    
    