# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:40:09 2017

@author: ovj


%reset does magic



"""
import os
os.chdir(r"C:\Users\ovj\Documents\PythonScripts\Data Analysis\python-transport-analysis")
import numpy as np
import re
import os
import codecs
import copy
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.colors as colors
import matplotlib.cm as cmaps
import scipy.constants as const
from scipy import interpolate
from SR830_TransportDataAnalysis_Functions import savitzky_golay

#==============================================================================
# Define Function for searching
#==============================================================================

""" Input Path here!!!"""

#dir = os.getcwd()

Wdir = r"C:\Users\ovj\Desktop"
dir = Wdir + '\input'


listFiles = os.listdir(dir)



#==============================================================================
# Define Sort-Parameters for 
#==============================================================================


sortBlue = "xx+32_"
sortRed = "xx-32_"


#==============================================================================
# Specify fit position for mobility extraction
# If a string is specified, the position of the CNP is searched (min(G)) a fit done at +- 10V of this value.
#==============================================================================

fitH, fitE = "search", "search"



#==============================================================================
# removing lockin & directories from listfiles
#==============================================================================
newfiles = []
for file in listFiles:
    if file.find("G(Vgs)") != -1 and file.find(".dat") != -1:  # naming Lujun
#    if file.find("G(Vgs)") != -1 and file.find(".dat") != -1:  # naming Jan
        newfiles.append(file)
listFiles=newfiles
del newfiles

filepath=os.path.join(dir, listFiles[0]) 
    
fname = filepath.split("\\")[-1]

#==============================================================================
# define colourmap
#==============================================================================

cstart = 0.0
cstop = 0.30
nLinesInPlot = len(listFiles)
cm_subsection = np.linspace(cstart, cstop, nLinesInPlot) 
del cstart, cstop, nLinesInPlot

colors = [ cmaps.terrain(x) for x in cm_subsection ]








##==============================================================================
## Colourmap Example
##==============================================================================
#import matplotlib.pyplot as plt
#
#from matplotlib import cm
#from numpy import linspace
#
#start = 0.0
#stop = 1.0
#number_of_lines= 1000
#cm_subsection = linspace(start, stop, number_of_lines) 
#
#colors = [ cmaps.jet(x) for x in cm_subsection ]
#
#for i, color in enumerate(colors):
#    plt.axhline(i, color=color)
#
#plt.ylabel('Line Number')
#plt.show()
#
#


#==============================================================================
# make figure for overall v-gate plot for all files in input
#==============================================================================
figname=fname.split(".dat")[0]


fig3 = plt.figure(figname+"_R_v_Gate_All", figsize = (12, 4)) # a new figure window (normal length of legend: 9,4)
ax3 = fig3.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#ax3.set_title(figname+"_all")
ax3.set_xlabel("Vgate [V]")
ax3.set_ylabel("Resistance [Ohm]")
fig3.subplots_adjust(left=0.1, right=0.4, top=0.95, bottom=0.12)

figa = plt.figure(figname+"_R_v_Gate_All_a", figsize = (12, 4)) # a new figure window
axa = figa.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#axa.set_title(figname+"_all_a")
axa.set_xlabel("Vgate [V]")
axa.set_ylabel("Resistance [Ohm]")
figa.subplots_adjust(left=0.1, right=0.4, top=0.95, bottom=0.12)

figb = plt.figure(figname+"_R_v_Gate_All_b", figsize = (12, 4)) # a new figure window
axb = figb.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#axb.set_title(figname+"_all_b")
axb.set_xlabel("Vgate [V]")
axb.set_ylabel("Resistance [Ohm]")
figb.subplots_adjust(left=0.1, right=0.4, top=0.95, bottom=0.12)

figsc = plt.figure(figname+"_R_v_Gate_All_scaled", figsize = (12, 4)) # a new figure window
axsc = figsc.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#axsc.set_title(figname+"_all_scaled")
axsc.set_xlabel("Vgate [V]")
axsc.set_ylabel("Resistance [Rmin(1)]")
figsc.subplots_adjust(left=0.1, right=0.4, top=0.95, bottom=0.12)

##==============================================================================
## old
##==============================================================================
#figsc = plt.figure(figname+"_R_v_Gate_All_scaled") # a new figure window
#axsc = figsc.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
##axsc.set_title(figname+"_all_scaled")
#axsc.set_xlabel("Vgate [V]")
#axsc.set_ylabel("Resistance [Rmin(1)]")



#==============================================================================
# Define output data-sets
#==============================================================================
listmob = np.zeros([len(listFiles),3])  # structure: [measN|mu_H|mu_E]
listVcnp = np.zeros([len(listFiles),2])  # structure: [measN|Vcnp]
listMaxR = np.zeros([len(listFiles),2])  # structure: [measN|maxR]
  
#==============================================================================
# Cycle over all input files                
#==============================================================================
file = listFiles[0]#debug             
cycleN = 0
                
for file in listFiles:
    print("Opening: " + file)
#    filepath=os.path.join(dir, file) 
    filepath=os.path.join(dir, file) 
    
    fname = filepath.split("\\")[-1]#Get the last part of the filepath as filename.
    figname=fname.split(".dat")[0]
#    header = [["Wavelength","Intensity", "Raman shift", "Intensity - Bgrd.", "Bgrd.", "Lorentz Fit"],["nm","CCD counts","rel. cm-1","CCD counts","CCD counts","CCD counts"]]
#    outdata = [] # will get the above structure
    
#    with open(filepath, 'r') as csvfile:  # because I've used with, I don't need to close it afterwards...
#        dialect = csv.Sniffer().sniff(csvfile.read(1024))
#        delim = dialect.delimiter
    
    a = np.loadtxt(filepath ,skiprows=1)#, delimiter=delim ) # No need to skip rows without header 
    Vgate = a[:,0]
    calcG = a[:,4] # conductance in units of G0
    calcR = 2*12960 / a[:,4] # Resistance in Ohm
    
    scandir = np.sign(Vgate[-1])                 
                       
    minG = min(calcG)
    maxR = max(calcR)
    #Vcnp=Vgate[np.where(calcG==calcG.min())[0][0]]
    Vcnp=round(Vgate[np.where(savitzky_golay(calcG,41,3)==min(savitzky_golay(calcG,41,3)))[0][0]],4)
    
    calcCurrent = a[:,3] # Units = A
    LockinX = a[:,1] # Units = V  - True if Output is set to X
    LockinYorPhase = a[:,2] # Units = V or degree - this depends on whether the output is set to Y or Display!
  
    Glength=3
    Gwidth=2
    area = Glength*Gwidth*1E-12 # L x W [m^2]
    cap = const.epsilon_0*area/(0.6*1E-6) # [F] put in sample-dependent capacitance data
    n = (Vgate-Vcnp)*cap/(area*1E4*const.e)
    


    if type(fitH)==str:
        #fitH1 = round(max(Vgate[np.where(calcG == minG)[0][0]]-10,-30),4)
        fitH1 = round(Vgate[np.where(scandir*savitzky_golay(calcG,201,3,1)==min(scandir*savitzky_golay(calcG,201,3,1)))[0][0]])
#        plt.plot(Vgate,savitzky_golay(calcG,201,3,1)*1000) # debug
        #fitE1 = round(min(Vgate[np.where(calcG == minG)[0][0]]+10,30),4)
        fitE1 = round(Vgate[np.where(scandir*savitzky_golay(calcG,201,3,1)==max(scandir*savitzky_golay(calcG,201,3,1)))[0][0]])
        
    
#==============================================================================
# Fitting                      
#==============================================================================
    startH, stopH = np.where(Vgate == fitH1)[0][0]-20, np.where(Vgate == fitH1)[0][0]+21               
    startE, stopE = np.where(Vgate == fitE1)[0][0]-20, np.where(Vgate == fitE1)[0][0]+21

    m_H, c_H = np.polyfit(n[startH:stopH], calcG[startH:stopH], 1)                  
    m_E, c_E = np.polyfit(n[startE:stopE], calcG[startE:stopE], 1)        

    mu_H = -1*(Glength/Gwidth)*(1/const.e)*m_H*((const.e)**2)/const.h
    mu_E = (Glength/Gwidth)*(1/const.e)*m_E*((const.e)**2)/const.h  

#==============================================================================
# Append calc. mob. to list
#==============================================================================
    listmob[cycleN] = [cycleN, mu_H, mu_E]
    listVcnp[cycleN] = [cycleN, Vcnp]
    listMaxR[cycleN] = [cycleN, maxR]
    
    

             
#    
#==============================================================================
#                       Create Save-Path
#==============================================================================
                          
    try:
        DateSampleName= fname.split("_")[1]+"_"+fname.split("_")[2]+"_"+fname.split("_")[3]
    except:
        DateSampleName=""
    
    savepath = filepath.split(fname)[0]
    
    
    

#==============================================================================
#        Plot G or R vs. Vg             
#==============================================================================
       
    fig1 = plt.figure(figname+"G_v_Vgate") # a new figure window
    ax1 = fig1.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
    #ax1.set_title(fname)
    
    ax1.set_xlabel("Vgate [V]")
#    ax1.set_ylabel("Conductance [e^2/h]")
#    ax1.plot(Vgate,calcG,'b.', label=figname)

    ax1.set_ylabel("Resistance [Ohm]")
    
    if scandir ==1.0:    
        ax1.plot(Vgate,calcR,'b>', markersize = 3, label=figname)  # indicate direction >
    else:
        ax1.plot(Vgate,calcR,'b<', markersize = 3, label=figname) # indicate direction <

    # Place a legend above this subplot, expanding itself to
    # fully use the given bounding box.
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
               ncol=1, mode="expand", borderaxespad=0.)
#    print("test0")


#==============================================================================
#   Plot G vs. n + Mobility
#==============================================================================

#
    figmob = plt.figure(figname+"G_v_n")
    axmob = figmob.add_subplot(111)
    
    axmob.set_xlabel("chargecarrier density [1/cm^2]")
    axmob.set_ylabel("Conductance [e^2/h]")
    
    if Vgate[0] <0:    
        markdir = ">"
        axmob.plot(n,calcG, marker= markdir, c='k', markersize = 2)  # indicate direction >
    else:
        markdir = "<"
        axmob.plot(n,calcG, marker= markdir, c='k', markersize = 2) # indicate direction <
    lhole = "Hole mobility = " + str(round(mu_H,1)) + " cm^2/Vs"
    lelec = "Electron mobility = " + str(round(mu_E,1)) + " cm^2/Vs"
    axmob.plot(n[startH:stopH],n[startH:stopH]*m_H+c_H,'ro', label=lhole)   
    axmob.plot(n[startE:stopE],n[startE:stopE]*m_E+c_E,'bo', label=lelec)   
    axmob.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
               ncol=1, mode="expand", borderaxespad=0.)

#==============================================================================
#        Plot R vs. Vgate & Phase
#==============================================================================
##       
#    fig2 = plt.figure(figname+"R_v_Vgate") # a new figure window
#    ax2 = fig2.add_subplot(111) # ax1 is an Axes element ("plotting Window"). Specify (nrows, ncols, axnum)
#    ax2.set_title(fname)
#    ax2.set_xlabel("Vgate [V]")
#    ax2.set_ylabel("Resistance [Ohm]")
#    ax2.plot(Vgate,calcR, 'r.')
#    fig2.savefig(savepath+figname+"_R_v_Vgate"+".png")
#    
##    
##    
#    fig2, ax2 = plt.subplots()
#    ax2.plot(Vgate,calcR, 'b-')
#    ax2.set_xlabel("Vgate [V]")
#    # Make the y-axis label, ticks and tick labels match the line color.
#    ax2.set_ylabel("Resistance [Ohm]", color='b')
#    ax2.tick_params("x", colors='b')
#    
#    ax4 = ax2.twinx()
#    ax4.plot(Vgate,LockinYorPhase, 'r.')
#    ax4.set_ylabel("Lock-In Phase [°]", color='r')
#    ax4.tick_params("y", colors='r')
#    
#    fig2.tight_layout()
#    plt.show()
#    
    
#==============================================================================
#     Add R vs. Vgate Yoko to overall V-gate-Plot
#==============================================================================
    if fname.find("a.dat") !=-1:
        savepatha = filepath.split(fname)[0]+DateSampleName+"_a_Plots\\" # create subdirectory
        if not os.path.exists(savepatha):
            os.makedirs(savepatha)
            
#==============================================================================
#       save each individual plot      
#==============================================================================
                
        fig1.savefig(savepatha+figname+"_G_v_Vgate"+".png")
        plt.close(fig1)
        
        figmob.savefig(savepatha+figname+"_G_v_n+mob"+".png")
        plt.close(figmob)
        
#==============================================================================
#       creat cumulative plots        
#==============================================================================


        if fname.find(sortBlue) != -1:
            ax3.plot(Vgate,calcR, marker = markdir, c='b', markersize = 3 , label = figname)
            axa.plot(Vgate,calcR, marker = markdir, c='b', markersize = 3 , label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, c='b', markersize = 3 , label = figname)            
        elif fname.find(sortRed) != -1:
            ax3.plot(Vgate,calcR, marker = markdir, c='r', markersize = 3 , label = figname)
            axa.plot(Vgate,calcR, marker = markdir, c='r', markersize = 3 , label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, c='r', markersize = 3 , label = figname)
        else:
            ax3.plot(Vgate,calcR, marker = markdir, markersize = 3 , color =  colors[cycleN], label = figname)
            axa.plot(Vgate,calcR, marker = markdir, markersize = 3 , color =  colors[cycleN], label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, markersize = 3 , color =  colors[cycleN], label = figname)
#            print("test1")
            
    elif fname.find("b.dat") !=-1:
        savepathb = filepath.split(fname)[0]+DateSampleName+"_b_Plots\\" # create subdirectory
        if not os.path.exists(savepathb):
            os.makedirs(savepathb)
            
#==============================================================================
#       save each individual plot      
#==============================================================================
                
        fig1.savefig(savepathb+figname+"_G_v_Vgate"+".png")
        plt.close(fig1)        

        figmob.savefig(savepathb+figname+"_G_v_n+mob"+".png")
        plt.close(figmob)


#==============================================================================
#       creat cumulative plots        
#==============================================================================
        
        if fname.find(sortBlue) != -1:
            ax3.plot(Vgate,calcR, marker = markdir, c='b', markersize = 3 ,  label = figname)
            axb.plot(Vgate,calcR, marker = markdir, c='b', markersize = 3 , label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, c='b', markersize = 3 , label = figname)      
        elif fname.find(sortRed) != -1:
            ax3.plot(Vgate,calcR, marker = markdir, c='r', markersize = 3 , label = figname)
            axb.plot(Vgate,calcR, marker = markdir, c='r', markersize = 3 , label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, c='r', markersize = 3 , label = figname)
        else:
            ax3.plot(Vgate,calcR, marker = markdir, markersize = 3 , color =  colors[cycleN], label = figname)
            axb.plot(Vgate,calcR, marker = markdir, markersize = 3 ,color =  colors[cycleN], label = figname)
            axsc.plot(Vgate,calcR/maxR, marker = markdir, markersize = 3 ,color =  colors[cycleN], label = figname)
#            plt.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
#               ncol=1, mode="expand", borderaxespad=0.)
#            print("test2")
    cycleN +=1                 
#    ax3.plot(Vgate,calcR)          




#==============================================================================
# Plot mobility vs. measurement N
#==============================================================================


figmob2 = plt.figure(DateSampleName+"_mobilities vs. measurement")
axmob2 = figmob2.add_subplot(111)
axmob2.set_xlabel("measurement No")
axmob2.set_ylabel("mobility [cm^2/Vs]")
axmob2.plot(listmob[:,0],listmob[:,1],'ro', label="hole mobility")
axmob2.plot(listmob[:,0],listmob[:,2],'bo', label="electron mobility")
axmob2.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
               ncol=1, mode="expand", borderaxespad=0.)

#==============================================================================
# Plot maxR vs. measurment N
#==============================================================================

figMaxR = plt.figure(DateSampleName+"_Max. R vs. measurement")
axMaxR = figMaxR.add_subplot(111)
axMaxR.set_xlabel("measurement No")
axMaxR.set_ylabel("Resistance [Ohm]")
axMaxR.plot(listVcnp[:,0],listVcnp[:,1],'ro', label="Max R")
axMaxR.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
               ncol=1, mode="expand", borderaxespad=0.)
#==============================================================================
# Plot Vcnp vs. measurement N
#==============================================================================

figVcnp = plt.figure(DateSampleName+"_Vcnp  vs. measurement")
axVcnp = figVcnp.add_subplot(111)
axVcnp.set_xlabel("measurement No")
axVcnp.set_ylabel("Vcnp [V]")
axVcnp.plot(listVcnp[:,0],listVcnp[:,1],'ro', label="Vcnp")
axVcnp.legend(bbox_to_anchor=(0., 1.02, 1., 0.102), loc=3,
               ncol=1, mode="expand", borderaxespad=0.)

#==============================================================================
# Place a legend above the multi-line-plots, expanding itself to
# fully use the given bounding box.
#==============================================================================
ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axa.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axb.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axsc.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


#==============================================================================
# Save multiple-lines plots
#==============================================================================

fig3.savefig(savepath+figname+".png")
figa.savefig(savepatha+figname+"_a.png")
figb.savefig(savepathb+figname+"_b.png")

#figsc.savefig(savepathb+figname+"_scaled.png")

figmob2.savefig(savepath+DateSampleName+"_mob_v_meas")
figMaxR.savefig(savepath+DateSampleName+"_mob_v_MaxR")
figVcnp.savefig(savepath+DateSampleName+"_mob_v_Vcnp")
plt.close(figsc) # if one doesn't need this every time
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
    
    