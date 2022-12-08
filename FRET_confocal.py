#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 18:30:37 2021

@author: Mathew
"""


import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import os
from matplotlib.colors import LogNorm

# This is the array of folders we want to analyse
pathlist=[]

# Where to store the overall file containing means etc. for each experiment.
path_root=r"/Users/Mathew/Documents/Current analysis/Hannah_Factor_H/30_11_22/"

# # # Foldert to analyse here:
pathlist.append(r"/Users/Mathew/Documents/Current analysis/Hannah_Factor_H/30_11_22/Prep_1_RFP_elution/")


file_stem="FH"          # This is the part of the filename that will be searched for in each folder.

# Thresholds and other parameters:
# For AND thresholding" 
channelA_thresh=10     # Threshold for Channel A (Green).
channelB_thresh=10     # Threshold for Channel B (Red).

# For sum thresholding:
threshold=50

# Type of thresholding: 'AND' or 'SUM'
thresh_type='SUM'

# Autofl and cross-talk
channelA_AF=0       # Autofluorescence
channelB_AF=0
xtalk=0             # Cross-talk from A to B


#  This is the code to load the files
def load_files(filename_contains,path):
    print(path)
    num=0
    channelA_sample=[]             # Where channel A data will be stored
    channelB_sample=[]             # Where channel B data will be stored
    for root, dirs, files in os.walk(path):
      for name in files:
     
              if filename_contains in name:
                  if 'pdf' not in name:
                   

                      num+=1
                      a=0
                      with open(path+name) as csvDataFile:                                                # Opens the file as a CSV
                            csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
                            for row in csvReader:
                                channelA_sample.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
                                channelB_sample.append(row[1])
                                a+=1

    # print("Loaded %s files in total, with a total of %s rows"%(num,rows))
        
    channelA_arr_sample=np.asarray(channelA_sample,dtype=np.float32)                              # Converts these to numpy arrays for vector calcs.
    channelB_arr_sample=np.asarray(channelB_sample,dtype=np.float32)
    return channelA_arr_sample,channelB_arr_sample,num


for path in pathlist:

    # This is to load all of the files
    channelA_arr,channelB_arr,num=load_files(file_stem,path)
    
    # Now need to account for autofluorescence and crosstalk etc. 
    
    channelB_arr=(channelB_arr-xtalk*channelA_arr)-channelB_AF
    channelA_arr=channelA_arr-channelA_AF
    
    
    # Plot the data
    
    channelB_arr_inv=channelB_arr*(-1)
    
    plt.plot(channelA_arr,color='green')
    plt.plot(channelB_arr_inv,color='red')
    plt.xlabel('Bin number')
    plt.ylabel('Intensity (photons)')
    plt.xlim(0,10000)
    plt.ylim(-100,100)
    plt.savefig(path+'/'+'example_trace.pdf')
    plt.show()
    
    
    #This part is for the thresholding (AND):
    if thresh_type=='AND':
        channelA_only_events=channelA_arr[(channelA_arr>channelA_thresh)]                       # Total A events
        channelB_only_events=channelB_arr[(channelB_arr>channelB_thresh)]                       # Total B events
        channelA_events=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # A coincident events             
        channelB_events=channelB_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # B coincident events
    elif thresh_type=='SUM':
        channelA_only_events=channelA_arr[(channelA_arr>channelA_thresh)]                       # Total A events
        channelB_only_events=channelB_arr[(channelB_arr>channelB_thresh)]                       # Total B events
        channelA_events=channelA_arr[((channelA_arr+channelB_arr)>threshold)]  # A coincident events             
        channelB_events=channelB_arr[((channelA_arr+channelB_arr)>threshold)]
    # Calculate FRET efficiency from data
        
    
    FRET_events=channelB_events/(channelB_events+channelA_events)
   
    
    # Fraction of events
    if thresh_type=='AND':
        frac=len(channelB_events)/(len(channelB_events)+len(channelA_only_events)+len(channelB_only_events))
        print(frac)
    # Now want histogram:
    
    FRET_hist=np.histogram(FRET_events, bins=20, range=(0,1), normed=None, weights=None, density=None)
  
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = "12"
    plt.figure(figsize=(8, 6))
    plt.hist(FRET_events, bins = 20,range=[0,1], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,label="Real Events")  
    plt.xlabel('Proximity Ratio')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET.pdf')
    plt.show()
    
