#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:57:38 2021

@author: Mathew
"""


from skimage.io import imread
import os
import pandas as pd
from picasso import render
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from sklearn.cluster import DBSCAN
from skimage import filters,measure
from skimage.filters import threshold_local


# Change below:

path="/Users/Mathew/Documents/Current analysis/Mathew_230824/GluN1/"
pre_string="C_PSD95-JF552-JFX650_GluN1_CA1sr004_"

antibody_string='AF488.tif'
donor_string='JF552.tif'
acceptor_string='JFX650.tif'
FRET_string='FRET.tif'




# Folders to analyse:
    
def load_image(toload):
    
    image=imread(toload)
    
    return image

def z_project(image):
    
    mean_int=np.mean(image,axis=0)
  
    return mean_int

# Subtract background:
def subtract_bg(image):
    background = threshold_local(image, 11, offset=np.percentile(image, 1), method='median')
    bg_corrected =image - background
    return bg_corrected

def threshold_image_otsu(input_image):
    threshold_value=filters.threshold_otsu(input_image)  
    
    # threshold_value=input_image.mean()+3*input_image.std()
    print(threshold_value)
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

def threshold_image_standard(input_image,thresh):
     
    binary_image=input_image>thresh

    return binary_image

# Threshold image using otsu method and output the filtered image along with the threshold value applied:
    
def threshold_image_fixed(input_image,threshold_number):
    threshold_value=threshold_number   
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Label and count the features in the thresholded image:
def label_image(input_image):
    labelled_image=measure.label(input_image)
    number_of_features=labelled_image.max()
 
    return number_of_features,labelled_image
    
# Function to show the particular image:
def show(input_image,color=''):
    if(color=='Red'):
        plt.imshow(input_image,cmap="Reds")
        plt.show()
    elif(color=='Blue'):
        plt.imshow(input_image,cmap="Blues")
        plt.show()
    elif(color=='Green'):
        plt.imshow(input_image,cmap="Greens")
        plt.show()
    else:
        plt.imshow(input_image)
        plt.show() 
    
        
# Take a labelled image and the original image and measure intensities, sizes etc.
def analyse_labelled_image(labelled_image,original_image):
    measure_image=measure.regionprops_table(labelled_image,intensity_image=original_image,properties=('area','perimeter','centroid','orientation','major_axis_length','minor_axis_length','mean_intensity','max_intensity'))
    measure_dataframe=pd.DataFrame.from_dict(measure_image)
    return measure_dataframe

# This is to look at coincidence purely in terms of pixels

def coincidence_analysis_pixels(binary_image1,binary_image2):
    pixel_overlap_image=binary_image1&binary_image2         
    pixel_overlap_count=pixel_overlap_image.sum()
    pixel_fraction=pixel_overlap_image.sum()/binary_image1.sum()
    
    return pixel_overlap_image,pixel_overlap_count,pixel_fraction

# Look at coincidence in terms of features. Needs binary image input 

def feature_coincidence(binary_image1,binary_image2):
    number_of_features,labelled_image1=label_image(binary_image1)          # Labelled image is required for this analysis
    coincident_image=binary_image1 & binary_image2        # Find pixel overlap between the two images
    coincident_labels=labelled_image1*coincident_image   # This gives a coincident image with the pixels being equal to label
    coinc_list, coinc_pixels = np.unique(coincident_labels, return_counts=True)     # This counts number of unique occureences in the image
    # Now for some statistics
    total_labels=labelled_image1.max()
    total_labels_coinc=len(coinc_list)
    fraction_coinc=total_labels_coinc/total_labels
    
    # Now look at the fraction of overlap in each feature
    # First of all, count the number of unique occurances in original image
    label_list, label_pixels = np.unique(labelled_image1, return_counts=True)
    fract_pixels_overlap=[]
    for i in range(len(coinc_list)):
        overlap_pixels=coinc_pixels[i]
        label=coinc_list[i]
        total_pixels=label_pixels[label]
        fract=1.0*overlap_pixels/total_pixels
        fract_pixels_overlap.append(fract)
    
    
    # Generate the images
    coinc_list[0]=1000000   # First value is zero- don't want to count these. 
    coincident_features_image=np.isin(labelled_image1,coinc_list)   # Generates binary image only from labels in coinc list
    coinc_list[0]=0
    non_coincident_features_image=~np.isin(labelled_image1,coinc_list)  # Generates image only from numbers not in coinc list.
    
    return coinc_list,coinc_pixels,fraction_coinc,coincident_features_image,non_coincident_features_image,fract_pixels_overlap


        
# :pad the images:
    
antibody=load_image(path+pre_string+antibody_string)
donor=load_image(path+pre_string+donor_string)
acceptor=load_image(path+pre_string+acceptor_string)
FRET=load_image(path+pre_string+FRET_string)

donor_flat=donor
acceptor_flat=acceptor
FRET_flat=FRET

# Threshold antibody:
    
ab_thresh,antibody_binary=threshold_image_otsu(antibody)

# Threshold acceptor channel:
    
ac_thresh,acceptor_binary=threshold_image_otsu(acceptor_flat)


# both_binary=donor_binary*acceptor_binary

# Make FRET_image:
fret_im=acceptor_flat/(acceptor_flat+donor_flat)

# Threshold for only channels where acceptor was present:
fret_im_thresh=acceptor_binary*fret_im


# Save the images
imsr2 = Image.fromarray(fret_im_thresh)
imsr2.save(path+pre_string+'_FRET.tif')

imsr2 = Image.fromarray(acceptor_binary)
imsr2.save(path+pre_string+'_Acceptor_Binary.tif')

imsr2 = Image.fromarray(antibody_binary)
imsr2.save(path+pre_string+'_AB_Binary.tif')

number,labelled=label_image(acceptor_binary)
print("%d feautres were detected in the image."%number)

imsr2 = Image.fromarray(labelled)
imsr2.save(path+pre_string+'_label.tif')

measurements=analyse_labelled_image(labelled,fret_im)


       
frets=measurements['mean_intensity']
plt.hist(frets, bins = 20,range=[0,1], rwidth=0.9,color='#ff0000')
plt.xlabel('FRET Efficiency',size=20)
plt.ylabel('Number of Features',size=20)
plt.savefig(path+pre_string+"_FRET_hist.pdf")
plt.show()

areas=measurements['area']
plt.hist(areas, bins = 30,range=[0,30], rwidth=0.9,color='#ff0000')
plt.xlabel('Area (pixels)',size=20)
plt.ylabel('Number of Features',size=20)
plt.savefig(path+pre_string+"_area_hist.pdf")
plt.show()



length=measurements['major_axis_length']
plt.hist(length, bins = 5,range=[0,10], rwidth=0.9,color='#ff0000')
plt.xlabel('Length',size=20)
plt.ylabel('Number of Features',size=20)
plt.title('Cluster lengths',size=20)
plt.savefig(path+pre_string+"Lengths.pdf")
plt.show()

measurements.to_csv(path + '/' + pre_string+'_Metrics.csv', sep = '\t')

measurements2=analyse_labelled_image(labelled,donor_flat)
measurements2.to_csv(path + '/' + pre_string+'_Donor_Metrics.csv', sep = '\t')

# Now need to look at coincident only:

coinc_binary=antibody_binary*acceptor_binary

number,labelled=label_image(coinc_binary)
print("%d feautres were detected in the image."%number)
measurements=analyse_labelled_image(labelled,fret_im)
measurements.to_csv(path + '/' + pre_string+'_Coinc_Metrics.csv', sep = '\t')
measurements2=analyse_labelled_image(labelled,donor_flat)
measurements2.to_csv(path + '/' + pre_string+'_Donor_Coinc_Metrics.csv', sep = '\t')
