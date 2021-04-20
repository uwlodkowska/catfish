#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

import mahotas as mh
import cv2, copy
from PIL import Image



from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, opening, square, disk
from skimage.draw import disk as pdisk
from skimage.color import label2rgb, rgb2gray

from scipy import ndimage as ndi

from skimage.morphology import watershed
from skimage.feature import peak_local_max

from skimage.filters import threshold_local, threshold_otsu, threshold_multiotsu, rank
from skimage.filters import sobel

from scipy.stats import skew
from scipy.signal import argrelextrema
from common import Approximation


def divide_region_approx_ellipse(region, last_id):
    new_nuclei = np.array([])
    alpha = (-180/math.pi)*(region.orientation)
    shift_x = 0.25*region.major_axis_length*math.sin(math.radians(-alpha))
    shift_y = 0.25*region.major_axis_length*math.cos(math.radians(-alpha))
    approxEllipse = mpatches.Ellipse((region.centroid[1]-shift_x,region.centroid[0]-shift_y), region.minor_axis_length,
                                                             region.major_axis_length/2,alpha,
                           fill=False, edgecolor='red', linewidth=1)  
    new_nucleus = Approximation(region, copy.copy(approxEllipse), last_id)
    
    new_nuclei = np.append(new_nuclei, new_nucleus) 
                
    approxEllipse = mpatches.Ellipse((region.centroid[1]+shift_x,region.centroid[0]+shift_y), 
                                                 region.minor_axis_length,
                                                             region.major_axis_length/2,alpha,
                           fill=False, edgecolor='red', linewidth=1)  
    new_nucleus = Approximation(region, copy.copy(approxEllipse), last_id+1)
    new_nuclei = np.append(new_nuclei, new_nucleus)
    return new_nuclei


def prepare_region_approximation(region, last_id):
    alpha = (-180/math.pi)*(region.orientation)
    if region.minor_axis_length == 0:
        ax_ratio = 0
    else:
        ax_ratio = region.major_axis_length/region.minor_axis_length
    new_nuclei = np.array([])
    #suwak + histogram do ustalania tej wartości?
    if ax_ratio * region.area > 1500:
        new_nuclei = divide_region_approx_ellipse(region, last_id)
    else:
        approxEllipse = mpatches.Ellipse((region.centroid[1],region.centroid[0]), region.minor_axis_length,
                                                             region.major_axis_length,alpha,
                           fill=False, edgecolor='red', linewidth=1)  
        new_nuclei = np.append(new_nuclei, Approximation(region, copy.copy(approxEllipse), last_id))

    
    return new_nuclei


def dapi_segmentation(image_ts, footprint_size, thresh_down, thresh_up, thresh_range, sobel_g=False): 
    ratio_thre = 1.6
    
    image_ts_orig = copy.copy(image_ts)
    
    image_ts = opening(image_ts, disk(2))
    

    thresh = threshold_otsu(image_ts[image_ts>0])/2
    

    local_thresh = threshold_local(image_ts, thresh_range, offset=0)
    binary_local = image_ts > local_thresh

    
    local_masked_global = np.ma.masked_array(data = binary_local, 
                                             mask = ~(image_ts > thresh), fill_value = 0).filled()

    local_masked_global = opening(local_masked_global, disk(5))
    elev_map = sobel(mh.gaussian_filter(image_ts,3))
    distance = ndi.distance_transform_edt(local_masked_global)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((footprint_size,footprint_size)),
                                labels=local_masked_global)
    markers = ndi.label(local_maxi)[0]
   
    labels = watershed(-elev_map, markers, mask=local_masked_global)
    label_image = label(labels, connectivity = 1)

    areas = np.array([])
    intensities = np.array([])
    ellipse_axes_ratios = np.array([])
    
    for region in regionprops(label_image, intensity_image=image_ts):
        areas = np.append(areas, region.area)
        intensities = np.append(intensities, region.mean_intensity)
        if region.minor_axis_length == 0:
            ax_ratio = 0
        else:
            ax_ratio = region.major_axis_length/region.minor_axis_length
        ellipse_axes_ratios = np.append(ellipse_axes_ratios, ax_ratio)

    nuclei = np.array([])
    for region in regionprops(label_image):
        #tu threshold do premyślenia
        if region.area > 10:
            new_nuclei = prepare_region_approximation(region, len(nuclei))
        nuclei = np.append(nuclei, new_nuclei)   

    
    return nuclei, ellipse_axes_ratios, areas, intensities



