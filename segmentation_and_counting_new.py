# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

import mahotas as mh

import copy

from skimage import io
from skimage.measure import label, regionprops
from skimage.morphology import closing, opening, square

from scipy import ndimage as ndi

from skimage.morphology import watershed
from skimage.feature import peak_local_max

from skimage.filters import threshold_local, threshold_otsu
from skimage.filters import sobel

path_orig = "/media/ula/D/mikroskop/CA1/MAX_CA1_glass_28_rat_1_hipp_left_skan_2019-08-05_flat.tif"

path_comm = "/media/ula/D/mikroskop/CA1/split_pages/MAX_CA1_glass_28_rat_1_hipp_left_skan_2019-08-05_page"

dapi_path = path_comm + "1.tif"
homer_path = path_comm + "2.tif"
arc_path = path_comm + "3.tif"





def get_frac_thresh(arr, frac):
    min_lum = np.min(arr)
    max_lum = np.max(arr)
    return min_lum + frac*(max_lum - min_lum)

'''
    class for objects representing identified nuclei
    
    approxEllipse - ellipse with the same parameters as those of regions
    identified by labeling function
'''
class Nucleus(): 
    def __init__(self, nucelusProps, approxEllipse,  idno):
        self.nucelusProps = nucelusProps
        self.approxEllipse = approxEllipse
        self.idno = idno
   

def ieg_segmentation(image_ts, footprint_size, thresh_down, thresh_up, 
                     thresh_range): 
    T_mean = image_ts.mean()
    opened = opening(image_ts > T_mean/4, square(3))
    bw = closing(opened, square(3))

    thre = get_frac_thresh(image_ts, 0.10)
    binary_local = image_ts > thre
    
    local_masked_global = np.ma.masked_array(data = binary_local, mask = ~bw, fill_value = 0).filled()

    elev_map = sobel(mh.gaussian_filter(image_ts,2))
    distance = ndi.distance_transform_edt(local_masked_global)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((footprint_size,footprint_size)),
                                labels=local_masked_global)
    markers = ndi.label(local_maxi)[0]

   
    labels = watershed(-elev_map, markers, mask=local_masked_global)

    label_image = label(labels, connectivity = 1)
    

    fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(label_image)
    ax[0].set_title("Labelled")
    
    ax[1].imshow(image_ts)
    ax[1].set_title("orig img")

    areas = np.array([])
    regions_to_reconsider = []

    for region in regionprops(label_image):
        areas = np.append(areas, region.area)

   
    nuclei = np.array([])

    for region in regionprops(label_image):
        minr, minc, maxr, maxc = region.bbox
        rect = mpatches.Circle((region.centroid[1], region.centroid[0]), region.major_axis_length, 
                               fill=False, edgecolor='red', linewidth=2)
        if region.area >= thresh_up:
            ax[0].add_patch(rect)
            regions_to_reconsider = np.append(regions_to_reconsider, region)
        elif region.area >= thresh_down:
            nuclei = np.append(nuclei, region)    
            ax[0].add_patch(rect)

    for a in ax:
        a.set_axis_off()

    fig.tight_layout()
    plt.show()
    return nuclei, regions_to_reconsider, label_image


'''
    Function for finding nuclei
    
    Parameters:
        image_ts - 16-bit tif image
        footprint_size - size of neighbourhood in peak detection
        thresh_down - min area [px] of labeled region that allows it to be 
                      identified as a nucelus
        thresh_up - max area [px] of labeled region that allows it to be 
                      identified as a nucelus
                      
        thresh_range - size of neighbourhood for local threshold calculation
        
    Returns:
      
'''
#TODO: usunąć footprint_size, thresholdy/wynieść do stałej???
def dapi_segmentation(image_ts, footprint_size, thresh_down, thresh_up, 
                      thresh_range): 
    
    thresh_global = threshold_otsu(image_ts[image_ts>0])/2
    

    local_thresh = threshold_local(image_ts, thresh_range, offset=0)
    binary_local = image_ts > local_thresh

    
    local_masked_global = np.ma.masked_array(data = binary_local, 
                                             mask = ~(image_ts > thresh_global), 
                                             fill_value = 0).filled()
    
    
    elev_map = sobel(mh.gaussian_filter(image_ts,4))
    distance = ndi.distance_transform_edt(local_masked_global)
    local_maxi = peak_local_max(distance, indices=False, 
                                footprint=np.ones((footprint_size,footprint_size)),
                                labels=local_masked_global)
    markers = ndi.label(local_maxi)[0]
   
    labels = watershed(-elev_map, markers, mask=local_masked_global)

    label_image = label(labels, connectivity = 1)
    

    fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(image_ts)
    ax[0].set_title("Labelled")
    
    ax[1].imshow(local_masked_global)
    ax[1].set_title("orig img")

    areas = np.array([])
    regions_to_reconsider = []

    for region in regionprops(label_image):
        areas = np.append(areas, region.area)


    nuclei = np.array([])
    
    for region in regionprops(label_image):
        minr, minc, maxr, maxc = region.bbox
        ellipse_angle = 90-(180/math.pi)*region.orientation
        approxEllipse = mpatches.Ellipse(region.centroid[::-1], 
                                         region.minor_axis_length, 
                                         region.major_axis_length,
                                         ellipse_angle,
                       fill=False, edgecolor='red', linewidth=1)
        if region.area >= thresh_up:
            ax[0].add_patch(approxEllipse)
            regions_to_reconsider = np.append(regions_to_reconsider, region)
        elif region.area >= thresh_down:
            new_nucleus = Nucleus(region, copy.copy(approxEllipse), len(nuclei))
            nuclei = np.append(nuclei, new_nucleus)      
            ax[0].add_patch(approxEllipse)
    
    for a in ax:
        a.set_axis_off()

    fig.tight_layout()
    plt.show()
    return nuclei, regions_to_reconsider, label_image



image = io.imread(dapi_path)




nuclei, regions_to_reconsider, thresholded = dapi_segmentation(image, 24, 99, 
                                                               1250, 39)
print(len(nuclei), len(regions_to_reconsider))
arc_im = io.imread(arc_path)
homer_im = io.imread(homer_path)

arc_locations, _,_ = ieg_segmentation(arc_im, 2, 7, 30, 21)
homer_locations, _,_ = ieg_segmentation(homer_im, 2, 7, 30, 21)



def find_ieg_colloc(ieg_locs, nuclei, to_reevaluate):
    ieg_positive = np.array([])
    fig, ax = plt.subplots(figsize=(6, 6))
    full_image = io.imread(path_orig)
    ax.imshow(full_image)
    for nucleus in nuclei:
        nucleus_patch = copy.copy(nucleus.approxEllipse)
        ieg_no = 0
        for ieg_dot in ieg_locs:      
            #if the centroid of area identified as ieg presence is inside of nucleus, add it
            #to the number of ieg collocalized with the nucleus
            
            tmp = nucleus.approxEllipse.contains_point(ieg_dot.centroid[::-1])
            
            if (tmp):
                ieg_no +=1
                if ieg_no == 1:
                    ieg_positive = np.append(ieg_positive, nucleus.idno)
                     
                    minr, minc, maxr, maxc = ieg_dot.bbox
                    
                    nucleus_patch.set_color('white')
        if ieg_no > 1:
            to_reevaluate = np.append(to_reevaluate, nucleus)
        ax.add_patch(nucleus_patch)   
    ax.set_axis_off()

    fig.tight_layout()
    plt.title("Nuclei with ")
    plt.savefig()
    return ieg_positive, to_reevaluate




to_reevalueate = []

arc_pos, arc_mult =  find_ieg_colloc(arc_locations, nuclei, to_reevalueate)
homer_pos, homer_mult =  find_ieg_colloc(homer_locations, nuclei, to_reevalueate)

print('arc len ', len(arc_pos), len(arc_mult))
print('homer len ', len(homer_pos),len(homer_mult))




len(np.intersect1d(arc_pos, homer_pos))









