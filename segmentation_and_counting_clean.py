
# coding: utf-8


#import matplotlib.pyplot as plt

#import matplotlib.patches as mpatches
import numpy as np

import mahotas as mh

from skimage import io
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, opening, square, disk
from skimage.color import label2rgb, rgb2gray

from scipy import ndimage as ndi

from skimage.morphology import watershed
from skimage.feature import peak_local_max

from skimage.filters import threshold_local, threshold_adaptive
from skimage.filters import sobel


#change to usr input
path_comm = "/media/ula/D/mikroskop/CA1/split_pages/MAX_CA1_glass_28_rat_1_hipp_left_skan_2019-08-05_page"


dapi_path = path_comm + "1.tif"
homer_path = path_comm + "2.tif"
arc_path = path_comm + "3.tif"



def get_frac_thresh(arr, frac):
    min_lum = np.min(arr)
    max_lum = np.max(arr)
    return min_lum + frac*(max_lum - min_lum)


class Nucleus(): 
    def __init__(self, nucelusProps, idno):
        self.nucelusProps = nucelusProps
        self.idno = idno
   

def ieg_segmentation(image_ts, footprint_size, thresh_down, thresh_up, thresh_range, sobel_g=False): 
    T_mean = image.mean()
    opened = opening(image_ts > T_mean, square(3))
    bw = closing(opened, square(3))

    thre = get_frac_thresh(image, 0.15)
    binary_local = image_ts > thre
    
    local_masked_global = np.ma.masked_array(data = binary_local, mask = ~bw, fill_value = 0).filled()

    elev_map = sobel(mh.gaussian_filter(image_ts,2))
    distance = ndi.distance_transform_edt(local_masked_global)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((footprint_size,footprint_size)),
                                labels=local_masked_global)
    markers = ndi.label(local_maxi)[0]

   
    labels = watershed(-elev_map, markers, mask=local_masked_global)

    label_image = label(labels, connectivity = 1)
    regions_to_reconsider = []
    nuclei = np.array([])
    
    for region in regionprops(label_image):
        minr, minc, maxr, maxc = region.bbox
        #rect = mpatches.Circle((region.centroid[1], region.centroid[0]), region.major_axis_length, 
        #                       fill=False, edgecolor='red', linewidth=2)
        if region.area >= thresh_up:
        #    ax[0].add_patch(rect)
            regions_to_reconsider = np.append(regions_to_reconsider, region)
        elif region.area >= thresh_down:
            nuclei = np.append(nuclei, region)

    return nuclei, regions_to_reconsider, label_image



def dapi_segmentation(image_ts, footprint_size, thresh_down, thresh_up, thresh_range, sobel_g=False): 
    T_mean = image.mean()
    opened = opening(image_ts > T_mean, square(5))
    bw = closing(opened, square(5))

    local_thresh = threshold_local(image_ts, thresh_range, offset=0)
    binary_local = image_ts > local_thresh

    
    local_masked_global = np.ma.masked_array(data = binary_local, mask = ~bw, fill_value = 0).filled()

    elev_map = sobel(mh.gaussian_filter(image_ts,1))
    distance = ndi.distance_transform_edt(local_masked_global)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((footprint_size,footprint_size)),
                                labels=local_masked_global)
    markers = ndi.label(local_maxi)[0]
   
    labels = watershed(-elev_map, markers, mask=local_masked_global)

    label_image = label(labels, connectivity = 2)
    

    
    areas = np.array([])
    regions_to_reconsider = []

    for region in regionprops(label_image):
        areas = np.append(areas, region.area)

    std = np.std(areas)
    mean = np.mean(areas)


    nuclei = np.array([])

    for region in regionprops(label_image):
        minr, minc, maxr, maxc = region.bbox
        
        if region.area >= thresh_up:
            regions_to_reconsider = np.append(regions_to_reconsider, region)
        elif region.area >= thresh_down:
            new_nucleus = Nucleus(region,len(nuclei))
            nuclei = np.append(nuclei, new_nucleus)      

    return nuclei, regions_to_reconsider, label_image


image = io.imread(dapi_path)

nuclei, regions_to_reconsider, thresholded = dapi_segmentation(image, 16, 50, 450, 39)#, sobel_g = True)
    
arc_im = io.imread(arc_path)
homer_im = io.imread(homer_path)

arc_locations, _,_ = ieg_segmentation(arc_im, 2, 7, 30, 21)
homer_locations, _,_ = ieg_segmentation(homer_im, 2, 7, 30, 21)



def prepare_coords_set(coords_tab):
    return set(tuple([m, n]) for m, n in zip(coords_tab[:,0], coords_tab[:,1]))

def find_ieg_colloc(ieg_locs, nuclei, to_reevaluate):
    ieg_positive = np.array([])

    for nucleus in nuclei:
        ieg_no = 0
        for ieg_dot in ieg_locs:
            nucleus_coords_set = prepare_coords_set(nucleus.nucelusProps.coords)
            ieg_coords_set = prepare_coords_set(ieg_dot.coords)      
            #if the centroid of area identified as immediate early gene presence is inside of nucleus, add it
            #to the number of ieg collocalized with the nucleus
            len_pre = len(nucleus_coords_set)
            
            set_pre = nucleus_coords_set.copy()

            nucleus_coords_set.add((int(np.floor(ieg_dot.centroid[0])), int(np.floor(ieg_dot.centroid[1]))))
            len_post = len(nucleus_coords_set)

            if (len_pre == len_post):
                ieg_no +=1
                if ieg_no == 1:
                    ieg_positive = np.append(ieg_positive, nucleus.idno)

        if ieg_no > 1:
            to_reevaluate = np.append(to_reevaluate, nucleus)
    return ieg_positive, to_reevaluate


to_reevalueate = []

arc_pos, arc_mult =  find_ieg_colloc(arc_locations, nuclei, to_reevalueate)
homer_pos, homer_mult =  find_ieg_colloc(homer_locations, nuclei, to_reevalueate)

print('nuclei with arc ', len(arc_pos))#, len(arc_mult))
print('nuclei with homer ', len(homer_pos))#,len(homer_mult))


double_ieg = np.intersect1d(arc_pos, homer_pos)
print("nuclei with both iegs ", len(double_ieg))


