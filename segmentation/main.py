#!/usr/bin/env python
# coding: utf-8

import os
import catfish_dapi as dapi
from skimage import io


from common import save_approximations_overlay

scan_path_root = "/home/ula/catfish/test_scans/"
dapi_suffix = "/1.tif"
overlay_suffix = "/overlay.tif"

#MAIN loop to walk through img folder, looking for FOLDERS
for filename in os.listdir(scan_path_root):
    if os.path.isdir(scan_path_root + filename):
        print(scan_path_root + filename)
        dapi_path = scan_path_root + filename + dapi_suffix
        dapi_img = io.imread(dapi_path)
        nuclei, ellipse_axes_ratios, areas, intensities = dapi.dapi_segmentation(dapi_img,10, 75, 1250, 53)
        save_approximations_overlay(dapi_img, nuclei, scan_path_root + filename + overlay_suffix)
