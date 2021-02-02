#!/bin/bash
 
/usr/bin/convert $1/$2.tif[0] $1/split_pages/$2_page1.tif
/usr/bin/convert $1/$2.tif[1] $1/split_pages/$2_page2.tif
/usr/bin/convert $1/$2.tif[2] $1/split_pages/$2_page3.tif 
