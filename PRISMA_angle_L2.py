# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
import h5py
import numpy as np
import os

### read PRSIMA vaa
file_path = '/config_folder/angle_config.txt'

with open(file_path, 'r') as file:
  lines = file.readlines()
  file_content = ''.join(lines)

#print(file_content)

for file in os.listdir(file_content):
  if file.endswith(".he5") & file.startswith("PRS"):
    inputfile = os.path.join(file_content,file)

f = h5py.File(inputfile, 'r')

var = os.path.basename(inputfile)[4:7]

# print(var)

xs = ['PRS_', var, '_HCO']
s = ''.join(xs)

viewing_zenith_angle_mat = f['HDFEOS/SWATHS'][s]['Geometric Fields']['Observing_Angle']
relative_azimuth_angle_mat = f['HDFEOS/SWATHS'][s]['Geometric Fields']['Rel_Azimuth_Angle']
num_rows, num_cols = viewing_zenith_angle_mat.shape

center_row = num_rows // 2
center_col = num_cols // 2

viewing_zenith_angle = viewing_zenith_angle_mat[center_row,center_col]
relative_azimuth_angle = relative_azimuth_angle_mat[center_row,center_col]
sun_azimuth_angle = f.attrs['Sun_azimuth_angle']
sun_zenith_angle = f.attrs['Sun_zenith_angle']
viewing_azimuth_angle = relative_azimuth_angle + sun_azimuth_angle
date = f.attrs['Product_StartTime']
if viewing_azimuth_angle > 180 :
  viewing_azimuth_angle = sun_azimuth_angle - relative_azimuth_angle

print('sun_azimuth_angle')
print(sun_azimuth_angle)
print('sun_zenith_angle')
print(sun_zenith_angle)
print('viewing_azimuth_angle')
print(viewing_azimuth_angle)
print('viewing_zenith_angle')
print(viewing_zenith_angle)
print('date')
print(date)
# print('relative_azimuth_angle')
# print(relative_azimuth_angle)

