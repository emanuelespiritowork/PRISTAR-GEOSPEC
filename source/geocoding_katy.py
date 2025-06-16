import numpy as np
import pandas as pd
from osgeo import gdal
import geopandas as gpd
import rasterio
#import h5py
import os
import numpy as np
from rasterio import Affine
from rasterio.warp import reproject
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
from scipy.ndimage import binary_dilation, binary_erosion
import matplotlib.pyplot as plt
from rasterio.features import geometry_mask
import fiona
import rasterio.mask
from rasterio.mask import mask
from rasterio import features
import matplotlib
from matplotlib import patches as mpatches, colors
from matplotlib.pyplot import imshow
from datetime import datetime
from pathlib import Path
import pyproj
import math
from glob import glob
from rasterio import plot
from geoarray import GeoArray
import spectral
from arosics import COREG_LOCAL
from rasterio.warp import calculate_default_transform, reproject, Resampling
os.chdir("C:/Users/Katayoun/OneDrive/Desktop")

#prisma_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\Braccagni\2023\L1\PRS_L1_STD_OFFL_20231004101115_20231004101119_0001_HCO_FULL.tif"
prisma_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20240928\PRS_L1_STD_OFFL_20240928101742_20240928101746_0001_HCO_FULL.tif"
prisma_dataset = gdal.Open(prisma_path, gdal.GA_ReadOnly) #Read
prisma_arr = prisma_dataset.ReadAsArray
prisma_prj = prisma_dataset.GetProjection()
with rasterio.open(prisma_path) as prs:
    print("Number of bands:", prs.count)
    print("Width:", prs.width)
    print("Height:", prs.height)
    print("CRS:", prs.crs)
    print("Transform (Affine):", prs.transform)
    prisma_arr = prs.read()
    print("shape:", prisma_arr.shape)

prisma_arr = np.transpose(prisma_arr, (1, 2, 0))
plt.imshow(prisma_arr[:, :, 51])
plt.title('prisma')
plt.show()
print("shape:", prisma_arr.shape)

cloud = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20240928\PRS_L1_STD_OFFL_20240928101742_20240928101746_0001_HCO_CLD.tif"
clouds = gdal.Open(cloud, gdal.GA_ReadOnly)
data_cloud = clouds.ReadAsArray()
#data_cloud[data_cloud == -999] = np.nan
#data_cloud[data_cloud > 0] = 1

kernel = np.array([[0, 0, 0],
                  [1, 1, 1],
                  [0, 0, 0]])
dilated_cld = binary_dilation(data_cloud, structure=kernel).astype(np.uint8)
eroded_cld = binary_erosion(data_cloud, structure=kernel).astype(np.uint8)

plt.figure(figsize=(7, 3))
plt.subplot(131)
plt.imshow(data_cloud)  #, cmap='gray'
plt.title('Original Image')
plt.subplot(132)
plt.imshow(dilated_cld)
plt.title('Dilated Image')
plt.subplot(133)
plt.imshow(eroded_cld)
plt.title('Eroded Image')
plt.show()

metadata = {
    'driver': 'GTiff',
    'count': 1,
    'dtype': 'uint8',
    'width': clouds.RasterXSize,
    'height': clouds.RasterYSize,
    'crs': prisma_dataset.GetProjection(),
    'transform': Affine.from_gdal(*clouds.GetGeoTransform())
}

eroded_cld[:, 0] = 1   # Set 1 in the first column for all bands
eroded_cld[:, -1] = 1  # Set 1 in the last column for all bands

output_file = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20240928\cloud.tif"
with rasterio.open(output_file, 'w', **metadata) as dst:
    dst.write(eroded_cld, 1)

eroded_cld = gdal.Open(output_file, gdal.GA_ReadOnly)
cloud = eroded_cld.ReadAsArray()
cloud_3d = cloud[:, :, np.newaxis]
merged_prisma_cld = np.concatenate((prisma_arr, cloud_3d), axis=2)
#combined_arr = np.dstack((prisma_arr, cloud_3d))
#merged_prisma_cld[merged_prisma_cld == -999] = np.nan

output_geotiff_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20241113\prs_cld.tif"
crs = prs.crs.to_wkt()
# Create a rasterio Profile
profile = {
    'driver': 'GTiff',
    'count': 231,  # Number of bands
    'dtype': merged_prisma_cld.dtype,  # Data type of the array
    'width': merged_prisma_cld.shape[1],  # Width of the image
    'height': merged_prisma_cld.shape[0],  # Height of the image
    'crs': crs,  # Coordinate Reference System
    'transform': prs.transform,
}

with rasterio.open(output_geotiff_path, 'w', **profile) as dst:
    for band in range(profile['count']):
        band_data = merged_prisma_cld[:, :, band]
        dst.write(band_data, indexes=band + 1)

print("GeoTIFF file has been created:", output_geotiff_path)

merged_prisma_cloud_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20240928\prs_cld.tif"
merged_prisma_cloud = gdal.Open(merged_prisma_cloud_path, gdal.GA_ReadOnly)
merged_prisma_cloud_ar = merged_prisma_cloud.ReadAsArray
merged_prisma_cloud_prj = merged_prisma_cloud.GetProjection()
with rasterio.open(merged_prisma_cloud_path) as mrg:
    print("Number of bands:", mrg.count)
    print("Width:", mrg.width)
    print("Height:", mrg.height)
    print("CRS:", mrg.crs)
    print("Transform (Affine):", mrg.transform)
    merged_prisma_cloud_ar = mrg.read()
merged_prisma_cloud_ar = np.transpose(merged_prisma_cloud_ar, (1, 2, 0))

S2_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20240928\S2_20240929_10m_B08.tif"
S2_dataset = gdal.Open(S2_path, gdal.GA_ReadOnly)
S2_data = S2_dataset.ReadAsArray
S2_prj = S2_dataset.GetProjection()
with rasterio.open(S2_path) as src:
    print("Number of bands:", src.count)
    print("Width:", src.width)
    print("Height:", src.height)
    print("CRS (Coordinate Reference System):", src.crs)
    print("Transform (Affine):", src.transform)
    S2_data = src.read()

output_warp = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20241113\prs_cld_warp.tif"
gdal.Warp(output_warp, merged_prisma_cloud, format='GTiff', dstSRS=S2_prj, resampleAlg="near", srcNodata=-999, xRes=30, yRes=30) #prisma_dataset


# S2_data_SR = (S2_data + (-1000)) / 10000
#
# metadata = {
#     'driver': 'GTiff',
#     'count': 1,
#     'dtype': 'float32',
#     'width': src.width,
#     'height': src.height,
#     'crs': src.crs,
#     'transform': src.transform
# }
#
# output_file = "//10.0.1.252/projects/2022_ASI-SCIA/5_SHARED/Torgnon_Strisciata_Mar2024/PlateauRosa/S2_20240129_B08_10m.tif"
# # Write processed data to GeoTIFF file
# with rasterio.open(output_file, 'w', **metadata) as dst:
#     for band_idx in range(S2_data_SR.shape[0]):
#         dst.write(S2_data_SR[band_idx, :, :], band_idx + 1)  # Write each band individually
#
#
# src_nodata = 0.0
# dst_nodata = -999
# prisma_path = "//10.0.1.252/projects/2022_ASI-SCIA/5_SHARED/Torgnon_Strisciata_Mar2024/PlateauRosa/prisma_cld_merged.tif"
# prisma_dataset = gdal.Open(prisma_path, gdal.GA_ReadOnly)
# prisma_arr = prisma_dataset.ReadAsArray
# prisma_prj = prisma_dataset.GetProjection()
# with rasterio.open(prisma_path) as prs:
#     print("Number of bands:", prs.count)
#     print("Width:", prs.width)
#     print("Height:", prs.height)
#     print("CRS:", prs.crs)
#     print("Transform (Affine):", prs.transform)
#     prisma_arr = prs.read()
# prisma_arr = np.transpose(prisma_arr, (1, 2, 0))
# plt.imshow(prisma_arr[:, :, 1])
# plt.title('prisma')
# plt.show()
# output_warp = "C:/Users/Katayoun/OneDrive/Desktop/prisma_warped.tif"
# gdal.Warp(output_warp, prisma_dataset, format='GTiff', dstSRS=S2_prj, resampleAlg="near", srcNodata=-999, xRes=30, yRes=30)

#merged_prisma_cloud
output_warp = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\Lodi\PRS_L1_STD_OFFL_20200419103351_20200419103355_0001\prs_warp.tif"
gdal.Warp(output_warp, merged_prisma_cloud, format='GTiff', dstSRS=S2_prj, resampleAlg="near", srcNodata=-999, xRes=30, yRes=30) #prisma_dataset

x = merged_prisma_cloud_ar - merged_prisma_cld
print(np.unique(x))

im_tar = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2024\L1\PRS_L1_STD_OFFL_20241113\prs_cld_warp.tif"
im_ref = S2_path
data_s2 = gdal.Open(im_ref, gdal.GA_ReadOnly)
S2_prj = data_s2.GetProjection()
data_prs = gdal.Open(im_tar, gdal.GA_ReadOnly)
prs_prj = data_prs.GetProjection()

def main():

    kwargs = {
        'grid_res': 10,
        'window_size': (128, 128),#(256, 256), #(128, 128),
        'path_out': 'C:/Users/Katayoun/OneDrive/Desktop/',
        'fmt_out': 'GTIFF',
        'max_iter': 8,
        'projectDir': 'coreg',
        'max_shift': 10,
        'resamp_alg_deshift': 'nearest',
        's_b4match': 51,
        'q': False,
    }

      # 'resamp_alg_calc': 'nearest',
      # 'tieP_filter_level': 2,


    CRL = COREG_LOCAL(im_ref, im_tar, **kwargs)
    corr = CRL.correct_shifts()
    CRL.view_CoRegPoints(figsize=(15, 15), shapes2plot='points', backgroundIm='tgt', savefigPath='C:/Users/Katayoun/OneDrive/Desktop/plot1.png') #hide_filtered=False, exclude_fillVals=False
    CRL.view_CoRegPoints(figsize=(15, 15), attribute2plot='ABS_SHIFT', exclude_fillVals=False, hide_filtered=False, savefigPath='C:/Users/Katayoun/OneDrive/Desktop/plot11.png')

    #CRL.view_CoRegPoints(figsize=(15, 15), attribute2plot='ABS_SHIFT', exclude_fillVals=False, hide_filtered=False)
    CRL_after_corr = COREG_LOCAL(im_ref, CRL.path_out, **kwargs)
    CRL_after_corr.view_CoRegPoints(figsize=(15, 15), shapes2plot='points', backgroundIm='tgt', savefigPath='C:/Users/Katayoun/OneDrive/Desktop/plot2.png')
    CRL_after_corr.view_CoRegPoints(figsize=(15, 15), shapes2plot='points', backgroundIm='tgt', hide_filtered=False, savefigPath='C:/Users/Katayoun/OneDrive/Desktop/plot22.png')

    CRL.CoRegPoints_table
    CRL.tiepoint_grid.to_PointShapefile(path_out='C:/Users/Katayoun/OneDrive/Desktop/tie_10.shp')

if __name__ == '__main__':
    main()

