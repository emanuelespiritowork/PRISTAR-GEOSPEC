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
prisma_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\PRS_L1_STD_OFFL_20250424100426_20250424100430_0001_HCO_FULL.tif"
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

