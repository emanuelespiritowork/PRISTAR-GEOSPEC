#def prs_cld_crs(prisma_FULL_path, S2_path):
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
#os.chdir("Z:/Progetto_PRISMA/PRISMA_code")

prisma_FULL_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\prove_per_pacchetto\PRS_L1_STD_OFFL_20250424100426_20250424100430_0001_HCO_FULL.tif"
S2_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\prove_per_pacchetto\S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif"
cloud_path = os.path.join(os.path.dirname(prisma_FULL_path), os.path.basename(prisma_FULL_path).replace("FULL","CLD"))

prisma_dataset = gdal.Open(prisma_FULL_path, gdal.GA_ReadOnly) #Read
prisma_arr = prisma_dataset.ReadAsArray
prisma_prj = prisma_dataset.GetProjection()
with rasterio.open(prisma_FULL_path) as prs:
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

#cloud = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\PRS_L1_STD_OFFL_20250424100426_20250424100430_0001_HCO_CLD.tif"
clouds = gdal.Open(cloud_path, gdal.GA_ReadOnly)
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

output_cloud = os.path.join(os.path.dirname(prisma_FULL_path),"cloud.tif") #r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\cloud.tif"
with rasterio.open(output_cloud, 'w', **metadata) as dst:
    dst.write(eroded_cld, 1)

eroded_cld = gdal.Open(output_cloud, gdal.GA_ReadOnly)
cloud = eroded_cld.ReadAsArray()
cloud_3d = cloud[:, :, np.newaxis]
merged_prisma_cld = np.concatenate((prisma_arr, cloud_3d), axis=2)
#combined_arr = np.dstack((prisma_arr, cloud_3d))
#merged_prisma_cld[merged_prisma_cld == -999] = np.nan

output_geotiff_path = os.path.join(os.path.dirname(prisma_FULL_path),"prs_cld.tif") #r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\prs_cld.tif"
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

merged_prisma_cloud_path = output_geotiff_path #r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\prs_cld.tif"
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

#S2_path = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif"
S2_dataset = gdal.Open(S2_path, gdal.GA_ReadOnly)
S2_prj = S2_dataset.GetProjection()

output_warp = os.path.join(os.path.dirname(prisma_FULL_path),"prs_cld_crs.tif") #r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\PRS_L1_STD_OFFL_20250424\prs_cld_crs.tif"
gdal.Warp(output_warp, merged_prisma_cloud, format='GTiff', dstSRS=S2_prj, resampleAlg="near", srcNodata=-999, xRes=30, yRes=30) #prisma_dataset
