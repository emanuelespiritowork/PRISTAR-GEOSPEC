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
os.chdir(r"\\10.0.1.243\nr_working\emanuele\Progetto_PRISMA\PRISMA_code")

path = "C:/Users/emast/Desktop/250606_agile/20230304/"
prisma_path = path+"PRS_L2C_STD_20230304102047_20230304102051_0001_HCO_FULL.tif"
cloud_path = path+"PRS_L2C_STD_20230304102047_20230304102051_0001_HCO_ANG.tif"
S2_path = path+"S2_20230309_B08_T32TQQ_ritagliato_QGIS.tif"

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


clouds = gdal.Open(cloud_path, gdal.GA_ReadOnly)
cloud = clouds.ReadAsArray()
cloud_3d = np.transpose(cloud, (1, 2, 0))
merged_prisma_cld = np.concatenate((prisma_arr, cloud_3d), axis=2)
#combined_arr = np.dstack((prisma_arr, cloud_3d))
#merged_prisma_cld[merged_prisma_cld == -999] = np.nan

output_geotiff_path = path+"prs_cld.tif"
crs = prs.crs.to_wkt()
# Create a rasterio Profile
profile = {
    'driver': 'GTiff',
    'count': 233,  # Number of bands
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

merged_prisma_cloud_path = path+"prs_cld.tif"
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

output_warp = path+"prs_cld_crs.tif"
gdal.Warp(output_warp, merged_prisma_cloud, format='GTiff', dstSRS=S2_prj, resampleAlg="near", srcNodata=-999, xRes=30, yRes=30) #prisma_dataset

