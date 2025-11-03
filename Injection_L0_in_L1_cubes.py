#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import warnings
warnings.filterwarnings('ignore')
 
import os
from os.path import join
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as npatches
from matplotlib.dates import DateFormatter
#import xarray as xr
import numpy as np
import glob
from glob import iglob
import imageio
import datetime
import time
from osgeo import gdal, osr
import h5py
import pandas as pd  


# In[ ]:


#=========================================================
#LEVEL 0
#=========================================================


# In[ ]:


# open the PRISMA file
path = "/space/put_PRISMA_he5_and_S2_tif_here/"
f_vnir = h5py.File(path+'data_SubAcq3_C_VNIR_SURFACE-OBS_Part0_S11.h5', 'r')#prendo .he5 che contiene VNIR ed S11 nel nome
f_swir = h5py.File(path+'data_SubAcq3_C_SWIR_SURFACE-OBS_Part0_S11.h5', 'r')#prendo .he5 che contiene SWIR ed S11 nel nome 


# In[ ]:


def printname(name):
    print(name)
f_vnir.visit(printname)
print('\n')
f_swir.visit(printname)


# In[ ]:


vnir_L0 = f_vnir['data']
swir_L0 = f_swir['data']

cw_list_vnir = f_vnir['Spectral_info/CW_list']
cw_list_swir = f_swir['Spectral_info/CW_list']

cw_matirx_vnir = f_vnir['Spectral_info/CW_matrix']
cw_matirx_swir = f_swir['Spectral_info/CW_matrix']


# In[ ]:


#maskout black edges (zeros-nodata)
vnir_L0 = vnir_L0[:1000,:,:]
swir_L0 = swir_L0[:1000,:,:]


# In[ ]:


indices_vnir = np.where(cw_list_vnir[...]>0)[0]
indices_swir = np.where(cw_list_swir[...]>0)[0]

#skip di mask dello swir
indices_swir = indices_swir[4:]


# In[ ]:


vnir_L0_filtered = vnir_L0[:, indices_vnir, :]
swir_L0_filtered = swir_L0[:, indices_swir, :]


# In[ ]:


vnir_L0_filtered = np.transpose(vnir_L0_filtered, (2,1,0))
swir_L0_filtered = np.transpose(swir_L0_filtered, (2,1,0))


# In[ ]:


plt.imshow(swir_L0_filtered[:,100,:])


# In[ ]:





# In[ ]:





# In[ ]:


#=========================================================
#'''REPLACE CUBES - L0 in L1'''
#=========================================================


# In[ ]:


# open the PRISMA file
#path = "C:/Users/Emanuele/Desktop/PRISMA_code/put_PRISMA_he5_and_S2_tif_here/"
for file in os.listdir(path):
        if file.endswith(".he5") & file.startswith("PRS"):
            inputfile = os.path.join(path,file)
f = h5py.File(inputfile, 'r+')#prendo il file che contiene he5 e PRS_L1

# reading name and value for root attributes (metadata contained in HDF5root)
for attribute in f.attrs:
    print(attribute,f.attrs[attribute])

#create the metadata dictionary
metadata_dic={}
for attribute in f.attrs:
    metadata_dic[attribute] = f.attrs[attribute]


# In[ ]:


central_wavelenght_VNIR = np.flip(metadata_dic['List_Cw_Vnir'])
central_wavelenght_SWIR = np.flip(metadata_dic['List_Cw_Swir'])
fwhm_wavelenght_VNIR = np.flip(metadata_dic['List_Fwhm_Vnir'])
fwhm_wavelenght_SWIR = np.flip(metadata_dic['List_Fwhm_Swir'])

ScaleFactor_PAN = np.flip(metadata_dic['ScaleFactor_Pan'])
ScaleFactor_SWIR = np.flip(metadata_dic['ScaleFactor_Swir'])
ScaleFactor_VNIR  = np.flip(metadata_dic['ScaleFactor_Vnir'])

Offset_PAN = np.flip(metadata_dic['Offset_Pan'])
Offset_SWIR = np.flip(metadata_dic['Offset_Swir'])
Offset_VNIR  = np.flip(metadata_dic['Offset_Vnir'])


# In[ ]:


def printname(name):
    print(name)
f.visit(printname)


# In[ ]:


swir = f['/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SWIR_Cube']
vnir = f['/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/VNIR_Cube']
pan = f['HDFEOS/SWATHS/PRS_L1_PCO/Data Fields/Cube']

lat = f['/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Latitude_VNIR']
lon = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Longitude_VNIR']
geo_field = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields']

lat_pan = f['/HDFEOS/SWATHS/PRS_L1_PCO/Geolocation Fields/Latitude']
lon_pan = f['HDFEOS/SWATHS/PRS_L1_PCO/Geolocation Fields/Longitude']
geo_field_pan = f['HDFEOS/SWATHS/PRS_L1_PCO/Geolocation Fields']


# In[ ]:


plt.imshow(swir[:, 100,:]) #show band 6


# In[ ]:


swir[...]= swir_L0_filtered * ScaleFactor_SWIR
vnir[...]= vnir_L0_filtered * ScaleFactor_VNIR


# In[ ]:


plt.imshow(swir[:, 100,:]) #show band 6
import re
#os.rename(inputfile, re.sub('.he5$',"_injected.he5",inputfile))
