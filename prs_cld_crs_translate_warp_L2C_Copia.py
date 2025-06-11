#!/usr/bin/env python
# coding: utf-8

# In[1]:


from arosics import COREG
from arosics import COREG_LOCAL
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, osr, ogr
import pandas as pd
import os
import geopandas as gpd

# In[5]:


path = 'C:/Users/emast/Desktop/250606_agile/20230407/'


# In[6]:


im_target = path+'prs_crs.tif'
im_reference = path+'S2_20230329_B8_ritagliato_QGIS.tif'


# In[12]:


gdal.__version__


# In[ ]:


'''
VNIR
B2-492.3 = B12-489.8
B3-558.9 = B21-559.0
B4-664.9 = B33-664.9
B5-703.8 = B37-703.7
B6-739.1 = B40-734.0
B7-779.7 = B44-775.3
B8-832.9 = B49-827.9
B8A-864.0 = B52 - 860.0

SWIR
B11-1610.4 = B64-1606.5
B12-2185.7 = B128-183.4

'''


# In[7]:


CRL = COREG_LOCAL(im_reference,
                  im_target, 
                  grid_res = 10,
                  window_size =(128,128),
                  max_shift = 10,
                  #align_grids =False,
                  #match_gsd=False,
                  resamp_alg_deshift = 'nearest',
                  resamp_alg_calc = 'nearest',
                  path_out = path+'Registered_AROSICS_Local_B8_B52.tif',
                  fmt_out = 'GTIFF',
                  max_iter = 8,
                  #r_b4match = 1,
                  s_b4match = 51)


# In[8]:


CRL.correct_shifts()


# In[ ]:


CRL.view_CoRegPoints(figsize=(15,15))


# In[11]:





# In[10]:


points = CRL.CoRegPoints_table
points


# In[9]:


CRL.tiepoint_grid.to_PointShapefile(path_out=path+'AROSICS_TiePoints_Local_B8_B52.shp')


# In[ ]:





# In[ ]:


#====================================================================================================================
'''OPTIMIZAZTION TIE POINTS'''
#====================================================================================================================


# In[ ]:


#path = 'D:\\Documenti\\CNR-IBE\\Dati\\PRISMA_pansharpening_project\\GROSSETO\\AROSICS_output\\20200801\\'


# In[ ]:


lista_files = []
dirFileList = os.listdir(path) #legge la cartella 1
os.chdir(path) #setta la cartella
for file in dirFileList: #per ciascun file nella cartella 1
    if os.path.splitext(file)[-1] == '.shp': #le cui ultime lettere sono ".img"
        lista_files.append(os.path.join(path, file)) #aggiungile alla lista vuot

lista_files


# In[ ]:


order_bands = []
i = 0
while i < len(lista_files):
    if i == 0:
        f = gpd.read_file(lista_files[i])
        order_bands.append(lista_files[i].split("_")[-2])
        f['S2_Band'] = lista_files[i].split("_")[-2]
    if i > 0:
        f2 = gpd.read_file(lista_files[i])
        order_bands.append(lista_files[i].split("_")[-2])
        f2['S2_Band'] = lista_files[i].split("_")[-2]
        f = f.merge(f2,on='POINT_ID')
    
    i+=1


# In[ ]:


#print columns
np.unique(f.columns)


# In[ ]:


b = dict()

b['RELIABILITY'] = np.array(f['RELIABILIT'])# np.array(f[['RELIABILIT_x','RELIABILIT_y']])
b['POINT_ID'] = np.array(f['POINT_ID'])
b['X_IM'] = np.array(f['X_IM'])# np.array(f[['X_IM_x','X_IM_y']])
b['X_MAP'] = np.array(f['X_MAP'])# np.array(f[['X_MAP_x','X_MAP_y']])
b['Y_IM'] = np.array(f['Y_IM'])# np.array(f[['Y_IM_x','Y_IM_y']])
b['Y_MAP'] = np.array(f['Y_MAP'])# np.array(f[['Y_MAP_x','Y_MAP_y']])
b['S2_Band'] = np.array(f['S2_Band'])# np.array(f[['S2_Band_x','S2_Band_y']])

b['ABS_SHIFT'] = np.array(f['ABS_SHIFT'])# np.array(f[['ABS_SHIFT_x', 'ABS_SHIFT_y']])
b['ANGLE'] = np.array(f['ANGLE'])# np.array(f[['ANGLE_x', 'ANGLE_y']])
b['L1_OUTLIER'] = np.array(f['L1_OUTLIER'])# np.array(f[['L1_OUTLIER_x','L1_OUTLIER_y' ]])
b['L2_OUTLIER'] = np.array(f['L2_OUTLIER'])# np.array(f[['L2_OUTLIER_x','L2_OUTLIER_y' ]])
b['L3_OUTLIER'] = np.array(f['L3_OUTLIER'])# np.array(f[['L3_OUTLIER_x','L3_OUTLIER_y' ]])
b['LAST_ERR'] = np.array(f['LAST_ERR'])# np.array(f[['LAST_ERR_x', 'LAST_ERR_y']])
b['OUTLIER'] = np.array(f['OUTLIER'])# np.array(f[['OUTLIER_x', 'OUTLIER_y']])
b['REF_BADDAT'] = np.array(f['REF_BADDAT'])# np.array(f[['REF_BADDAT_x','REF_BADDAT_y' ]])
b['SSIM_AFTER'] = np.array(f['SSIM_AFTER'])# np.array(f[['SSIM_AFTER_x','SSIM_AFTER_y' ]])
b['SSIM_BEFOR'] = np.array(f['SSIM_BEFOR'])# np.array(f[['SSIM_BEFOR_x','SSIM_BEFOR_y' ]])
b['SSIM_IMPRO'] = np.array(f['SSIM_IMPRO'])# np.array(f[['SSIM_IMPRO_x','SSIM_IMPRO_y' ]])
b['TGT_BADDAT'] = np.array(f['TGT_BADDAT'])# np.array(f[['TGT_BADDAT_x','TGT_BADDAT_y' ]])
b['X_SHIFT_M'] = np.array(f['X_SHIFT_M'])# np.array(f[['X_SHIFT_M_x', 'X_SHIFT_M_y']])
b['X_SHIFT_PX'] = np.array(f['X_SHIFT_PX'])# np.array(f[['X_SHIFT_PX_x','X_SHIFT_PX_y' ]])
b['X_WIN_SIZE'] = np.array(f['X_WIN_SIZE'])# np.array(f[['X_WIN_SIZE_x','X_WIN_SIZE_y' ]])
b['Y_SHIFT_M'] = np.array(f['Y_SHIFT_M'])# np.array(f[['Y_SHIFT_M_x', 'Y_SHIFT_M_y']])
b['Y_SHIFT_PX'] = np.array(f['Y_SHIFT_PX'])# np.array(f[['Y_SHIFT_PX_x','Y_SHIFT_PX_y' ]])
b['Y_WIN_SIZE'] = np.array(f['Y_WIN_SIZE'])# np.array(f[['Y_WIN_SIZE_x','Y_WIN_SIZE_y' ]])

b['geometry'] = np.array(f['geometry'])# np.array(f[['geometry_x','geometry_y']])


# In[ ]:


#applichiamo 

x = np.argsort(b['RELIABILITY']) #feature prescelta per argsort()


RELIABILITY_sorted = np.take_along_axis(b['RELIABILITY'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
X_IM = np.take_along_axis(b['X_IM'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
X_MAP = np.take_along_axis(b['X_MAP'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
Y_IM = np.take_along_axis(b['Y_IM'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
Y_MAP  = np.take_along_axis(b['Y_MAP'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
S2_Band  = np.take_along_axis(b['S2_Band'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
geometry = np.take_along_axis(b['geometry'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)

ABS_SHIFT = np.take_along_axis(b['ABS_SHIFT'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
ANGLE= np.take_along_axis(b['ANGLE'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
L1_OUTLIER= np.take_along_axis(b['L1_OUTLIER'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
L2_OUTLIER= np.take_along_axis(b['L2_OUTLIER'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
L3_OUTLIER= np.take_along_axis(b['L3_OUTLIER'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
LAST_ERR= np.take_along_axis(b['LAST_ERR'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
OUTLIER= np.take_along_axis(b['OUTLIER'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
REF_BADDAT= np.take_along_axis(b['REF_BADDAT'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
SSIM_AFTER= np.take_along_axis(b['SSIM_AFTER'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
SSIM_BEFOR= np.take_along_axis(b['SSIM_BEFOR'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
SSIM_IMPRO= np.take_along_axis(b['SSIM_IMPRO'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
TGT_BADDAT= np.take_along_axis(b['TGT_BADDAT'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
X_SHIFT_M= np.take_along_axis(b['X_SHIFT_M'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
X_SHIFT_PX= np.take_along_axis(b['X_SHIFT_PX'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
X_WIN_SIZE= np.take_along_axis(b['X_WIN_SIZE'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
Y_SHIFT_M= np.take_along_axis(b['Y_SHIFT_M'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
Y_SHIFT_PX= np.take_along_axis(b['Y_SHIFT_PX'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)
Y_WIN_SIZE= np.take_along_axis(b['Y_WIN_SIZE'], x, axis=0) #ordinare in funzione di x (quindi della feature prescelta)

# In[ ]:


'''GET GCPS FROM AROSICS Tie Points TABLE'''

#apply shifting to the points
#applying the modificaiton during the dataframe construction


'''
se volessimo usare un dataframe esterno. Ma attenzione perchè delle volte fallisce gdal.Warp
a causa del dtype del dataframe

points = gpd.read_file(path+'AROSICS_TiePoints_Local_HigherReliability_20200801.shp')

gcps = pd.DataFrame({'X': pd.to_numeric(points['X_MAP'])+pd.to_numeric(points['X_SHIFT_M']), 
                     'Y': pd.to_numeric(points['Y_MAP'])+pd.to_numeric(points['Y_SHIFT_M']),
                     'Col': pd.to_numeric(points['X_IM']), 
                     'Row': pd.to_numeric(points['Y_IM'])}) 
'''

#quindi, utilizziamo direttamente le colonne dell'arrey sorted:
gcps = pd.DataFrame({'X': X_MAP[:,]+X_SHIFT_M[:,],
                     'Y': Y_MAP[:,]+Y_SHIFT_M[:,],
                     'Col': X_IM[:,],
                     'Row': Y_IM[:,]})

#gcps for PANCHROMATIC image: since it is at 5m instead of 30m, the position of the X/Y-IM must be multiplied for 6
gcps_PAN = pd.DataFrame({'X': X_MAP[:,]+X_SHIFT_M[:,],
                     'Y': Y_MAP[:,]+Y_SHIFT_M[:,],
                     'Col': X_IM[:,]*6,
                     'Row': Y_IM[:,]*6})


#create GDAL GCPs
gcps_gdal = [gdal.GCP(row['X'], row['Y'],0, row['Col'], row['Row']) for index, row in gcps.iterrows()] ### "gcps" or "gcps_PAN"


# In[ ]:


kwargs = {
    'format': 'GTiff'}
    #'outputType': gdal.GDT_UInt16}

output_image = path+'prs_crs_translate.tif'
ds_gcp = gdal.Translate(output_image, 
                        path+'prs_crs.tif',
                        outputSRS='EPSG:32632', 
                        GCPs=gcps_gdal,
                        **kwargs)


options = gdal.WarpOptions(dstSRS='EPSG:32632', polynomialOrder=2, targetAlignedPixels=True, xRes=30, yRes =30)
ds = gdal.Warp(path+'prs_crs_translate_warp.tif', ds_gcp, dstNodata = np.nan, options=options, resampleAlg="near")
ds_gcp = None
ds = None


# In[ ]:




