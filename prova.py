
from arosics import COREG
from arosics import COREG_LOCAL
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, osr, ogr
import pandas as pd
import os
import geopandas as gpd

# In[5]:

# In[6]:

im_target = r"\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\prove_per_pacchetto\prs_cld_crs.tif"
im_reference = r'\\10.0.1.243\nr_data\3_rs_data\PRISMA\JDS\2025\L1\prove_per_pacchetto\S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif'
path = os.path.dirname(im_target)

print(im_target)
print(im_reference)
print(path)

AROSICS_out = os.path.join(path,'Registered_AROSICS_Local_B8_B52.tif')

CRL = COREG_LOCAL(im_reference,
                  im_target,
                  grid_res = 10,
                  window_size =(128,128),
                  max_shift = 10,
                  #align_grids =False,
                  #match_gsd=False,
                  resamp_alg_deshift = 'nearest',
                  resamp_alg_calc = 'nearest',
                  path_out = AROSICS_out,
                  fmt_out = 'GTIFF',
                  max_iter = 8,
                  #r_b4match = 1,
                  q = False,
                  s_b4match = 51)