#!/usr/bin/env python3

"""
Temporal coregistration of Earth Observation data using AROSICS algorithm
"""

__author__ = 'Federico Filipponi'
__copyright__ = 'CC BY-NC-SA 2021-2025 Federico Filipponi'
__credits__ = ['Federico Filipponi']
__license__ = 'GPLv3'
__version__ = '1.1'
__date__ = '2025-05-02'
__maintainer__ = 'Federico Filipponi'
__contact__ = 'federico.filipponi@gmail.com'
__status__ = 'Development'

"""
##################################
# GPLv3 disclaimer:

Copyright (CC BY-NC-SA) 2021-2025 Federico Filipponi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""

###########################################################################
# load libraries
import os
import sys
import argparse
import warnings
import json
from arosics import COREG_LOCAL
from arosics.DeShifter import DESHIFTER
from py_tools_ds.geo.map_info import geotransform2mapinfo
from osgeo import gdal
import pandas as pd
import geopandas

###########################################################################
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
 
###########################################################################
# define functions

def to_GCP_gpd(self):
    # get copy of tie points grid without no data
    try:
        GDF = self.CoRegPoints_table.loc[self.CoRegPoints_table.ABS_SHIFT != self.outFillVal, :].copy()
    except AttributeError:
        # self.CoRegPoints_table has no attribute 'ABS_SHIFT' because all points have been excluded
        return []
    if getattr(GDF, 'empty'):  # GDF.empty returns AttributeError
        return []
    else:
        # exclude all points flagged as outliers
        if 'OUTLIER' in GDF.columns:
            GDF = GDF[GDF.OUTLIER.__eq__(False)].copy()
        avail_TP = len(GDF)
        if not avail_TP:
            # no point passed all validity checks
            return []
        if avail_TP > 7000:
            GDF = GDF.sample(7000)
            warnings.warn('By far not more than 7000 tie points can be used for warping within a limited '
                          'computation time (due to a GDAL bottleneck). Thus these 7000 points are randomly chosen '
                          'out of the %s available tie points.' % avail_TP)
        # calculate GCPs
        GDF['X_MAP_new'] = GDF.X_MAP + GDF.X_SHIFT_M
        GDF['Y_MAP_new'] = GDF.Y_MAP + GDF.Y_SHIFT_M
        return GDF

def coreg_info_base(crl) -> dict:
    """Return a dictionary containing general raster information to correct the detected local displacements of the target image."""
    if crl._coreg_info:
        return crl._coreg_info
    else:
        if not crl._tiepoint_grid:
            crl.calculate_spatial_shifts()
        TPG = crl._tiepoint_grid

        base_coreg_info = {
            # 'GCPList': gcp,
            'mean_shifts_px': {'x': TPG.mean_x_shift_px if TPG.GCPList else None,
                               'y': TPG.mean_y_shift_px if TPG.GCPList else None},
            'mean_shifts_map': {'x': TPG.mean_x_shift_map if TPG.GCPList else None,
                                'y': TPG.mean_y_shift_map if TPG.GCPList else None},
            'updated map info means': crl._get_updated_map_info_meanShifts() if TPG.GCPList else None,
            'original map info': geotransform2mapinfo(crl.imref.gt, crl.imref.prj),
            'reference projection': crl.imref.prj,
            'reference geotransform': crl.imref.gt,
            'reference grid': [[crl.imref.gt[0], crl.imref.gt[0] + crl.imref.gt[1]],
                               [crl.imref.gt[3], crl.imref.gt[3] + crl.imref.gt[5]]],
            'reference extent': {'cols': crl.imref.xgsd, 'rows': crl.imref.ygsd},  # FIXME not needed anymore
            'success': crl.success
        }
        return base_coreg_info

def GetExtentGPS(ds):
    """ Return DataFrame of GPS of corner coordinates from a gdal Dataset extent """
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    df = pd.DataFrame([[0,0,xmin,ymax],[0,height,xmin,ymin],[width,0,xmax,ymax],[width,height,xmax,ymin]], columns=['X_IM', 'Y_IM', 'X_MAP_new', 'Y_MAP_new'])
    #return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)
    return df
 
###########################################################################
# parse command line arguments

if len(sys.argv) == 1:
    prog = os.path.basename(sys.argv[0])
    print( '      '+sys.argv[0]+' [options]')
    print("      Help : ", prog, " --help")
    print("        or : ", prog, " -h")
    print("")
    print("      Coregister image using reference image, save GCP and coregistration information:")
    print("      python  %s --reference reference_raster.tif --target input_raster_to_coregister.tif --input input_raster_to_coregister.tif --output output_coregistered_raster.tif --output_gcp gcp.txt --output_cor_info coregistration_info.json"%sys.argv[0])
    print("")
    print("      Coregister image using reference image and save GCP, coregistration information and extra information (resampled GCP, tie-point grid in GeoPakage format):")
    print("      python  %s --reference reference_raster.tif --target input_raster_to_coregister.tif --input input_raster_to_coregister.tif --output output_coregistered_raster.tif --output_gcp gcp.txt --output_cor_info coregistration_info.json --output_gcp_resampled gcp_resampled.txt --output_tpg tiepoint_grid.gpkg"%sys.argv[0])
    print("")
    print("      Only compute GCP and coregistration information:")
    print("      python  %s --reference reference_raster.tif --target input_raster_to_coregister.tif --output_gcp gcp.txt --output_cor_info coregistration_info.json"%sys.argv[0])
    print("")
    print("      Coregister image using GCP and coregistration information:")
    print("      python  %s --input input_raster_to_coregister.tif --output output_coregistered_raster.tif --input_gcp gcp.txt --input_cor_info coregistration_info.json"%sys.argv[0])
    print("")
    sys.exit(-1)
else:
    
    parser = argparse.ArgumentParser()
    
    usage = "usage: %prog [options] "
    # parser = OptionParser(usage=usage)
    
    parser.add_argument("-i","--input", dest="input_file", action="store", type=str,  \
            help="Path of input raster file to be registered", default=None)
    parser.add_argument("-o","--output", dest="output_file", action="store", type=str,  \
            help="Output coregistered raster file path", default=None)
    parser.add_argument("-j","--input_cor_info", dest="input_cor_info", action="store", type=str,  \
            help="Input coregistration information", default=None)
    parser.add_argument("-k","--input_gcp", dest="input_gcp", action="store", type=str,  \
            help="Input GCP file name", default=None)
    parser.add_argument("-r","--reference", dest="reference", action="store", type=str,  \
            help="Reference raster file path",default=None)
    parser.add_argument("-t","--target", dest="target", action="store", type=str,  \
            help="Target raster file path", default=None)
    parser.add_argument("-l","--output_cor_info", dest="output_cor_info", action="store", type=str,  \
            help="Output coregistration information", default=None)
    parser.add_argument("-m","--output_gcp", dest="output_gcp", action="store", type=str,  \
            help="Output GCP file name", default=None)
    parser.add_argument("-n","--output_gcp_resampled", dest="output_gcp_resampled", action="store", type=str,  \
            help="Output resampled GCP file name", default=None)
    parser.add_argument("-g","--output_tpg", dest="output_tpg", action="store", type=str,  \
            help="Output GeoPackage file name containing tie-point grid", default=None)
    parser.add_argument("-q","--cores", dest="cores", action="store", type=int,  \
            help="Number of cores to use (default all cores)", default=None)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true",  \
            help="Overwrite GCP and coregistration info file is already existing", default=False)
    parser.add_argument("-v","--verbose", dest="verbose", action="store_true",  \
            help="Verbose mode", default=False)
    
    # check input arguments
    # (options, args) = parser.parse_args()
    options = parser.parse_args()

    coregister_image = False
    if options.reference == None or options.target == None:
        if options.input_gcp == None or options.input_cor_info == None:
            print("ERROR: Either input reference ('-r') or target raster ('-t') file or input GCP ('-k') or coregistration info ('-j') not provided.")
            sys.exit(1)
        else:
            if not os.path.exists(options.input_gcp):
                print("ERROR: Input GCP file set using argument '-k' does not exist. Set another input GCP file path.")
                sys.exit(1)
            if not os.path.exists(options.input_cor_info):
                rint("ERROR: Input coregistration info file set using argument '-j' does not exist. Set another input coregistration info file path.")
                sys.exit(1)
            compute_gcp = False
            if options.input_file == None or options.output_file == None:
                print("ERROR: Either input ('-i') or output ('-o') raster file not provided.")
                sys.exit(1)
            else:
                if not os.path.exists(options.input_file):
                    print("ERROR: Input raster file to be registered set using argument '-i' does not exist. Set another input raster file path.")
                    sys.exit(1)
                if os.path.exists(options.output_file):
                    print("ERROR: Output coregistered raster file set using argument '-o' already exists. Set another output raster file path.")
                    sys.exit(1)
                else:
                    output_dir = os.path.dirname(options.output_file)
                    if not os.path.exists(output_dir):
                        try:
                            os.mkdir(output_dir)
                        except OSError:
                            os.makedirs(output_dir, exist_ok=True)
                        except OSError:
                            print ("Creation of the output directory %s failed" % output_dir)
                            sys.exit(1)
                    outdir_writable = os.access(output_dir, os.W_OK)
                    if outdir_writable==False:
                        print("ERROR: Output file path set using argument '-o' is not writable.")
                        sys.exit(1)
                coregister_image = True
    else:
        if not os.path.exists(options.reference):
            print("ERROR: Input reference raster file to be registered set using argument '-r' does not exist. Set another input reference raster file path.")
            sys.exit(1)
        if not os.path.exists(options.target):
            print("ERROR: Input targer raster file set using argument '-t' does not exist. Set another input target raster file path.")
            sys.exit(1)
        if options.input_file == None or options.output_file == None:
            if options.output_gcp == None or options.output_cor_info == None:
                print("ERROR: Either input ('-i') or output ('-o') raster file or output GCP ('-m') or coregistration info ('-l') provided.")
                sys.exit(1)
            else:
                if os.path.exists(options.output_gcp) and options.overwrite == False:
                    print("ERROR: Output GCP file set using argument '-m' already exists. Set another output GCP file path.")
                    sys.exit(1)
                else:
                    output_dir = os.path.dirname(options.output_gcp)
                    if not os.path.exists(output_dir):
                        try:
                            os.mkdir(output_dir)
                        except OSError:
                            os.makedirs(output_dir, exist_ok=True)
                        except OSError:
                            print ("Creation of the output directory %s failed" % output_dir)
                            sys.exit(1)
                    outdir_writable = os.access(output_dir, os.W_OK)
                    if outdir_writable==False:
                        print("ERROR: Output GCP file path set using argument '-m' is not writable.")
                        sys.exit(1)
                if os.path.exists(options.output_cor_info) and options.overwrite == False:
                    print("ERROR: Output coregistration info file set using argument '-l' already exists. Set another output coregistration info file path.")
                    sys.exit(1)
                else:
                    output_dir = os.path.dirname(options.output_cor_info)
                    if not os.path.exists(output_dir):
                        try:
                            os.mkdir(output_dir)
                        except OSError:
                            os.makedirs(output_dir, exist_ok=True)
                        except OSError:
                            print ("Creation of the output directory %s failed" % output_dir)
                            sys.exit(1)
                    outdir_writable = os.access(output_dir, os.W_OK)
                    if outdir_writable==False:
                        print("ERROR: Output coregistration info file path set using argument '-l' is not writable.")
                        sys.exit(1)
            coregister_image = False
        else:
            if not os.path.exists(options.input_file):
                print("ERROR: Input raster file to be registered set using argument '-i' does not exist. Set another input raster file path.")
                sys.exit(1)
            if os.path.exists(options.output_file):
                print("ERROR: Output coregistered raster file set using argument '-o' already exists. Set another output raster file path.")
                sys.exit(1)
            else:
                output_dir = os.path.dirname(options.output_file)
                if not os.path.exists(output_dir):
                    try:
                        os.mkdir(output_dir)
                    except OSError:
                        os.makedirs(output_dir, exist_ok=True)
                    except OSError:
                        print ("Creation of the output directory %s failed" % output_dir)
                        sys.exit(1)
                outdir_writable = os.access(output_dir, os.W_OK)
                if outdir_writable==False:
                    print("ERROR: Output file path set using argument '-o' is not writable.")
                    sys.exit(1)
            coregister_image = True
        compute_gcp = True

###########################################################################
# set quiet and verbose arguments
verbose = options.verbose
quiet = not(verbose)

###########################################################################
# import coregistration information previously calculated

if compute_gcp == False:
    
    # import GCP from file
    GDF = pd.read_table(options.input_gcp, header=None, index_col=None, sep=" ")
    GDF.columns = ['X_IM', 'Y_IM', 'X_MAP_new', 'Y_MAP_new']

    # convert DataFrame GCP to GDAL GCP
    GDF['GCP'] = GDF.apply(lambda GDF_row: gdal.GCP(GDF_row.X_MAP_new,
                                                    GDF_row.Y_MAP_new,
                                                    0,
                                                    GDF_row.X_IM,
                                                    GDF_row.Y_IM),
                           axis=1)

    # initialize dictionary with full information for coregistration
    cor_info = {
        'GCPList': GDF.GCP.tolist()
    }

    # import coregistration information from json file
    cor_info_file = open(options.input_cor_info, 'rb')
    # update dictionary
    cor_info.update(json.load(cor_info_file))
    cor_info_file.close()

###########################################################################
# compute coregistration information from input image
if compute_gcp == True:

    # set coregistration arguments
    kwargs = {
        'grid_res'     : 10,#modificato
        'window_size'  : (1280,1280),#modificato
        'max_shift'    : 1000,#modificato (era 10 per Katy),
        'max_points'   : 100,
        'resamp_alg_deshift': "nearest",#non c'era
        'resamp_alg_calc': 'nearest',#non c'era
        'path_out'     : 'auto',
        'projectDir'   : 'my_project',
        'min_reliability' : 40,
        'ignore_errors' : True,
        'CPUs'         : options.cores,
        'q'            : quiet,
    }
    
    # compute local shifts
    CRL = COREG_LOCAL(options.reference,options.target,**kwargs)
    
    # clean GCP and convert to GeoDataFrame Object
    GDF = to_GCP_gpd(CRL)
    
    r_ds = gdal.Open(options.target)
    
    # get baseline from metadata
    r_ds_meta = r_ds.GetMetadata()
    if r_ds_meta.get('PROCESSING_BASELINE') == None:
        tgt_baseline = 'N0000'
    else:
        tgt_baseline = r_ds_meta.get('PROCESSING_BASELINE')
        tgt_baseline = ''.join(('N',tgt_baseline[0:2],tgt_baseline[3:5]))
   
    # convert GeoDataFrame GCP to GDAL GCP
    if len(GDF) == 0:
        # print("ERROR: No available GCP after cleaning. Exiting.")
        # sys.exit(1)
        
        print("WARNING: No available GCP after cleaning. Generating GCP using extent coordinates.")
        options.output_tpg = None

        # r_ds = gdal.Open(options.target)
        GDF = GetExtentGPS(r_ds)

    GDF['GCP'] = GDF.apply(lambda GDF_row: gdal.GCP(GDF_row.X_MAP_new,
                                                        GDF_row.Y_MAP_new,
                                                        0,
                                                        GDF_row.X_IM,
                                                        GDF_row.Y_IM),
                               axis=1)
    # GCPL = GDF.GCP.tolist()

    # generate dictionary with full information for coregistration
    cor_base_info = coreg_info_base(CRL)
    cor_info = {
        'GCPList': GDF.GCP.tolist()
    }
    cor_info.update(cor_base_info)
    
    cor_baseline = {
        'PROCESSING_BASELINE': tgt_baseline
    }
    cor_base_info.update(cor_baseline)

    # export GCP to text file
    if options.output_gcp != None:

        out_create = True
        if os.path.exists(options.output_gcp):
            if options.overwrite == False:
                print("WARNING: Output GCP file set using argument '-m' already exists. Skipping GCP file creation.")
                out_create = False
            else:
                print("WARNING: Output GCP file set using argument '-m' already exists. Overwriting file.")

            if out_create == True:
                output_dir = os.path.dirname(options.output_gcp)
                if not os.path.exists(output_dir):
                    try:
                        os.mkdir(output_dir)
                    except OSError:
                        os.makedirs(output_dir, exist_ok=True)
                    except OSError:
                        print ("Creation of the output directory %s failed" % output_dir)
                        out_create = False
                outdir_writable = os.access(output_dir, os.W_OK)
                if outdir_writable==False:
                    print("WARNING: Output GCP file path set using argument '-m' is not writable. Skipping GCP file creation.")
                    out_create = False

        if out_create == True:
            # subset GDF fields
            GCP_table = GDF[['X_IM', 'Y_IM', 'X_MAP_new', 'Y_MAP_new']]

            # export to text file
            GCP_table.to_csv(options.output_gcp, header=False, index=False, sep=" ")

        # # alternative file export
        # # create GCP string
        # GCP_string = GCP_table.to_string(index=False, header=False, float_format='{:.6f}'.format)
        # # write string to file
        # out_gcp_file = open(out_gcp, 'w')
        # n = out_gcp_file.write(GCP_string)
        # out_gcp_file.close()

    # export resampled GCP (factor 2) to text file
    if options.output_gcp_resampled != None:

        out_create = True
        if os.path.exists(options.output_gcp_resampled):
            if options.overwrite == False:
                print("WARNING: Output resampled GCP file set using argument '-n' already exists. Skipping GCP file creation.")
                out_create = False
            else:
                print("WARNING: Output resampled GCP file set using argument '-n' already exists. Overwriting file.")

            if out_create == True:
                output_dir = os.path.dirname(options.output_gcp_resampled)
                if not os.path.exists(output_dir):
                    try:
                        os.mkdir(output_dir)
                    except OSError:
                        os.makedirs(output_dir, exist_ok=True)
                    except OSError:
                        print ("Creation of the output directory %s failed" % output_dir)
                        out_create = False
                outdir_writable = os.access(output_dir, os.W_OK)
                if outdir_writable==False:
                    print("WARNING: Output resampled GCP file path set using argument '-n' is not writable. Skipping GCP file creation.")
                    out_create = False

        if out_create == True:
            # subset GDF fields
            GCP_table_resampled = GDF[['X_IM', 'Y_IM', 'X_MAP_new', 'Y_MAP_new']]
            GCP_table_resampled['X_IM'] = GCP_table_resampled['X_IM'] / 2
            GCP_table_resampled['Y_IM'] = GCP_table_resampled['Y_IM'] / 2
            GCP_table_resampled['X_IM'] = GCP_table_resampled['X_IM'].astype('int')
            GCP_table_resampled['Y_IM'] = GCP_table_resampled['Y_IM'].astype('int')

            # export to text file
            GCP_table_resampled.to_csv(options.output_gcp_resampled, header=False, index=False, sep=" ")
    # (ema) export coregistration points reliability
    CRL.view_CoRegPoints(figsize = (15,15), attribute2plot = "ABS_SHIFT", exclude_fillVals=False,hide_filtered=False,
    savefigPath = os.path.join(os.path.dirname(options.output_cor_info),"fig.png"))

    CRL.view_CoRegPoints(figsize = (15,15), shapes2plot = "points", backgroundIm = "tgt",
    savefigPath = os.path.join(os.path.dirname(options.output_cor_info),"fig2.png"))

    # export coregistration information to json file
    if options.output_cor_info != None:

        out_create = True
        if os.path.exists(options.output_cor_info):
            if options.overwrite == False:
                print("WARNING: Output coregistration info file set using argument '-l' already exists. Skipping coregistration info file creation.")
                out_create = False
            else:
                print("WARNING: Output coregistration info file set using argument '-l' already exists. Overwriting file.")

            if out_create == True:
                output_dir = os.path.dirname(options.output_cor_info)
                if not os.path.exists(output_dir):
                    try:
                        os.mkdir(output_dir)
                    except OSError:
                        os.makedirs(output_dir, exist_ok=True)
                    except OSError:
                        print ("Creation of the output directory %s failed" % output_dir)
                        out_create = False
                outdir_writable = os.access(output_dir, os.W_OK)
                if outdir_writable==False:
                    print("WARNING: Output coregistration info file path set using argument '-l' is not writable. Skipping coregistration info file creation.")
                    out_create = False

        if out_create == True:
            cor_info_handler = open(options.output_cor_info, 'w')
            json.dump(cor_base_info, cor_info_handler)
            cor_info_handler.close()

    # export tie-point grid vector data to GPKG file
    if options.output_tpg != None:
        out_create = True
        if os.path.exists(options.output_tpg):
            if options.overwrite == False:
                print("WARNING: Output tie-point grid vector file set using argument '-g' already exists. Skipping tie-point grid vector file creation.")
                out_create = False
            else:
                print("WARNING: Output tie-point grid vector file set using argument '-g' already exists. Overwriting file.")

            if out_create == True:
                output_dir = os.path.dirname(options.output_tpg)
                if not os.path.exists(output_dir):
                    try:
                        os.mkdir(output_dir)
                    except OSError:
                        os.makedirs(output_dir, exist_ok=True)
                    except OSError:
                        print ("Creation of the output directory %s failed" % output_dir)
                        out_create = False
                outdir_writable = os.access(output_dir, os.W_OK)
                if outdir_writable==False:
                    print("WARNING: Output tie-point grid vector file path set using argument '-g' is not writable. Skipping tie-point grid vector file creation.")
                    out_create = False

        # export fail because some fields are of type 'bool'
        # solution 1: change field data type
        # GDF['REF_BADDATA'] = GDF['REF_BADDATA'].astype('str')
        # GDF['TGT_BADDATA'] = GDF['TGT_BADDATA'].astype('str')
        # GDF['SSIM_IMPROVED'] = GDF['SSIM_IMPROVED'].astype('str')
        # GDF['L1_OUTLIER'] = GDF['L1_OUTLIER'].astype('str')
        # GDF['L2_OUTLIER'] = GDF['L2_OUTLIER'].astype('str')
        # GDF['L3_OUTLIER'] = GDF['L3_OUTLIER'].astype('str')
        # GDF['OUTLIER'] = GDF['OUTLIER'].astype('str')
        # solution 2: remove boolean fields
        GDF = GDF.drop(columns=['REF_BADDATA','TGT_BADDATA','SSIM_IMPROVED','L1_OUTLIER','L2_OUTLIER','L3_OUTLIER','OUTLIER'])
        GDF = GDF.drop(columns=['GCP'])
        # export data to vector file
        GDF.to_file(driver = 'GPKG', filename= options.output_tpg)

        # import vector data as GeoDataFrame
        # GDF = geopandas.read_file(vector_file)

# print(cor_info)

###########################################################################
# coregister image

if coregister_image == True:

    # get NoData value from input raster
    img_file = gdal.Open(options.input_file)
    img_nodata_value = img_file.GetRasterBand(1).GetNoDataValue()

    # set coregistration parameters
    kwargs = {
        'path_out'     : options.output_file,
        'resamp_alg'   : 'nearest',#ERA CUBIC
        'fmt_out'      : 'GTIFF',
        #'out_crea_options' : 'COMPRESS=LZW',
        'nodata'       : img_nodata_value,
        'align_grids'  : True,
        'match_gsd'    : True,
        'CPUs'         : options.cores,
        'q'            : quiet,
        'v'            : verbose,
    }

    # local image coregistration
    DESHIFTER(options.input_file, cor_info, **kwargs).correct_shifts()
    # DESHIFTER(im2shift=options.input_file, coreg_results=cor_info, path_out=options.output_file, resamp_alg='bilinear', fmt_out='GTIFF', out_crea_options="COMPRESS=LZW", align_grids=True, match_gsd=True, q=False, v=True)

sys.exit(0)

