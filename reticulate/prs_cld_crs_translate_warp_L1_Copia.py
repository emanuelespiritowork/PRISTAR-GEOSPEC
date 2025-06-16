def prs_cld_crs_translate_warp(path,im_target,im_reference):
    from arosics import COREG
    from arosics import COREG_LOCAL
    import matplotlib.pyplot as plt
    import numpy as np
    from osgeo import gdal, osr, ogr
    import pandas as pd
    import os
    import geopandas as gpd
    CRL = COREG_LOCAL(im_reference,im_target,grid_res = 10,window_size =(128,128),max_shift = 10,resamp_alg_deshift = 'nearest',resamp_alg_calc = 'nearest',path_out = path+'Registered_AROSICS_Local_B8_B52.tif',fmt_out = 'GTIFF',max_iter = 8,s_b4match = 51,q=True)
    CRL.correct_shifts()
    CRL.view_CoRegPoints(figsize=(15,15))
    points = CRL.CoRegPoints_table
    points
    CRL.tiepoint_grid.to_PointShapefile(path_out=path+'AROSICS_TiePoints_Local_B8_B52.shp')
    lista_files = []
    dirFileList = os.listdir(path) #legge la cartella 1
    os.chdir(path) #setta la cartella
    for file in dirFileList: #per ciascun file nella cartella 1
        if os.path.splitext(file)[-1] == '.shp': #le cui ultime lettere sono ".img"
            lista_files.append(os.path.join(path, file)) #aggiungile alla lista vuot
    lista_files
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
    #print columns
    np.unique(f.columns)
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
    gcps = pd.DataFrame({'X': X_MAP[:,]+X_SHIFT_M[:,],
                         'Y': Y_MAP[:,]+Y_SHIFT_M[:,],
                         'Col': X_IM[:,],
                         'Row': Y_IM[:,]})
    gcps_gdal = [gdal.GCP(row['X'], row['Y'],0, row['Col'], row['Row']) for index, row in gcps.iterrows()] ### "gcps" or "gcps_PAN"
    kwargs = {
        'format': 'GTiff',
        'outputType': gdal.GDT_UInt16}
    output_image = path+'prs_cld_crs_translate.tif'
    ds_gcp = gdal.Translate(output_image,
                            path+'prs_cld_crs.tif',
                            outputSRS='EPSG:32632',
                            GCPs=gcps_gdal,
                            **kwargs)
    options = gdal.WarpOptions(dstSRS='EPSG:32632', polynomialOrder=2, targetAlignedPixels=True, xRes=30, yRes =30)
    ds = gdal.Warp(path+'prs_cld_crs_translate_warp.tif', ds_gcp, dstNodata = np.nan, options=options, resampleAlg="near")
    return 0




