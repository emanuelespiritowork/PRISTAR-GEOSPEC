
#_____________________________________________________________________
#prismaread ----
#_____________________________________________________________________
prismaread_function <- function(product_type, he5_file){
  
  if(product_type == "L1"){
    CLOUD <- T
    ATCOR <- T
  }
  if(product_type == "L2"){
    CLOUD <- F
    ATCOR <- F
  }
  
  prismaread::pr_convert(
    ANGLES = T,
    in_file = he5_file,
    out_folder = paste0(base::dirname(he5_file),"/"),
    out_format = "GTiff",
    base_georef = T,
    VNIR = F,
    SWIR = F,
    FULL = T,
    fill_gaps = T,
    source = "HCO",
    join_priority = "SWIR",
    LATLON = T,
    PAN = F,
    CLOUD = CLOUD,
    overwrite = T,
    ATCOR = ATCOR
  )
}

#_____________________________________________________________________
#atcor and prisma_angle ----
#_____________________________________________________________________
atcor_parameters <- function(angle_file_path){
  prismaread_angle_file <- utils::read.table(angle_file_path, header =T)
  PRISMA_angle_command <- base::paste0("python"," ",getwd(),"/PRISMA_angle.py")
  sensor_angles <- base::system(PRISMA_angle_command, intern = TRUE)
  
  sensor_zenith <- strsplit(sensor_angles,",")[[1]][1]
  sensor_zenith <- gsub('^.', '', sensor_zenith)
  
  sensor_azimuth <- strsplit(sensor_angles,",")[[1]][2]
  sensor_azimuth <- gsub('^.|.$', '', sensor_azimuth)
  
  prismaread_angle_file$sensor_zenith <- as.numeric(sensor_zenith)
  prismaread_angle_file$sensor_azimuth <- as.numeric(sensor_azimuth)
  
  write.csv(prismaread_angle_file,paste0(base::dirname(angle_file_path),"/ATCOR/all_angles_file.csv"), quote = F, row.names = F)
}

cloud_mask <- function(cloud_path, full_path){
  cloud <- terra::rast(cloud_path)
  cloud <- terra::subst(cloud, NA, 1)
  terra::plot(cloud, range = c(0,1))
  
  cloud_dil <- spatialist::erodil_raster(raster = cloud, 
                                         width = c(3,3), 
                                         type = "box", 
                                         erosion = T, 
                                         dilation = T,
                                         erosion_first = T, 
                                         nt = 1)
  
  terra::plot(cloud_dil)
  
  full <- terra::rast(full_path)
  
  terra::add(full) <- cloud_dil
  
  terra::writeRaster(full, 
                     base::gsub("_HCO_FULL.tif","_HCO_FULL_CLD.tif",full_path),
                     overwrite = T)
}

#_____________________________________________________________________
#coreg ----
#_____________________________________________________________________
coregistration_to_s2 <- function(s2_file,coreg_input_path,coreg_proj_path,coreg_out_folder){
  #create single layer image to coregister and change crs to EPSG:32632
  target_epsg <- paste0("epsg:",terra::crs(terra::rast(s2_file),T,T,T)[3]$code)
  
  prisma_projected <- terra::project(x = terra::rast(coreg_input_path),
                                     y = target_epsg,
                                     method = "near")
  terra::writeRaster(prisma_projected, 
                     coreg_proj_path,
                     overwrite = T)
  prisma_b52 <- terra::subset(x = prisma_projected, subset = 52)
  terra::writeRaster(prisma_b52, 
                     base::gsub("proj","proj_52",coreg_proj_path),
                     overwrite = T)
  
  # set arguments
  single_band_reference_image <- s2_file
  single_band_image_to_coregister <- base::list.files(base::dirname(s2_file),"\\proj_52.tif$", full.names = T)
  multiband_image_to_coregister <- base::list.files(base::dirname(s2_file),"\\proj.tif$", full.names = T)
  arosics_local_path <- base::paste0("python"," ",base::getwd(),"/arosics_local.py")
  output_directory <- coreg_out_folder
  output_file <- base::paste0(coreg_out_folder,"/prs_crs_translate_warp.tif")
  
  base::dir.create(coreg_out_folder, recursive = T, showWarnings = F)
  
  # setup AROSICS run command
  #IMPROVE: potrei usare AROSICS da riga di comando?
  arosics_run_command <- base::paste(arosics_local_path, "-v -r", single_band_reference_image, "-t", single_band_image_to_coregister, "-l", base::normalizePath(path=base::paste(output_directory, "/", "AROSICS_coregistration_info.json", sep=""), winslash="/", mustWork=FALSE), "-m", base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE), "-g", base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE), sep=" ")
  # run AROSICS to get GCPs
  base::system(arosics_run_command)
  # import GCPs
  gcp_lines <- base::readLines(con=base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE))
  # generate GCP line to be used in GDAL
  gcp_args <- base::paste(base::paste("-gcp", gcp_lines, sep=" "), collapse=" ")
  # create VRT with GCP
  GDAL_VRT_run_command <- base::paste("gdal_translate -q -of VRT", gcp_args, multiband_image_to_coregister, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=FALSE), sep=" ")
  base::system(GDAL_VRT_run_command)
  # warp input image using second order polynomial
  GDAL_WARP_run_command <- base::paste("gdalwarp -q -of GTiff -r near -order 2 -tap -tr 30 30 -t_srs", target_epsg, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=TRUE), output_file, sep=" ")
  base::system(GDAL_WARP_run_command)
  
}

#_____________________________________________________________________
#regrid ----
#_____________________________________________________________________
regrid_function <- function(master_image_path, name_of_current_output_folder, regrid_input_path){
  terra::extend(x = terra::rast(master_image_path),
                y = terra::rast(regrid_input_path),
                fill = 0,
                filename = base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"),
                overwrite = T)
  
  terra::resample(x = terra::rast(regrid_input_path),
                  y = terra::rast(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif")),
                  method = resample_type,
                  threads = T,
                  by_util = T,
                  filename = base::paste0(name_of_current_output_folder,"/PRISMA_resample.tif"),
                  overwrite = T)
  
  file.remove(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"))
  
}

#_____________________________________________________________________
#crop ----
#_____________________________________________________________________
crop_function <- function(master_image_path, name_of_current_output_folder, crop_input_path){
  terra::crop(x = terra::rast(crop_input_path),
              y = terra::rast(master_image_path),
              filename = base::paste0(name_of_current_output_folder,"/PRISMA_crop.tif"),
              overwrite = T)
}

#_____________________________________________________________________
#smooth ----
#_____________________________________________________________________
smooth_spectra <- function(terra_image_path,smoothing_out,cloud_smooth){
  PRISMA_config <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_spectral_configuration.csv")) %>%
    tidytable::mutate(band_row = tidytable::row_number()) 
  
  PRISMA_bad_bands_table <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_band_selections.csv")) %>%
    tidytable::filter(BB_SUPER_V3 == 1)
  
  input_bad_bands <- PRISMA_bad_bands_table$band
  
  input_wvl <- PRISMA_config$center
  
  output_wvl <- PRISMA_config %>%
    tidytable::filter(BND_SEL  == 1) %>%
    tidytable::pull(center)
  
  print("Define spline function")
  
  spline_fun <- function(pixel, band_center_input, bad_bands_pos, band_center_output, df = 40) {
    # togliamo le bande cattive
    ref_valid <- pixel[-bad_bands_pos]
    wvl_valid <- band_center_input[-bad_bands_pos]
    
    if (base::any(is.na(ref_valid))) {
      return(base::rep(NA_real_, base::length(band_center_output)))
    }
    
    sp <- stats::smooth.spline(x = wvl_valid, y = ref_valid, df = df)
    y_smooth <- stats::predict(sp, x = band_center_output)$y
    # limitiamo a zero
    y_smooth[y_smooth < 0] <- 0
    return(y_smooth)
  }
  
  print("Read terra image")
  
  terra_image <- terra::rast(terra_image_path)
  
  if(cloud_smooth){
    terra_image_sub <- terra::subset(terra_image,subset = 231,negate=T)
  }else{
    terra_image_sub <- terra_image
  }
  
  print("Apply smoothing")
  
  terra::terraOptions(memmin = 30, print=T, progress = 1, memfrac = 0.8, verbose = T)
  
  terra::gdalCache(1000000)
  
  terra::app(
    x = terra_image_sub,
    fun = spline_fun,
    band_center_input = input_wvl, 
    bad_bands_pos = input_bad_bands, 
    band_center_output = output_wvl,
    df=40,
    
    cores = 7,                     
    filename = smoothing_out,
    overwrite = TRUE,
    wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
}