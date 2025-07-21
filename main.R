#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI

#citations: 

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics
#https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busatto @ CNR-IREA for prismaread package 
#https://github.com/IREA-CNR-MI/prismaread
#distributed under GPL-3.0 license

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL
#https://github.com/GFZ/arosics
#distributed under Apache-2.0 license
#https://gdal.org/en/stable
#distributed under MIT license

#Lorenzo Parigi @ CNR-IREA for smoothing procedure

#Riccardo Canazza for advice in regrid procedure

#Yulun Wu @ University of Ottawa for PRISMA_angle.py code
#https://github.com/yulunwu8/tmart/blob/main/tmart/AEC/read_PRISMA_vaa.py
#distributed under GPL-3.0 license

#_____________________________________________________________________
#inputs ----
#_____________________________________________________________________
#for normal users: modify only things in this section 

#_____________________________________________________________________
#setup folders ----
#_____________________________________________________________________
out_folder  <-  "/space/put_PRISMA_he5_and_S2_tif_here/"
in_file <- base::list.files(path = out_folder, pattern = "\\.he5$", ignore.case = T, full.names = T)
s2_file <- base::list.files(path = out_folder, pattern = "\\.tif$", ignore.case = T, full.names = T)
base::setwd(base::dirname(out_folder))
product_type <- base::substring(base::basename(in_file),5,6)
coreg_out_folder <- base::paste0(out_folder,"1_coreg/")
regrid_out_folder <- base::paste0(out_folder,"2_regrid/")
smoothing_out_folder <- base::paste0(out_folder,"3_smoothing/")

#_____________________________________________________________________
#prismaread ----
#_____________________________________________________________________
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
  in_file = in_file,
  out_folder = out_folder,
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

#_____________________________________________________________________
#write ATCOR parameters ----
#_____________________________________________________________________
if(product_type == "L1" & F){
  prismaread_angle_file_path <- base::paste0(out_folder,gsub(".he5","_HCO.ang",basename(in_file)))
  prismaread_angle_file <- utils::read.table(prismaread_angle_file_path, header =T)
  PRISMA_angle_command <- base::paste0("python"," ",getwd(),"/PRISMA_angle.py")
  sensor_angles <- base::system(PRISMA_angle_command, intern = TRUE) #per questo ci serve h5py in python
  prismaread_angle_file$sensor_zenith <- sensor_angles[1]
  prismaread_angle_file$sensor_azimuth <- sensor_angles[2]
}

#_____________________________________________________________________
#preparing ----
#_____________________________________________________________________
#work on cloud mask for pixel strips
if(product_type == "L1"){
  cloud <- terra::rast(base::paste0(out_folder,base::gsub(".he5","_HCO_CLD.tif",base::basename(in_file))))
  cloud <- terra::subst(cloud, NA, 1)
  terra::plot(cloud, range = c(0,1))
  
  full <- terra::rast(base::paste0(out_folder,base::gsub(".he5","_HCO_FULL.tif",base::basename(in_file))))
  
  cloud_dil <- spatialist::erodil_raster(raster = cloud, #devo aggiungere spatialist al Rstudio server
                                         width = c(3,3), 
                                         type = "box", 
                                         erosion = T, 
                                         dilation = T,
                                         erosion_first = T, 
                                         nt = 1)
  
  terra::plot(cloud_dil)
  
  terra::add(full) <- cloud_dil
  
  terra::writeRaster(full, base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_CLD.tif",base::basename(in_file))),
                     overwrite = T)
  
  prisma_input <- base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_CLD.tif",base::basename(in_file)))
}
if(product_type == "L2"){
  prisma_input <- base::paste0(out_folder,base::gsub(".he5","_HCO_FULL.tif",base::basename(in_file)))
}

#create single layer image to coregister and change crs to EPSG:32632
prisma_projected <- terra::project(x = terra::rast(prisma_input),
                                       y = "epsg:32632",
                                       method = "near")
terra::writeRaster(prisma_projected, 
                   base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_proj.tif",base::basename(in_file))),
                   overwrite = T)
prisma_b52 <- terra::subset(x = prisma_projected, subset = 52)
terra::writeRaster(prisma_b52, 
                   base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_proj_52.tif",base::basename(in_file))),
                   overwrite = T)

#_____________________________________________________________________
#coregistration ----
#_____________________________________________________________________
# set arguments
single_band_reference_image <- s2_file
single_band_image_to_coregister <- base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_proj_52.tif",base::basename(in_file)))
multiband_image_to_coregister <- base::paste0(out_folder,base::gsub(".he5","_HCO_FULL_proj.tif",base::basename(in_file)))
arosics_local_path <- base::paste0("python"," ",base::getwd(),"/arosics_local.py")
target_epsg <- "'EPSG:32632'"
output_directory <- coreg_out_folder
output_file <- base::paste0(coreg_out_folder,"prs_crs_translate_warp.tif")

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
target_epsg <- "'EPSG:32632'"
GDAL_WARP_run_command <- base::paste("gdalwarp -q -of GTiff -r near -order 2 -tap -tr 30 30 -t_srs", target_epsg, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=TRUE), output_file, sep=" ")
base::system(GDAL_WARP_run_command)

#_____________________________________________________________________
#regrid ----
#_____________________________________________________________________
#IMPROVE: should I use gdalUtils::gdalwarp or sen2r::gdalwarp_grid??

base::dir.create(regrid_out_folder, recursive = T, showWarnings = F)

if(product_type == "L1"){
  resample_type <- "near"
}
if(product_type == "L2"){
  resample_type <- "cubicspline"
}

#NOTE: master image from \\10.0.1.243\projects\2022_ASI-PRIS4VEG\3-DATA\images\PRISMA_img_master\PRS_L2D_STD_20200407_HCO_JDS_EXT_FULL_30m_smooth_v1_170b
master_image_path <- base::paste0(base::getwd(),"/regrid_master_image_52.tif")

coreg_out <- base::paste0(coreg_out_folder,"prs_crs_translate_warp.tif")

terra::extend(x = terra::rast(master_image_path),
              y = terra::rast(coreg_out),
              fill = 0,
              filename = base::paste0(regrid_out_folder,"PRISMA_extend.tif"),
              overwrite = T)

terra::resample(x = terra::rast(coreg_out),
                y = terra::rast(base::paste0(regrid_out_folder,"PRISMA_extend.tif")),
                method = resample_type,
                threads = T,
                by_util = T,
                filename = base::paste0(regrid_out_folder,"PRISMA_resample.tif"),
                overwrite = T)

#_____________________________________________________________________
#atcor ----
#_____________________________________________________________________
if(product_type == "L1"){
  print("here we need to put atcor")
}

#_____________________________________________________________________
#smoothing ----
#_____________________________________________________________________

if(product_type == "L2"){
  library(tidytable)
  
  base::dir.create(smoothing_out_folder, recursive = T, showWarnings = F)
  
  input_image_path <- base::paste0(regrid_out_folder,"PRISMA_resample.tif")
  
  smoothing_out <-  base::paste0(smoothing_out_folder, "PRISMA_smoothed.tif")
  
  PRISMA_config <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_spectral_configuration.csv")) %>%
    tidytable::mutate(band_row = tidytable::row_number()) 
  
  PRISMA_bad_bands_table <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_band_selections.csv")) %>%
    tidytable::filter(BB_SUPER_V3 == 1)
  
  input_bad_bands <- PRISMA_bad_bands_table$band
  
  input_wvl <- PRISMA_config$center
  
  output_wvl <- PRISMA_config %>%
    tidytable::filter(BND_SEL  == 1) %>%
    tidytable::pull(center)
  
  
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
  
  
  terra_image <- terra::rast(input_image_path)
  
  if(product_type == "L1"){
    terra_image_sub <- terra::subset(terra_image,subset = 231,negate=T)
  }
  if(product_type == "L2"){
    terra_image_sub <- terra_image
  }
  
  output_image <- terra::app(
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
if(product_type == "L1"){
  print("come here when you are done with ATCOR")
}




