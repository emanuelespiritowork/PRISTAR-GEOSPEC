#main

#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI
out_folder  <-  "/space/put_PRISMA_he5_and_S2_tif_here/"
in_file <- list.files(path = out_folder, pattern = "\\.he5$", ignore.case = T, full.names = T)
s2_file <- list.files(path = out_folder, pattern = "\\.tif$", ignore.case = T, full.names = T)
#in_file <- "/space/prove_per_pacchetto/PRS_L1_STD_OFFL_20250424100426_20250424100430_0001.he5"
#s2_file <- "/space/prove_per_pacchetto/S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif"
setwd(dirname(out_folder))


#in_file <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/PRS_L1_STD_OFFL_20250424100426_20250424100430_0001.he5"
#out_folder <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/"
#s2_file <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif"

######################################################################
#deduce folders ----
######################################################################
product_type <- base::substring(base::basename(in_file),5,6)
coreg_out_folder <- base::paste0(out_folder,"coreg/")
smoothing_out_folder <- base::paste0(out_folder,"smoothing/")
regrid_out_folder <- base::paste0(out_folder,"regrid/")

######################################################################
#prismaread ----
######################################################################
#https://github.com/IREA-CNR-MI/prismaread

if(product_type == "L1"){
  CLOUD <- T
  ATCOR <- T
}
if(product_type == "L2"){
  CLOUD <- F
  ATCOR <- F
}

prismaread::pr_convert(
  in_file = in_file,
  out_folder = out_folder,
  out_format = "GTiff",
  base_georef = T,
  VNIR = T,
  SWIR = T,
  FULL = T,
  ANGLES = T,
  fill_gaps = T,
  source = "HCO",
  join_priority = "SWIR",
  LATLON = T,
  PAN = T,
  CLOUD = CLOUD,
  overwrite = T,
  ATCOR = ATCOR
)

######################################################################
#coregistration ----
######################################################################

#create single layer image to coregister
prisma_input <- paste0(out_folder,gsub(".he5","_HCO_FULL.tif",basename(in_file)))
prisma_projected <- terra::project(x = terra::rast(prisma_input),
                                       y = "epsg:32632",
                                       method = "near")
terra::writeRaster(prisma_projected, paste0(out_folder,gsub(".he5","_HCO_FULL_proj.tif",basename(in_file))))
prisma_b52 <- terra::subset(x = prisma_projected, subset = 52)
terra::plot(prisma_projected)
terra::plot(prisma_b52)
#change crs to EPSG:32632

terra::writeRaster(prisma_b52, paste0(out_folder,gsub(".he5","_HCO_FULL_proj_52.tif",basename(in_file))))

# written by Federico Filipponi @ CNR-IGAG

# set arguments (questo lo devi modificare con i percorsi definiti nel docker container)
single_band_reference_image <- s2_file
single_band_image_to_coregister <- paste0(out_folder,gsub(".he5","_HCO_FULL_proj_52.tif",basename(in_file)))
multiband_image_to_coregister <- paste0(out_folder,gsub(".he5","_HCO_FULL_proj.tif",basename(in_file)))
arosics_local_path <- paste0("python"," ",getwd(),"/arosics_local.py")
target_epsg <- "'EPSG:32632'"
output_directory <- coreg_out_folder
output_file <- paste0(coreg_out_folder,"prs_crs_translate_warp.tif")

dir.create(coreg_out_folder, recursive = T, showWarnings = F)


# setup AROSICS run command (potrei usare AROSICS da riga di comando?)
arosics_run_command <- paste(arosics_local_path, "-v -r", single_band_reference_image, "-t", single_band_image_to_coregister, "-l", normalizePath(path=paste(output_directory, "/", "AROSICS_coregistration_info.json", sep=""), winslash="/", mustWork=FALSE), "-m", normalizePath(path=paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE), "-g", normalizePath(path=paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE), sep=" ")
# run AROSICS to get GCPs
system(arosics_run_command)
# import GCPs
gcp_lines <- readLines(con=normalizePath(path=paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE))
# generate GCP line to be used in GDAL
gcp_args <- paste(paste("-gcp", gcp_lines, sep=" "), collapse=" ")
# create VRT with GCP
GDAL_VRT_run_command <- paste("gdal_translate -q -of VRT", gcp_args, multiband_image_to_coregister, normalizePath(path=paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=FALSE), sep=" ")
system(GDAL_VRT_run_command)
# warp input image using second order polynomial
target_epsg <- "'EPSG:32632'"
GDAL_WARP_run_command <- paste("gdalwarp -q -of GTiff -r near -order 2 -tap -tr 30 30 -t_srs", target_epsg, normalizePath(path=paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=TRUE), output_file, sep=" ")
system(GDAL_WARP_run_command)


######################################################################
#smoothing ----
######################################################################
#written by Lorenzo Parigi @ CNR-IREA

library(tidytable)

base::dir.create(smoothing_out_folder, recursive = T, showWarnings = F)

input_image_path <- base::paste0(coreg_out_folder,"prs_crs_translate_warp.tif")

smoothing_out <-  paste0(smoothing_out_folder, "PRISMA_smoothed.tif")

PRISMA_config <- tidytable::fread(paste0(getwd(),"/PRISMA_spectral_configuration.csv")) %>%
  tidytable::mutate(band_row = tidytable::row_number()) 

PRISMA_bad_bands_table <- tidytable::fread(paste0(getwd(),"/PRISMA_band_selections.csv")) %>%
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

######################################################################
#regrid ----
######################################################################
#thanks to advice of Riccardo Canazza @ CNR-IREA

#should I use gdalUtils::gdalwarp or sen2r::gdalwarp_grid??

base::dir.create(regrid_out_folder, recursive = T, showWarnings = F)

if(product_type == "L1"){
  resample_type <- "near"
}
if(product_type == "L2"){
  resample_type <- "cubicspline"
}

#master image from \\10.0.1.243\projects\2022_ASI-PRIS4VEG\3-DATA\images\PRISMA_img_master\PRS_L2D_STD_20200407_HCO_JDS_EXT_FULL_30m_smooth_v1_170b
master_image_path <- paste0(getwd(),"/regrid_master_image_52.tif")

terra::extend(x = terra::rast(master_image_path),
              y = terra::rast(smoothing_out),
              fill = 0,
              filename = paste0(regrid_out_folder,"PRISMA_extend.tif"),
              overwrite = T)

terra::resample(x = terra::rast(smoothing_out),
                y = terra::rast(paste0(regrid_out_folder,"PRISMA_extend.tif")),
                method = resample_type,
                threads = T,
                by_util = T,
                filename = paste0(regrid_out_folder,"PRISMA_resample.tif"),
                overwrite = T)
