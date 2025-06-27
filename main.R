#main

#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI

in_file <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2023/L2C/prove_per_pacchetto/PRS_L2C_STD_20230304102047_20230304102051_0001.he5"
out_folder  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2023/L2C/prove_per_pacchetto/"
s2_file <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2023/L2C/prove_per_pacchetto/S2_20230309_B08_T32TQQ_ritagliato_QGIS.tif"
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
#set working directory of python shell ----
######################################################################
command <- paste0("cd", " ", here::here())#this needs the script to be in a mount drived, not in the server
#alternatives at https://superuser.com/questions/1014248/change-working-directory-to-network-share
system(command)

#this way the prompt is set to a specific working directory
#when executing python, its working directory will be the same as this?

######################################################################
#change CRS ----
######################################################################
#first I write the inputs of prs_cld_crs.py into a file
prisma_input <- paste0(out_folder,gsub(".he5","_HCO_FULL.tif",basename(in_file)))
s2_input <- s2_file
prs_cld_crs <- c(prisma_input,s2_input)
utils::write.csv(prs_cld_crs, 
          paste0(here::here(),"/inputs_of_python_code.csv"),
          row.names = F,
          col.names = F, 
          quote = F)


if(product_type == "L1"){
  prs_crs_python_path <- paste0(here::here(),"/prs_cld_crs_L1.py")
}
if(product_type == "L2"){
  prs_crs_python_path <- paste0(here::here(),"/prs_cld_crs_L2C.py")
}
command <- paste0("python ",prs_crs_python_path)
system(command)


######################################################################
#warp ----
######################################################################
#https://doi.org/10.1016/j.isprsjprs.2024.07.003
#https://doi.org/10.5281/zenodo.11547257
if(product_type == "L1"){
  prs_crs_warp_python_path <- paste0(here::here(),"/prs_cld_crs_translate_warp_L1.py")
}
if(product_type == "L2"){
  prs_crs_warp_python_path <- paste0(here::here(),"/prs_cld_crs_translate_warp_L2C.py")
}
command <- paste0("python ",prs_crs_warp_python_path)
system(command)


######################################################################
#smoothing ----
######################################################################
#thanks to Lorenzo Parigi @ CNR-IREA

library(tidytable)

base::dir.create(smoothing_out_folder, recursive = T, showWarnings = F)

if(product_type == "L1"){
  input_image_path <- base::paste0(coreg_out_folder,"prs_cld_crs_translate_warp.tif")
}
if(product_type == "L2"){
  input_image_path <- base::paste0(coreg_out_folder,"prs_crs_translate_warp.tif")
}

smoothing_out <-  paste0(smoothing_out_folder, "PRISMA_smoothed.tif")

PRISMA_config <- tidytable::fread(here::here("PRISMA_spectral_configuration.csv")) %>%
  tidytable::mutate(band_row = tidytable::row_number()) 

PRISMA_bad_bands_table <- tidytable::fread(here::here("PRISMA_band_selections.csv")) %>%
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
#thanks to Riccardo Canazza @ CNR-IREA

#should I use gdalUtils::gdalwarp or sen2r::gdalwarp_grid??

base::dir.create(regrid_out_folder, recursive = T, showWarnings = F)

if(product_type == "L1"){
  resample_type <- "nearest"
}
if(product_type == "L2"){
  resample_type <- "cubicspline"
}

master_image_path <- "//10.0.1.243/projects/2022_ASI-PRIS4VEG/3-DATA/images/PRISMA_img_master/PRS_L2D_STD_20200407_HCO_JDS_EXT_FULL_30m_smooth_v1_170b"

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









######################################################################
#python settings ----
######################################################################

reticulate::conda_create(envname = "C:/prova/Rconda",
                         packages = c("gdal==3.6.1","arosics==1.10.2","rasterio==1.3.4"),                         packages = c("gdal==3.6.1","arosics==1.10.2","rasterio==1.3.4"),
                         python_version = "3.8.20")

Sys.setenv(RETICULATE_PYTHON = "C:/prova/Rconda/python.exe")

#way to force the environment
#usethis::edit_r_environ(scope = "project")

reticulate::py_discover_config(use_environment = "C:/prova/Rconda/python.exe")

reticulate::use_python(python = "C:/prova/Rconda/python.exe")

reticulate::py_config()

#try the environment
#reticulate::repl_python()

