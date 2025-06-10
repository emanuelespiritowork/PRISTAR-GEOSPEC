#main

#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI

in_file  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/PRS_L1_STD_OFFL_20250424100426_20250424100430_0001.he5"
out_folder  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/"
s2_file <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/prove_per_pacchetto/S2_20250422T101051_B08_T32TQQ_ritagliato_coordinate.tif"

product_type <- substring(basename(in_file),5,6)

######################################################################
#prismaread ----
######################################################################
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
#python settings ----
######################################################################

reticulate::py_require(python_version = "3.8.20",
                       packages = c("geos>=3.7.2","gdal==3.6.1","arosics","rasterio"),
                       action = "set")

reticulate::py_config()

######################################################################
#python settings ----
######################################################################
#library(reticulate)
library(here)

reticulate::install_miniconda()#install miniconda

reticulate::conda_binary()#to check conda path

reticulate::conda_create(envname = "C:/prova/Rconda",
                         packages = c("gdal==3.6.1","arosics==1.10.2","rasterio==1.3.4"),
                         python_version = "3.8.20")

reticulate::use_condaenv(condaenv = "C:/prova/Rconda", required = F)

######################################################################
#change CRS ----
######################################################################
if(product_type == "L1"){
  reticulate::source_python(here::here("prs_cld_crs_L1.py"))
}
if(product_type == "L2"){
  reticulate::source_python(here::here("prs_cld_crs_L1.py"))
}

prs_cld_crs(paste0(out_folder,gsub(".he5","_HCO_FULL.tif",basename(in_file))),
            s2_file)





