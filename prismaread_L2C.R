library(remotes)
library(rgdal)
library(prismaread)

in_file  <-  "C:/Users/emast/Desktop/250606_agile/20230304/PRS_L2C_STD_20230304102047_20230304102051_0001.he5"
out_folder  <-  "C:/Users/emast/Desktop/250606_agile/20230304/"

pr_convert(
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
  overwrite = T,
  ATCOR = F
)
