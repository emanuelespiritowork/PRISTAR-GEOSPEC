library(remotes)
library(rgdal)
library(prismaread)

in_file  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/PRS_L1_STD_OFFL_20250424/PRS_L1_STD_OFFL_20250424100426_20250424100430_0001.he5"
out_folder  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2025/L1/PRS_L1_STD_OFFL_20250424/"

pr_convert(
  in_file = in_file,
  out_folder = out_folder,
  out_format = "GTiff",
  base_georef = T,
  VNIR = T,
  SWIR = T,
  FULL = T,
  ANGLES = T,
  fill_gaps = TRUE,
  source = "HCO",
  join_priority = "SWIR",
  LATLON = TRUE,
  PAN = T,
  CLOUD = T,
  overwrite = T
)

