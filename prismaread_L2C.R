library(remotes)
library(rgdal)
library(prismaread)

in_file  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2023/L2C/PRS_L2C_STD_20230407/PRS_L2C_STD_20230407100729_20230407100733_0001.he5"
out_folder  <-  "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2023/L2C/PRS_L2C_STD_20230407/"

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
