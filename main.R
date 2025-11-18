#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI

#citations: 

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics and for reading L0 products #https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busetto @ CNR-IREA for prismaread package #https://github.com/IREA-CNR-MI/prismaread #distributed under GPL-3.0 license

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL, for the maintanance of the Docker Container and any hardware-related solution #https://github.com/GFZ/arosics #distributed under Apache-2.0 license #https://gdal.org/en/stable #distributed under MIT license

#Lorenzo Parigi @ CNR-IREA for smoothing procedure

#Riccardo Canazza for advice in regrid procedure

#Yulun Wu @ University of Ottawa for PRISMA_angle.py code #https://github.com/yulunwu8/tmart/blob/main/tmart/AEC/read_PRISMA_vaa.py #distributed under GPL-3.0 license

#_____________________________________________________________________
#inputs ----
#_____________________________________________________________________
#for normal users: modify only things in this section
regrid_option <- "N" #can be N for near, B for bilinear, C for cubic
full_230_bands <- T
validation_for_coreg <- F
PRS_band_for_coreg <- 52

#for expert users:
#procedure_order <- c("inject","read","cloud","coreg","atcor","regrid","crop","smooth")
procedure_order <- c("coreg")
#elements: inject, read, atcor, cloud, coreg, regrid, crop, smooth, ortho

#_____________________________________________________________________
#setup folders ----
#_____________________________________________________________________
source("/space/functions.R")
out_folder  <-  "/space/put_PRISMA_he5_and_S2_tif_here/"
s2_folder <- out_folder
he5_folder <- out_folder
base::setwd(base::dirname(out_folder))

he5_file <- base::list.files(path = out_folder, pattern = "^PRS.*\\.he5$", ignore.case = T, full.names = T)
if(!identical(he5_file,character(0)))
{
  if(grepl("injected",base::basename(he5_file))){
    product_type <- "L0"
  }else{
    product_type <- base::substring(base::basename(he5_file),5,6)
  }
}else{
  product_type <- base::substring(base::basename(he5_file),5,6)
}

#_____________________________________________________________________
#work ----
#_____________________________________________________________________
number_of_operations <- length(procedure_order)
name_of_current_output_folder <- ""
#current_operation <- procedure_order[1]
for(index_of_operations in 1:number_of_operations){
  current_operation <- procedure_order[index_of_operations]
  
  #procedures that are the start of a sequence
  if(current_operation == "inject"){
    print("INJECT")
    if(identical(he5_file,character(0))){
      print("ERRORE")
    }else{
      
      #copy previous he5 file
      dir.create(paste0(out_folder,"/original_L1/"),recursive = T, showWarnings = F)
      invisible(file.copy(he5_file,gsub(out_folder,paste0(out_folder,"/original_L1/"),he5_file)))
      
      injection_command <- base::paste0("python"," ",getwd(),"/Injection_L0_in_L1_cubes.py")
      get_result <- base::system(injection_command, intern = TRUE)
      exit_result <- get_result[length(get_result)]
    }
    
    if(exit_result == "Correctly_finished"){
      he5_file <- base::list.files(path = out_folder, pattern = "^PRS.*\\.he5$", ignore.case = T, full.names = T)
      he5_file_injected <- gsub(".he5$","_injected.he5",he5_file)
      file.rename(he5_file, he5_file_injected)
      he5_file <- he5_file_injected
      product_type <- "L0"
    }else{
      print("ERRORE")
    }
    
  }
  if(current_operation == "read"){
    print("READ")
    if(identical(he5_file,character(0))){
      print("ERRORE")
    }else{
      prismaread_function(product_type, he5_file)
    }
  }
  if(current_operation == "cloud"){
    print("CLOUD")
    cloud_path <- base::list.files(path = he5_folder, pattern = "\\_HCO_CLD.tif$", full.names = T)
    full_path <- base::list.files(path = he5_folder, pattern = "\\_HCO_FULL.tif$", full.names = T)
    if(identical(cloud_path,character(0)) | identical(full_path,character(0))){
      print("ERRORE")
    }else{
      cloud_mask(cloud_path, full_path)
    }
  }
  if(current_operation == "atcor"){
    print("ATCOR")
    angle_file_path <- base::list.files(path = he5_folder, pattern = "\\HCO.ang$", full.names = T)
    if(identical(he5_file,character(0)) | identical(angle_file_path,character(0))){
      print("ERRORE")
    }else{
      atcor_parameters(angle_file_path)
      #when atcor will be implemented in this procedure, here you will find and atcor will be in 
      #the next operations category 
    }
  }
  
  #procedures that can be swapt
  if(current_operation == "coreg" | current_operation == "ortho"){
    #chain part
    print("COREG/ORTHO")
    vector_chain <- check_folder_chain(name_of_current_output_folder, out_folder, current_operation)
    name_of_current_output_folder <- vector_chain[1]
    out_folder <- vector_chain[2]
    coreg_input_path <- check_file_chain(out_folder)
    
    #S2 part
    s2_file <- base::list.files(path = s2_folder, pattern = glob2rx("S2*.tif$"), ignore.case = T, full.names = T)
    
    #DEM part
    base::dir.create(base::paste0(base::getwd(),"/DEM/"), recursive = T, showWarnings = F)
    dem <- F
    dem_path <- NULL
    target_epsg <- paste0("epsg:",terra::crs(terra::rast(s2_file),T,T,T)[3]$code)
    if(!identical(list.files(base::paste0(base::getwd(),"/DEM/")),character(0))){
      dem <- T
      if(!identical(list.files(base::paste0(base::getwd(),"/DEM/"), pattern = "\\_projected.tif$", ignore.case = T, full.names = T), character(0))){
        dem_path <- list.files(base::paste0(base::getwd(),"/DEM/"), pattern = "\\_projected.tif$", ignore.case = T, full.names = T)
      }else{
        dem_path_not_projected <- list.files(base::paste0(base::getwd(),"/DEM/"), pattern = "\\.tif$", ignore.case = T, full.names = T)
        dem_projected <- terra::project(x = terra::rast(dem_path_not_projected),
                                        y = target_epsg,
                                        method = "near")
        file.remove(dem_path_not_projected)
        dem_path <- gsub("\\.tif","\\_projected.tif",dem_path_not_projected)
        terra::writeRaster(dem_projected, dem_path)
      }
      
    }
    
    #if I do not want to use dem, do not use it!
    if(current_operation == "coreg"){
      dem <- F
    }
    
    #coregistration/orthoprojection part
    
    coregistration_to_s2(s2_file,coreg_input_path,name_of_current_output_folder,dem,dem_path,product_type,PRS_band_for_coreg)
    
    
    if(validation_for_coreg){
      base::dir.create(paste0(name_of_current_output_folder,"/validation"), recursive = T, showWarnings = F)
      dem <- F
      output_file <- paste0(name_of_current_output_folder,"/validation/prs_crs_translate_warp.tif")
      if(current_operation == "coreg"){
        input_file <- paste0(name_of_current_output_folder,"/prs_crs_translate_warp.tif")
      }else{
        input_file <- paste0(name_of_current_output_folder,"/raster_focal.tif")
      }
      coregistration_to_s2(s2_file,input_file,paste0(name_of_current_output_folder,"/validation"),dem,dem_path,product_type,PRS_band_for_coreg)
      file.remove(output_file)
    }
    
    
  }
  if(current_operation == "regrid"){
    #chain part
    print("REGRID")
    vector_chain <- check_folder_chain(name_of_current_output_folder, out_folder, current_operation)
    name_of_current_output_folder <- vector_chain[1]
    out_folder <- vector_chain[2]
    regrid_input_path <- check_file_chain(out_folder)
    
    #regrid part
    
    if(regrid_option == "C"){
      resample_type <- "cubicspline"
    }else if(regrid_option == "B"){
      resample_type <- "bilinear"
    }else{
      resample_type <- "near"
    }
    
    #NOTE: for Jolanda di Savoia use master image 
    #from //10.0.1.243/projects/2022_ASI-PRIS4VEG/3-DATA/images/PRISMA_img_master/PRS_L2D_STD_20200407_HCO_JDS_EXT_FULL_30m_smooth_v1_170b
    #for Piacenza the master image is not available yet
    
    master_image_path <- base::list.files(base::paste0(base::getwd(),"/master_image_for_regridding/"), full.names = T, pattern = "\\.tif$")
    
    #regrid_input_path <- base::list.files(out_folder, full.names = T, pattern = "\\.bsq$")
    
    regrid_function(master_image_path, name_of_current_output_folder, regrid_input_path, resample_type)
  }
  
  if(current_operation == "crop"){
    #chain part
    print("CROP")
    vector_chain <- check_folder_chain(name_of_current_output_folder, out_folder, current_operation)
    name_of_current_output_folder <- vector_chain[1]
    out_folder <- vector_chain[2]
    crop_input_path <- check_file_chain(out_folder)
    
    #crop part
    
    #NOTE: for Jolanda di Savoia use master image 
    #from //10.0.1.243/projects/2022_ASI-PRIS4VEG/3-DATA/images/PRISMA_img_master/PRS_L2D_STD_20200407_HCO_JDS_EXT_FULL_30m_smooth_v1_170b
    #for Piacenza the master image is not available yet
    
    master_image_path <- base::list.files(base::paste0(base::getwd(),"/master_image_for_regridding/"), full.names = T, pattern = "\\.tif$")
    
    #regrid_input_path <- base::list.files(out_folder, full.names = T, pattern = "\\.bsq$")
    
    crop_function(master_image_path, name_of_current_output_folder, crop_input_path)
  }
  if(current_operation == "smooth"){
    #chain part
    print("SMOOTH")
    vector_chain <- check_folder_chain(name_of_current_output_folder, out_folder, current_operation)
    name_of_current_output_folder <- vector_chain[1]
    out_folder <- vector_chain[2]
    terra_image_path <- check_file_chain(out_folder)
    
    #smooth part
    smoothing_out <-  base::paste0(name_of_current_output_folder, "/PRISMA_smoothed.tif")
    
    if("cloud" %in% procedure_order | !identical(list.files(he5_folder,pattern = "*_FULL_CLD.tif$"),character(0))){
      cloud_smooth <- T
    }else{
      cloud_smooth <- F
    }
    
    
    
    library(tidytable)
    
    PRISMA_config <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_spectral_configuration.csv")) %>%
      tidytable::mutate(band_row = tidytable::row_number()) 
    
    PRISMA_bad_bands_table <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_band_selections.csv")) %>%
      tidytable::filter(BB_SUPER_V3 == 1)
    
    #this would be the idea but when I put the smoothing code into a new function the terra::app does not work
    #smooth_spectra(terra_image_path,PRISMA_config,PRISMA_bad_bands_table,smoothing_out,cloud_smooth)
    
    input_bad_bands <- PRISMA_bad_bands_table$band
    
    input_wvl <- PRISMA_config$center
    
    selection_vector <- 1
    
    if(full_230_bands){
      selection_vector <- c(0,1)
    }
    
    #which output bands
    output_wvl <- PRISMA_config %>%
      tidytable::filter(BND_SEL %in% selection_vector) %>%
      tidytable::pull(center)
    
    #print("Define spline function")
    
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
    
    #print("Read terra image")
    
    terra_image <- terra::rast(terra_image_path)
    
    if(cloud_smooth){
      terra_image_sub <- terra::subset(terra_image,subset = 231,negate=T)
    }else{
      terra_image_sub <- terra_image
    }
    
    #print("Apply smoothing")
    
    #terra::terraOptions(memmin = 30, print=T, progress = 1, memfrac = 0.8, verbose = T)
    
    #terra::gdalCache(1000000)
    
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
  
}






