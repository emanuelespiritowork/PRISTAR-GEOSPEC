#designed in Jun 2025
#author: Emanuele Spirito (spirito.e@irea.cnr.it)
#site: CNR-IREA-MI

#citations: 

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics and for reading L0 products #https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busetto @ CNR-IREA for prismaread package #https://github.com/IREA-CNR-MI/prismaread #distributed under GPL-3.0 license

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL, for the maintanance of the Docker Container and any hardware-related solution #https://github.com/GFZ/arosics #distributed under Apache-2.0 license #https://gdal.org/en/stable #distributed under MIT license

#Lorenzo Parigi @ CNR-IREA for smoothing procedure, new Docker structure, API and isofit incorporation 

#Riccardo Canazza for advice  in regrid procedure

#Yulun Wu @ University of Ottawa for PRISMA_angle.py code #https://github.com/yulunwu8/tmart/blob/main/tmart/AEC/read_PRISMA_vaa.py #distributed under GPL-3.0 license

#_____________________________________________________________________
#User inputs ----
#_____________________________________________________________________
#for normal users: modify only things in this section
regrid_option <- "N" #can be N for near, B for bilinear, C for cubic
full_230_bands <- T
PRS_band_for_coreg <- 52
shift <- F
shift_x <- -8000
shift_y <- 0
n_threads <- 7
dem_root_path <- "//10.0.1.243/nr_data/4_rs_product/DTM/Italia/Tinitaly/data"

#for expert users:
#procedure_order <- c("inject","read","cloud","coreg","atcor","regrid","crop","smooth","addmetadata")
# procedure_order <- c("inject","read","coreg")
procedure_order <- c("read","coreg","isofit")
#elements: inject, read, cloud, coreg, regrid, crop, smooth, ortho,addmetadata, isofit

#_____________________________________________________________________
# Main -----
#_____________________________________________________________________
source("/config_folder/functions.R")
#input folder identification
root_folders <- list.dirs(path = "/space",
                          recursive = F)

PRISTAR_processing <- function(root_folder){
  out_folder <-  paste(root_folder, "PRISTAR_processing", sep = "/")
  dir.create(out_folder,
             recursive = F,
             showWarnings = T)
  
  
  
  #1.1 identify file paths ----
  he5_path <- base::list.files(path = root_folder, pattern = "^PRS.*\\.he5$", ignore.case = T, full.names = T)
  s2_path <- base::list.files(path = root_folder, pattern = glob2rx("S2*.tif$"), ignore.case = T, full.names = T)
  dem_path <- base::list.files(path = "/config_folder/DEM", pattern = "dtm", full.names = T)
  master_image_path <- base::list.files("/config_folder/master_image_for_regridding/", full.names = T, pattern = "\\.tif$")
  
  #1.2 identify product_type ----
  product_type <- identify_product_type(he5_path)
  
  #_____________________________________________________________________
  #1.3 create workflow ======
  #_____________________________________________________________________
  
  ##1.3.0 check consistency of chain ----
  if(!check_consistency(procedure_order)){
    stop("Not correct sequence in the procedure_order parameters.")
  }
  
  ##1.3.1 execute unchained operations ----
  select_unchained_operations <- procedure_order[procedure_order %in% c("read","cloud","atcor","inject")]
  select_chained_operations <- procedure_order[!(procedure_order %in% c("read","cloud","atcor","inject"))]
  cloud_present_in_stack <- F
  
  number_of_unchained_operations <- length(select_unchained_operations)
  if(number_of_unchained_operations > 0){
    unchained_out_folder <- paste(out_folder, "0_read", sep = "/")
    dir.create(unchained_out_folder, showWarnings = F, recursive = F)
  }
  for(index_of_unchained_operations in 1:number_of_unchained_operations){
    current_operation <- select_unchained_operations[index_of_unchained_operations]
    
    ###1.3.1.1 "inject" operation ----
    if(current_operation == "inject"){
      print("INJECT")
      if(identical(he5_path,character(0))){
        stop(paste0("No he5 file found in ", root_folder))
      }else{
        
        #copy previous he5 file
        dir.create(paste0(out_folder,"/0_original_L1/"),recursive = F, showWarnings = F)
        invisible(file.copy(from = he5_path,
                            to = paste0(paste0(out_folder,"/0_original_L1/"),basename(he5_path))))
        
        injection_command <- base::paste0("python"," ","/space/Injection_L0_in_L1_cubes.py")
        get_result <- base::system(injection_command, intern = TRUE)
        exit_result <- get_result[length(get_result)]
      }
      
      if(exit_result == "Correctly_finished"){
        he5_file_injected <- gsub(".he5$","_injected.he5",he5_path)
        file.rename(he5_path, he5_file_injected)
        he5_path <- he5_file_injected
        product_type <- "L0"
      }else{
        stop("Error during injection. See error of python console")
      }
      
    }
    
    ##1.3.1.2 "read" operation ----
    if(current_operation == "read"){
      print("READ")
      if(identical(he5_path,character(0))){
        stop(paste0("No he5 file found in ", root_folder))
      }else{
        prismaread_function(product_type = product_type, 
                            he5_path = he5_path, 
                            unchained_out_folder = unchained_out_folder, 
                            root_folder = root_folder)
      }
    }
    
    ##1.3.1.3 "cloud" operation ----
    if(current_operation == "cloud"){
      print("CLOUD")
      
      if(product_type == "L0"){
        stop("Generation of cloud mask not relative to the L0 product but to the originary L1 product")
      }
      
      if(product_type == "L2"){
        stop("No generation of cloud mask for L2 product")
      }
      
      cloud_path <- base::list.files(path = unchained_out_folder, pattern = "\\_HCO_CLD.tif$", full.names = T, recursive = F)
      full_path <- base::list.files(path = unchained_out_folder, pattern = "\\_HCO_FULL.tif$", full.names = T, recursive = F)
      if(identical(cloud_path,character(0)) | identical(full_path,character(0))){
        stop(paste0("Lacking either _FULL.tif or _CLD.tif in the folder ",unchained_out_folder))
      }else{
        cloud_mask(cloud_path = cloud_path, 
                   full_path = full_path)
      }
      
      cloud_present_in_stack <- T
    }
    
  }
  
  ##1.3.2 execute chained operations ----
  PRISMA_wvl_info <- tidytable::fread(list.files(path = unchained_out_folder,
                                                 pattern = "*.wvl$",
                                                 full.names = T))
  
  PRISMA_angle_info <- tidytable::fread(list.files(path = paste0(unchained_out_folder,"/ATCOR/"),
                                                   pattern = "all_angles_file.csv",
                                                   full.names = T))
  
  PRISMA_config <- tidytable::fread(base::paste0("/config_folder/PRISMA_spectral_configuration.csv")) |>
    tidytable::mutate(band_row = tidytable::row_number()) 
  
  number_of_chained_operations <- length(select_chained_operations)
  name_of_current_output_folder <- ""
  
  if(number_of_chained_operations > 0){
    for(index_of_chained_operations in 1:number_of_chained_operations){
      current_operation <- select_chained_operations[index_of_chained_operations]
      
      ###1.3.2.0 check vector_chain ----
      name_of_current_output_folder <- check_folder_chain(name_of_current_output_folder = name_of_current_output_folder, 
                                                          out_folder = out_folder, 
                                                          current_operation = current_operation)
      if(index_of_chained_operations == 1){
        input_file_path <- get_starting_prisma_image(unchained_out_folder = unchained_out_folder,
                                                     cloud_present_in_stack = cloud_present_in_stack)
      }else{
        input_file_path <- output_file_path
      }
      
      output_file_path <- naming_convention(out_folder = out_folder,
                                            input_file_path = input_file_path, 
                                            name_of_current_output_folder = name_of_current_output_folder, 
                                            current_operation = current_operation, 
                                            index_of_chained_operations = index_of_chained_operations,
                                            PRISMA_angle_info = PRISMA_angle_info, 
                                            product_type = product_type)
      
      ##1.3.2.1 "coreg" or "ortho" operation ----
      if(current_operation == "coreg" | current_operation == "ortho"){
        if(current_operation == "coreg"){
          dem <- F
          print("COREG")
        }else{
          dem <- T
          print("ORTHO")
        }
        
        coregistration_to_s2(s2_path = s2_path,
                             input_file_path = input_file_path,
                             output_file_path = output_file_path,
                             dem = dem,
                             dem_path = dem_path,
                             product_type = product_type,
                             PRS_band_for_coreg = PRS_band_for_coreg,
                             shift = shift,
                             shift_x = shift_x,
                             shift_y = shift_y)
        
        # if(validation_for_coreg){
        #   base::dir.create(paste0(name_of_current_output_folder,"/validation"), recursive = T, showWarnings = F)
        #   dem <- F
        #   output_file <- paste0(name_of_current_output_folder,"/validation/prs_crs_translate_warp.tif")
        #   if(current_operation == "coreg"){
        #     input_file <- paste0(name_of_current_output_folder,"/prs_crs_translate_warp.tif")
        #   }else{
        #     input_file <- paste0(name_of_current_output_folder,"/raster_focal.tif")
        #   }
        #   coregistration_to_s2(s2_file,input_file,paste0(name_of_current_output_folder,"/validation"),dem,dem_path,product_type,PRS_band_for_coreg,all_wvl)
        #   file.remove(output_file)
        # }
      }
      
      ##1.3.2.2 "regrid" operation ----
      if(current_operation == "regrid"){
        print("REGRID")
        
        if(regrid_option == "C"){
          resample_type <- "cubicspline"
        }else if(regrid_option == "B"){
          resample_type <- "bilinear"
        }else{
          resample_type <- "near"
        }
        
        regrid_function(master_image_path = master_image_path,
                        input_file_path = input_file_path, 
                        output_file_path = output_file_path, 
                        resample_type = resample_type)
      }
      
      ##1.3.2.3 "crop" operation ----
      if(current_operation == "crop"){
        print("CROP")
        crop_function(master_image_path = master_image_path,
                      input_file_path = input_file_path, 
                      output_file_path = output_file_path)
      }
      
      ##1.3.2.4 "smooth" operation ----
      if(current_operation == "smooth"){
        print("SMOOTH")
        smooth_spectra(input_file_path = input_file_path,
                       PRISMA_config = PRISMA_config,
                       PRISMA_wvl_info = PRISMA_wvl_info,
                       cloud_present_in_stack = cloud_present_in_stack,
                       full_230_bands = full_230_bands,
                       n_threads = n_threads,
                       output_file_path = output_file_path)
      }
      
      ##1.3.2.5 "addmetadata" operation ----
      if(current_operation == "addmetadata"){
        print("ADD PRISMA METADATA")
        add_PRISMA_metadata(input_file_path = input_file_path,
                            PRISMA_wvl_info = PRISMA_wvl_info ,
                            PRISMA_angle_info = PRISMA_angle_info,
                            PRISMA_config = PRISMA_config,
                            cloud_present_in_stack = cloud_present_in_stack,
                            output_file_path = output_file_path,
                            full_230_bands = full_230_bands)
        
      }
      
      ##1.3.2.6 "isofit" operation ----
      if(current_operation == "isofit"){
        print("ISOFIT ATMOSPHERIC CORRECTION")
        
        if(product_type == "L0"){
          stop("Generation of geometry angles not of the L0 product but to the originary L1 product")
        }
        
        if(product_type == "L2"){
          stop("The L2 product is already in reflectance, so do not need any atcor.")
        }
        
        isofit_atcor(output_file_path = output_file_path,
               input_file_path = input_file_path,
               PRISMA_wvl_info = PRISMA_wvl_info, 
               root_folder = root_folder,
               he5_path = he5_path,
               PRISMA_angle_info = PRISMA_angle_info)
        
      }
      
    }
  }
  
  
  
  
  
  
  
  
}

lapply(root_folders, PRISTAR_processing)

