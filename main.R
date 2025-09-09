#designed in Jun 2025
#author: Emanuele Spirito
#site: CNR-IREA-MI

#citations: 

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics
#https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busetto @ CNR-IREA for prismaread package 
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
regrid_option <- "N" #can be N for near, B for bilinear, C for cubic

#for expert users:
#procedure_order <- c("read","cloud","coreg","atcor","regrid","smooth")
procedure_order <- c("regrid")
#elements: read, atcor, cloud, coreg, regrid, smooth

#_____________________________________________________________________
#setup folders ----
#_____________________________________________________________________
source("/space/functions.R")
out_folder  <-  "/space/put_PRISMA_he5_and_S2_tif_here/"
s2_folder <- out_folder
he5_folder <- out_folder
base::setwd(base::dirname(out_folder))

he5_file <- base::list.files(path = out_folder, pattern = "\\.he5$", ignore.case = T, full.names = T)
product_type <- base::substring(base::basename(he5_file),5,6)

#_____________________________________________________________________
#work ----
#_____________________________________________________________________
number_of_operations <- length(procedure_order)
name_of_current_output_folder <- ""
for(index_of_operations in 1:number_of_operations){
  current_operation <- procedure_order[index_of_operations]
  
  #procedures that are the start of a sequence
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
  if(current_operation == "coreg"){
    #chain part
    print("COREG")
    if(name_of_current_output_folder == ""){
      name_of_current_output_folder <- paste0(out_folder,current_operation)
    }else{
      out_folder <- name_of_current_output_folder
      name_of_current_output_folder <- paste0(name_of_current_output_folder,"_",current_operation)
    }
    
    base::dir.create(name_of_current_output_folder, recursive = T, showWarnings = F)
    
    #coregistration part
    s2_file <- base::list.files(path = s2_folder, pattern = glob2rx("S2*.tif$"), ignore.case = T, full.names = T)
    coreg_input_path <- base::list.files(path = out_folder, pattern = "\\HCO_FULL_CLD.tif$", full.names = T)
    coreg_proj_path <- base::gsub("FULL_CLD","FULL_CLD_proj",coreg_input_path)
    if(identical(coreg_input_path,character(0))){
      coreg_input_path <- base::list.files(path = out_folder, pattern = "\\HCO_FULL.tif$", full.names = T)
      coreg_proj_path <- base::gsub("FULL","FULL_proj",coreg_input_path)
      if(identical(coreg_input_path,character(0))){
        print("I take ELSE")
        coreg_input_path <- base::list.files(path = out_folder, pattern = glob2rx("*.tif$"), ignore.case = T, full.names = T)
        coreg_input_path <- coreg_input_path[!substr(basename(coreg_input_path),0,2) == "S2"]
        coreg_proj_path <- gsub(".tif","_proj.tif",coreg_input_path)
      }else{
        print("I take FULL")
      }
    }else{
      print("I take FULL_CLD")
    }
    
    coregistration_to_s2(s2_file,coreg_input_path,coreg_proj_path,name_of_current_output_folder)
  }
  if(current_operation == "regrid"){
    #chain part
    print("REGRID")
    if(name_of_current_output_folder == ""){
      name_of_current_output_folder <- paste0(out_folder,current_operation)
    }else{
      out_folder <- name_of_current_output_folder
      name_of_current_output_folder <- paste0(name_of_current_output_folder,"_",current_operation)
    }
    
    base::dir.create(name_of_current_output_folder, recursive = T, showWarnings = F)
    
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
    
    regrid_input_path <- base::list.files(out_folder, full.names = T, pattern = "\\.tif$")
    #regrid_input_path <- base::list.files(out_folder, full.names = T, pattern = "\\.bsq$")
    
    regrid_function(master_image_path, name_of_current_output_folder, regrid_input_path)
  }
  if(current_operation == "smooth"){
    #chain part
    print("SMOOTH")
    if(name_of_current_output_folder == ""){
      name_of_current_output_folder <- paste0(out_folder,current_operation)
    }else{
      out_folder <- name_of_current_output_folder
      name_of_current_output_folder <- paste0(name_of_current_output_folder,"_",current_operation)
    }
    
    base::dir.create(name_of_current_output_folder, recursive = T, showWarnings = F)
    
    #smooth part
    smoothing_out <-  base::paste0(name_of_current_output_folder, "/PRISMA_smoothed.tif")
    
    if("cloud" %in% procedure_order){
      cloud_smooth <- T
    }else{
      cloud_smooth <- F
    }
    
    terra_image_path <- base::list.files(out_folder, pattern = "\\.tif$", full.names = T)
    
    smooth_spectra(terra_image_path,smoothing_out,cloud_smooth)
  }
  
}






