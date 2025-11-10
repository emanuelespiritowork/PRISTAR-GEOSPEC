
#_____________________________________________________________________
#prismaread ----
#_____________________________________________________________________
prismaread_function <- function(product_type, he5_file){
  
  if(product_type == "L1"){
    CLOUD <- T
    ATCOR <- T
  }
  if(product_type == "L2"){
    CLOUD <- F
    ATCOR <- F
  }
  if(product_type == "L0"){
    CLOUD <- F
    ATCOR <- T
  }
  
  prismaread::pr_convert(
    ANGLES = T,
    in_file = he5_file,
    out_folder = paste0(base::dirname(he5_file),"/"),
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
}

#_____________________________________________________________________
#atcor and prisma_angle ----
#_____________________________________________________________________
atcor_parameters <- function(angle_file_path){
  
  prismaread_angle_file <- utils::read.table(angle_file_path, header =T)
  PRISMA_angle_command <- base::paste0("python"," ",getwd(),"/PRISMA_angle.py")
  sensor_angles <- base::system(PRISMA_angle_command, intern = TRUE)
  
  sensor_zenith <- strsplit(sensor_angles,",")[[1]][1]
  sensor_zenith <- gsub('^.', '', sensor_zenith)
  
  sensor_azimuth <- strsplit(sensor_angles,",")[[1]][2]
  sensor_azimuth <- gsub('^.|.$', '', sensor_azimuth)
  
  prismaread_angle_file$sensor_zenith <- as.numeric(sensor_zenith)
  prismaread_angle_file$sensor_azimuth <- as.numeric(sensor_azimuth)
  
  write.csv(prismaread_angle_file,paste0(base::dirname(angle_file_path),"/ATCOR/all_angles_file.csv"), quote = F, row.names = F)
}

cloud_mask <- function(cloud_path, full_path){
  cloud <- terra::rast(cloud_path)
  cloud <- terra::subst(cloud, NA, 1)
  terra::plot(cloud, range = c(0,1))
  
  cloud_dil <- spatialist::erodil_raster(raster = cloud, 
                                         width = c(3,3), 
                                         type = "box", 
                                         erosion = T, 
                                         dilation = T,
                                         erosion_first = T, 
                                         nt = 1)
  
  terra::plot(cloud_dil)
  
  full <- terra::rast(full_path)
  
  terra::add(full) <- cloud_dil
  
  terra::writeRaster(full, 
                     base::gsub("_HCO_FULL.tif","_HCO_FULL_CLD.tif",full_path),
                     overwrite = T)
}

#_____________________________________________________________________
#coreg ----
#_____________________________________________________________________
coregistration_to_s2 <- function(s2_file,coreg_input_path,coreg_proj_path,coreg_out_folder,dem,dem_path,product_type){
  #print("DEM is")
  #print(dem)
  #create single layer image to coregister and change crs to EPSG:32632
  target_epsg <- paste0("epsg:",terra::crs(terra::rast(s2_file),T,T,T)[3]$code)
  
  prisma_projected <- terra::project(x = terra::rast(coreg_input_path),
                                     y = target_epsg,
                                     method = "near")
  saveRDS(prisma_projected, paste0(coreg_out_folder,"/prisma_projected.rds"))
  terra::writeRaster(prisma_projected, 
                     coreg_proj_path,
                     overwrite = T)
  terra::subset(x = prisma_projected, 
                subset = 52, 
                filename = base::gsub("proj","proj_52",coreg_proj_path), 
                overwrite = T)
  
  # set arguments
  single_band_reference_image <- s2_file
  single_band_image_to_coregister <- base::list.files(base::dirname(coreg_input_path),"\\proj_52.tif$", full.names = T)
  multiband_image_to_coregister <- base::list.files(base::dirname(coreg_input_path),"\\proj.tif$", full.names = T)
  if(product_type == "L0"){
    arosics_local_path <- base::paste0("python"," ",base::getwd(),"/arosics_local_L0.py")
  }else{
    arosics_local_path <- base::paste0("python"," ",base::getwd(),"/arosics_local.py")
  }
  
  output_directory <- coreg_out_folder

  # setup AROSICS run command
  #IMPROVE: potrei usare AROSICS da riga di comando?
  arosics_run_command <- base::paste(arosics_local_path, "-v -r", single_band_reference_image, "-t", single_band_image_to_coregister, "-l", base::normalizePath(path=base::paste(output_directory, "/", "AROSICS_coregistration_info.json", sep=""), winslash="/", mustWork=FALSE), "-m", base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE), "-g", base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE), sep=" ")
  # run AROSICS to get GCPs
  base::system(arosics_run_command)
  # import GCPs
  gcp_lines <- base::readLines(con=base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE))
  # generate GCP line to be used in GDAL
  
  #add z coordinates to the GCP if there is a DEM
  if(dem){
    dem_projected <- terra::rast(dem_path)
    library(tidytable)
    gcp_x_y <- read.table(text = as.character(gcp_lines), header = FALSE) %>%
      dplyr::rename(easting = V3, northing = V4) %>%
      dplyr::rename(X_map = V1, Y_map = V2) %>%
      tidytable::mutate(ID = row_number()) %>%
      data.table::as.data.table() %>%
      terra::vect(., geom = c("easting", "northing"), crs = target_epsg, keepgeom = T)
    
    get_altitudes <- terra::extract(x = dem_projected, y = gcp_x_y, ID = T)
    names(get_altitudes) <- c("ID","z")
    gcp <- get_altitudes %>% left_join(gcp_x_y %>% data.table::as.data.table(), by = "ID") %>%
      select(-ID) %>%
      select(X_map, Y_map, easting, northing, z) %>%
      data.table::as.data.table()
    write.table(gcp,base::paste(output_directory, "/", "GCP.txt", sep=""), append = F, row.names = F, col.names = F)
    list_gcp <- lapply( 1:nrow(gcp), function(n) {paste(as.character(gcp[n,]), collapse = " ")} )
    gcp_args <- base::paste(base::paste("-gcp", list_gcp, sep=" "), collapse=" ")
  }else{
    gcp_args <- base::paste(base::paste("-gcp", gcp_lines, sep=" "), collapse=" ")
  }
  
  #orthorectification function using second-order RPF with GCPs 
  if(dem){
    gcp_prova <- gcp %>% rename(j = X_map, i = Y_map, x = easting, y = northing) %>% #HERE I revert (X,Y) to (i,j) since PRISMA is taking in descending orbit #LET USER CHOOSE
      mutate(x2 = x^2, y2 = y^2, z2 = z^2, xy = x*y, yz = y*z, xz = x*z)
    
    linear_model_i <- lm(i ~ x + y + z + xy + xz + yz + x2 + y2 + z2 + I(i*x) + I(i*y) + I(i*z) + I(i*xy) + I(i*xz) + I(i*yz) +
                           I(i*x2) + I(i*y2) + I(i*z2), data = gcp_prova)
    
    linear_model_j <- lm(j ~ x + y + z + xy + xz + yz + x2 + y2 + z2 + I(j*x) + I(j*y) + I(j*z) + I(j*xy) + I(j*xz) + I(j*yz) +
                           I(j*x2) + I(j*y2) + I(j*z2), data = gcp_prova)
    
    predicted_i <- function(x,y,z){
      coefficients <- linear_model_i$coefficients
      num_coefficients <- as.numeric(coefficients)
      numerator <- num_coefficients[1] + num_coefficients[2]*x + num_coefficients[3]*y + num_coefficients[4]*z +
        num_coefficients[5]*x*y + num_coefficients[6]*x*z + num_coefficients[7]*y*z + num_coefficients[8]*x^2 + 
        num_coefficients[9]*y^2 + num_coefficients[10]*z^2
      denominator <- 1 - (num_coefficients[11]*x + num_coefficients[12]*y + num_coefficients[13]*z + 
                            num_coefficients[14] * x*y + num_coefficients[15]*x*z + num_coefficients[16]*y*z + 
                            num_coefficients[17] * x^2 + num_coefficients[18]*y^2 + num_coefficients[19]*z^2)
      return(numerator/denominator)
    }
    
    predicted_j <- function(x,y,z){
      coefficients <- linear_model_j$coefficients
      num_coefficients <- as.numeric(coefficients)
      numerator <- num_coefficients[1] + num_coefficients[2]*x + num_coefficients[3]*y + num_coefficients[4]*z +
        num_coefficients[5]*x*y + num_coefficients[6]*x*z + num_coefficients[7]*y*z + num_coefficients[8]*x^2 + 
        num_coefficients[9]*y^2 + num_coefficients[10]*z^2
      denominator <- 1 - (num_coefficients[11]*x + num_coefficients[12]*y + num_coefficients[13]*z + 
                            num_coefficients[14] * x*y + num_coefficients[15]*x*z + num_coefficients[16]*y*z + 
                            num_coefficients[17] * x^2 + num_coefficients[18]*y^2 + num_coefficients[19]*z^2)
      return(numerator/denominator)
    }
    
    #####
    #####get orthoprojection of PRISMA into S2 grid
    #####
    s2_raster <- terra::rast(s2_file)
    library(tidytable)
    s2_df <- as.data.frame(s2_raster, xy =T)[,c(1,2)] %>% mutate(flag = row_number() %% 9) %>%
      filter(flag == 0) %>% select(-flag) %>% as.data.frame() 
    rm(s2_raster)
    invisible(gc())
    terra::gdalCache(10000000000000)
    terra::terraOptions(memmax = 24, memfrac = 0.9)
    s2_grid_df <- terra::extract(dem_projected, s2_df, xy = T) %>% select(-ID) %>% 
      rename(z = colnames(.)[colnames(.) != "x" & colnames(.) != "y"])
    rm(s2_df)
    invisible(gc())
    s2_grid <- s2_grid_df %>% mutate(i = predicted_i(x,y,z), j = predicted_j(x,y,z))
    rm(s2_grid_df)
    invisible(gc())
    
    #filter the s2_grid by extent in the PRISMA image back
    #then sample the PRISMA image back using nearest neighbor
    PRISMA_grid_x_range <- nrow(prisma_projected)
    PRISMA_grid_y_range <- ncol(prisma_projected)
    s2_grid_filtered <- s2_grid %>% 
      mutate(i = round(i), j = round(j)) %>%
      filter(i > 0 & j > 0 & i < PRISMA_grid_x_range & j < PRISMA_grid_y_range) %>%
      group_by(i,j) %>%
      mutate(z_mean = mean(z), x_mean = mean(x), y_mean = mean(y)) %>%
      ungroup() %>%
      select(-x,-y,-z) %>%
      rename(z = z_mean, x = x_mean, y = y_mean) %>%
      unique()
    rm(s2_grid)
    invisible(gc())
    
    
    saveRDS(s2_grid_filtered, paste0(coreg_out_folder,"/s2_grid_filtered.rds"))
    #as.numeric(prisma_projected[raster::cellFromRowCol(prisma_projected, c(800,700,600), c(800,700,600))])
    #s2_grid_filtered$value <- as.numeric(prisma_projected[raster::cellFromRowCol(prisma_projected, s2_grid_filtered$i, s2_grid_filtered$j)])
    #s2_prova <- s2_grid_filtered %>% filter(i > 800 & i < 820 & j > 800 & j < 820)
    
    #nome_colonna <- paste(colnames(s2_grid_sampled %>% select(-x,-y)),"_mean", sep = "")
    #nome_vecchio_colonna <- colnames(s2_grid_sampled %>% select(-x,-y))
    
    #s2_grid_filtered <- s2_grid_filtered %>% filter(i > 600 & i < 900 & j > 600 & j < 900)
    ###############
    s2_grid_spectra <- prisma_projected[raster::cellFromRowCol(prisma_projected, s2_grid_filtered$i, s2_grid_filtered$j)]
    saveRDS(s2_grid_spectra,paste0(coreg_out_folder,"/s2_grid_spectra.rds"))
    invisible(gc())
    s2_grid_sampled <- cbind(s2_grid_filtered$i,s2_grid_filtered$j,s2_grid_spectra) %>%
      rename(i = any_of(c("s2_grid_filtered$i")), j = any_of(c("s2_grid_filtered$j"))) %>%
      left_join(s2_grid_filtered,.) %>%
      select(-i,-j,-z)
    saveRDS(s2_grid_sampled,paste0(coreg_out_folder,"/s2_grid_sampled.rds"))
    rm(s2_grid_spectra)
    rm(s2_grid_filtered)
    invisible(gc())
    #library(tidytable)
    s2_grid_sampled <- readRDS(paste0(coreg_out_folder,"/s2_grid_sampled.rds"))
    s2_grid_sampled_dt <- s2_grid_sampled %>% data.table::as.data.table() 
    rm(s2_grid_sampled)
    invisible(gc())
    
    #with all pixels
    #empty_raster <- terra::rast(crs = target_epsg, extent = terra::ext(min(s2_grid_sampled_dt$x),max(s2_grid_sampled_dt$x),min(s2_grid_sampled_dt$y),max(s2_grid_sampled_dt$y)),
    #                            resolution = c(30,30), nlyrs = 230)
    
    #with x y distribution
    quantile <- 0.01
    nlyrs <- terra::nlyr(prisma_projected)
    empty_raster <- terra::rast(crs = target_epsg, extent = terra::ext(quantile(s2_grid_sampled_dt$x,quantile),quantile(s2_grid_sampled_dt$x,1-quantile),quantile(s2_grid_sampled_dt$y,quantile),quantile(s2_grid_sampled_dt$y,1-quantile)),
                                resolution = terra::res(prisma_projected), nlyrs = nlyrs)
    
    #empty_raster <- terra::rast(crs = target_epsg, extent = terra::ext(640000,645000,5170000,5180000),
    #                            resolution = c(30,30), nlyrs = 230)
    s2_grid_sampled_vect <- terra::vect(s2_grid_sampled_dt, crs = target_epsg, geom = c("x","y"))
    
    names <- colnames(s2_grid_sampled_dt)[colnames(s2_grid_sampled_dt) != "x" & colnames(s2_grid_sampled_dt) != "y"]#UPDATE: use not column position but names
    rm(s2_grid_sampled_dt)
    invisible(gc())
    output_raster <- terra::rasterize(s2_grid_sampled_vect, empty_raster, field = names)
    rm(s2_grid_sampled_vect)
    rm(empty_raster)
    invisible(gc())
    #terra::writeRaster(output_raster, paste0(coreg_out_folder,"/raster_not_focal.tif"), overwrite = T)
    
    #fill NA with bilinear value
    output_focal <- terra::focal(output_raster, w=3, fun=mean, na.policy="only", na.rm=T) #FIX: use focal only if it finds at least two values 
    output_focal_again <- terra::focal(output_focal, w=3, fun=mean, na.policy="only", na.rm=T)
    rm(output_focal)
    invisible(gc())
    terra::writeRaster(output_focal_again, paste0(coreg_out_folder,"/raster_focal.tif"), overwrite = T)
    rm(output_focal_again)
    invisible(gc())
    
    #s2_grid_x_range <- nrow(s2_raster)
    #s2_grid_y_range <- ncol(s2_raster)
    #s2_grid_number_of_cells <- s2_grid_x_range*s2_grid_y_range
    #s2_grid_list_of_cells <- lapply(1:s2_grid_number_of_cells, function(pixel){raster::xyFromCell(s2_raster,pixel)})
    
      #raster::xyFromCell(dem_projected, raster::cellFromRowCol(dem_projected, 800, 800))
    #file.remove(paste0(coreg_out_folder,"/raster_not_focal.tif"))
  }else{
    output_file <- base::paste0(output_directory,"/prs_crs_translate_warp.tif")
    # create VRT with GCP
    GDAL_VRT_run_command <- base::paste("gdal_translate -q -of VRT", gcp_args, multiband_image_to_coregister, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=FALSE), sep=" ")
    writeLines(GDAL_VRT_run_command, base::paste(output_directory, "/", "GDAL_VRT_run_command.sh", sep=""))
    base::system(paste0("sh ", base::paste(output_directory, "/", "GDAL_VRT_run_command.sh", sep="")))
    
    #GDAL_VRT_run_command <- base::paste("gdal_translate -q -of VRT", gcp_args, multiband_image_to_coregister, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=FALSE), sep=" ")
    #base::system(GDAL_VRT_run_command)
    # warp input image using second order polynomial
    GDAL_WARP_run_command <- base::paste("gdalwarp -q -of GTiff -r near -order 2 -tap -tr 30 30 -t_srs", target_epsg, base::normalizePath(path=base::paste(output_directory, "/", "multiband_file_with_GCP.vrt", sep=""), winslash="/", mustWork=TRUE), output_file, sep=" ")
    base::system(GDAL_WARP_run_command)
  }
  
  
  
  #AAAAA ALTRA IDEA: se utilizzassi prima i tie points come quota per fare correzione dem e poi usassi invece riconoscimento 2D?
  #in questo modo prima correggo gli artefatti dovuti alla visualizzazione del satellite e poi coregistro le due immagini
  
  #convert DEM from geoid to ellipsoidal
  #gdalwarp -s_srs "+proj=longlat +datum=WGS84 +no_defs +geoidgrids=/space/put_PRISMA_he5_and_S2_tif_here/DEM/egm96_15.gtx" -t_srs "+proj=longlat +datum=WGS84 +no_def" /space/put_PRISMA_he5_and_S2_tif_here/DEM/Senales_2024_02.tif /space/put_PRISMA_he5_and_S2_tif_here/DEM/Senales_2024_02_ellipsoid.tif
  
  
}

#_____________________________________________________________________
#regrid ----
#_____________________________________________________________________
regrid_function <- function(master_image_path, name_of_current_output_folder, regrid_input_path, resample_type){
  terra::extend(x = terra::rast(master_image_path),
                y = terra::rast(regrid_input_path),
                fill = 0,
                filename = base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"),
                overwrite = T)
  
  terra::resample(x = terra::rast(regrid_input_path),
                  y = terra::rast(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif")),
                  method = resample_type,
                  threads = T,
                  by_util = T,
                  filename = base::paste0(name_of_current_output_folder,"/PRISMA_resample.tif"),
                  overwrite = T)
  
  file.remove(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"))
  
}

#_____________________________________________________________________
#crop ----
#_____________________________________________________________________
crop_function <- function(master_image_path, name_of_current_output_folder, crop_input_path){
  terra::crop(x = terra::rast(crop_input_path),
              y = terra::rast(master_image_path),
              filename = base::paste0(name_of_current_output_folder,"/PRISMA_crop.tif"),
              overwrite = T)
}

#_____________________________________________________________________
#smooth ----
#_____________________________________________________________________
smooth_spectra <- function(terra_image_path,smoothing_out,cloud_smooth){
  PRISMA_config <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_spectral_configuration.csv")) %>%
    tidytable::mutate(band_row = tidytable::row_number()) 
  
  PRISMA_bad_bands_table <- tidytable::fread(base::paste0(base::getwd(),"/PRISMA_band_selections.csv")) %>%
    tidytable::filter(BB_SUPER_V3 == 1)
  
  input_bad_bands <- PRISMA_bad_bands_table$band
  
  input_wvl <- PRISMA_config$center
  
  output_wvl <- PRISMA_config %>%
    tidytable::filter(BND_SEL  == 1) %>%
    tidytable::pull(center)
  
  print("Define spline function")
  
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
  
  print("Read terra image")
  
  terra_image <- terra::rast(terra_image_path)
  
  if(cloud_smooth){
    terra_image_sub <- terra::subset(terra_image,subset = 231,negate=T)
  }else{
    terra_image_sub <- terra_image
  }
  
  print("Apply smoothing")
  
  terra::terraOptions(memmin = 30, print=T, progress = 1, memfrac = 0.8, verbose = T)
  
  terra::gdalCache(1000000)
  
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