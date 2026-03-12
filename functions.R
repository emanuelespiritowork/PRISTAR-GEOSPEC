#_____________________________________________________________________
#check_consistency ----
#_____________________________________________________________________
check_consistency <- function(procedure_order
                              ){
  continue <- T
  er <- NULL
  
  if("inject" %in% procedure_order){
    if(procedure_order[[1]] != "inject"){
      er <- 0
      continue <- F
      stop("Inject must be the first parameter of this sequence.")
    }
    
    if(length(procedure_order) > 1){
      if(procedure_order[[2]] != "read"){
        er <- 1
        continue <- F
        stop("Read must appear after injection.")
      }
    }
  }else{
    if(procedure_order[[1]] != "read"){
      er <-  2
      continue <- F
      stop("Read must be the first parameter of this sequence.")
    }
  }
  
  return(continue)
}


#_____________________________________________________________________
#identify_product_type ----
#_____________________________________________________________________
identify_product_type <- function(he5_path
                                  ){
  if(!identical(he5_path,character(0)))
  {
    if(grepl("injected",base::basename(he5_path))){
      product_type <- "L0"
    }else{
      product_type <- base::substring(base::basename(he5_path),5,6)
    }
  }else{
    product_type <- base::substring(base::basename(he5_path),5,6)
  }
  
  return(product_type)
}

#_____________________________________________________________________
#check_folder_chain ----
#_____________________________________________________________________
check_folder_chain <- function(name_of_current_output_folder, 
                               out_folder, 
                               current_operation
                               ){
  if(name_of_current_output_folder == ""){
    name_of_current_output_folder <- paste(out_folder,current_operation, sep= "/")
  }else{
    out_folder <- name_of_current_output_folder
    name_of_current_output_folder <- paste0(name_of_current_output_folder,"_",current_operation)
  }
  
  base::dir.create(name_of_current_output_folder, recursive = T, showWarnings = F)
  
  return(name_of_current_output_folder)
}


#_____________________________________________________________________
#check_file_chain ----
#_____________________________________________________________________
check_file_chain <- function(out_folder, 
                             name_of_current_output_folder
                             ){
  
  if(name_of_current_output_folder == ""){
    folder <- out_folder
    file_input_path <- base::list.files(path = folder, pattern = "\\HCO_FULL_CLD.tif$", full.names = T, recursive = T)
    if(identical(file_input_path,character(0))){
      file_input_path <- base::list.files(path = folder, pattern = "\\HCO_FULL.tif$", full.names = T, recursive = T)
      if(identical(file_input_path,character(0))){
        print("I take ELSE")
          file_input_path <- base::list.files(path = folder, pattern = glob2rx("*.tif$"), ignore.case = T, full.names = T)
          file_input_path <- file_input_path[!(substr(basename(file_input_path),0,2) == "S2") & !(substr(basename(file_input_path),0,2) == "s2")]
      }else{
        print("I take FULL")
      }
    }else{
      print("I take FULL_CLD")
    }
  }else{
    folder <- name_of_current_output_folder
    file_input_path <- base::list.files(path = folder, pattern = ".tif$", full.names = T, recursive = T)
  }
  
  
  
  return(file_input_path)
}


#_____________________________________________________________________
#prismaread ----
#_____________________________________________________________________
prismaread_function <- function(product_type, 
                                he5_file, 
                                out_folder
                                ){
  
  if(product_type == "L1"){
    CLOUD <- T
    ATCOR <- T
    ANGLES <- T
    LATLON <- T
  }
  if(product_type == "L2"){
    CLOUD <- F
    ATCOR <- F
    ANGLES <- T
    LATLON <- T
  }
  if(product_type == "L0"){
    CLOUD <- F
    ATCOR <- F
    ANGLES <- F
    LATLON <- F
  }
  
  prismaread::pr_convert(
    ANGLES = ANGLES,
    in_file = he5_file,
    out_folder = out_folder,
    out_format = "GTiff",
    base_georef = T,
    VNIR = F,
    SWIR = F,
    FULL = T,
    fill_gaps = T,
    source = "HCO",
    join_priority = "SWIR",
    LATLON = LATLON,
    PAN = F,
    CLOUD = CLOUD,
    overwrite = T,
    ATCOR = ATCOR
  )
}

#_____________________________________________________________________
#atcor and prisma_angle ----
#_____________________________________________________________________
atcor_parameters <- function(angle_file_path
                             ){
  
  prismaread_angle_file <- utils::read.table(angle_file_path, header =T)
  PRISMA_angle_command <- base::paste0("python"," ","/space/PRISMA_angle.py")
  sensor_angles <- base::system(PRISMA_angle_command, intern = TRUE)[[2]]
  
  sensor_zenith <- strsplit(sensor_angles,",")[[1]][1]
  sensor_zenith <- gsub('^.', '', sensor_zenith)
  
  sensor_azimuth <- strsplit(sensor_angles,",")[[1]][2]
  sensor_azimuth <- gsub('^.|.$', '', sensor_azimuth)
  
  prismaread_angle_file$sensor_zenith <- as.numeric(sensor_zenith)
  prismaread_angle_file$sensor_azimuth <- as.numeric(sensor_azimuth)
  
  if(prismaread_angle_file$sensor_azimuth < 0){
    prismaread_angle_file$sensor_azimuth <- 360 + prismaread_angle_file$sensor_azimuth
  }
  
  suppressWarnings(dir.create(paste0(base::dirname(angle_file_path),"/ATCOR/")))
  
  write.table(prismaread_angle_file,paste0(base::dirname(angle_file_path),"/ATCOR/all_angles_file.csv"), 
            quote = F, 
            row.names = F,
            append = F,
            sep = ",")
  
  print("Please read the ATCOR_readme.txt file for info on angles.")
  
  atcor_readme <- "Azimuth angles range = [0,360]Deg. Zenith angles range = [0,90]Deg"
  
  write(atcor_readme,paste0(base::dirname(angle_file_path),"/ATCOR/atcor_readme.txt"))
}

#_____________________________________________________________________
#cloud mask ----
#_____________________________________________________________________
cloud_mask <- function(cloud_path, 
                       full_path
                       ){
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
  
  terra::writeRaster(x = full, 
                     filename = base::gsub("*.tif$","_CLD.tif",full_path),
                     # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                     overwrite = T
                     )
}

#_____________________________________________________________________
#coreg ----
#_____________________________________________________________________
predicted <- function(linear_model,
                      x,
                      y,
                      z
                      ){
  coefficients <- linear_model$coefficients
  num_coefficients <- as.numeric(coefficients)
  numerator <- num_coefficients[1] + num_coefficients[2]*x + num_coefficients[3]*y + num_coefficients[4]*z +
    num_coefficients[5]*x*y + num_coefficients[6]*x*z + num_coefficients[7]*y*z + num_coefficients[8]*x^2 + 
    num_coefficients[9]*y^2 + num_coefficients[10]*z^2
  denominator <- 1 - (num_coefficients[11]*x + num_coefficients[12]*y + num_coefficients[13]*z + 
                        num_coefficients[14] * x*y + num_coefficients[15]*x*z + num_coefficients[16]*y*z + 
                        num_coefficients[17] * x^2 + num_coefficients[18]*y^2 + num_coefficients[19]*z^2)
  return(numerator/denominator)
}

coregistration_to_s2 <- function(s2_file,
                                 coreg_input_path,
                                 coreg_out_folder,
                                 dem,
                                 dem_path,
                                 product_type,
                                 PRS_band_for_coreg,
                                 shift = F, 
                                 shift_x = 0, 
                                 shift_y = 0
                                 ){
  #create single layer image to coregister and change crs to EPSG:32632
  coreg_proj_path <- gsub("*.tif$","_proj.tif",coreg_input_path)
  coreg_proj_path_52 <- base::gsub("proj","proj_52",coreg_proj_path)
  
  target_epsg <- paste0("epsg:",terra::crs(terra::rast(s2_file),T,T,T)[3]$code)
  
  #project_dem if needed
  if(dem){
    dem_file <- terra::rast(dem_path)
    if(paste0("epsg:",terra::crs(dem_file, describe = T)$code) != target_epsg){
      reproject <- terra::project(dem_file, y = target_epsg)
      file.remove(dem_path)
      dem_path <- gsub(".tif","_reprojected.tif",dem_path)
      terra::writeRaster(x = reproject,
                         filename = dem_path,
                         # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                         overwrite=T
      )
    }
  }
  
  
  prisma_projected <- terra::project(x = terra::rast(coreg_input_path),
                                     y = target_epsg,
                                     method = "near",
                                     threads = T)
  
  prisma_projected_52 <- terra::subset(x = prisma_projected, 
                                       subset = as.numeric(PRS_band_for_coreg))
  
  if(shift){
    print("shift for all bands")
    shifted <- raster::shift(x = prisma_projected,
                             dx = shift_x,
                             dy = shift_y)
    
    terra::writeRaster(x = shifted, 
                       # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                       filename = gsub("*.tif$","_traslated.tif",coreg_proj_path),
                       overwrite = T
                       )
    
    coreg_proj_path <- gsub("*.tif$","_traslated.tif",coreg_proj_path)
    
    print("shift for 52 band")
    shifted <- raster::shift(x = prisma_projected_52,
                             dx = shift_x,
                             dy = shift_y)
    
    terra::writeRaster(x = shifted, 
                       # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                       filename = gsub("*.tif$","_traslated.tif",coreg_proj_path_52),
                       overwrite = T
                       )
    
    coreg_proj_path_52 <- gsub("*.tif$","_traslated.tif",coreg_proj_path_52)
    
    
  }else{
    terra::writeRaster(x = prisma_projected, 
                       filename = coreg_proj_path,
                       # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                       overwrite = T
                       )
    
    terra::writeRaster(x = prisma_projected_52, 
                       filename = coreg_proj_path_52,
                       # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                       overwrite = T)
  }
  
  # set arguments
  single_band_reference_image <- s2_file
  single_band_image_to_coregister <- coreg_proj_path_52
  multiband_image_to_coregister <- coreg_proj_path
  if(product_type == "L0"){
    arosics_local_path <- base::paste0("python"," ","/space/arosics_local_L0.py")
  }else{
    arosics_local_path <- base::paste0("python"," ","/space/arosics_local.py")
  }
  
  output_directory <- coreg_out_folder

  # setup AROSICS run command
  #IMPROVE: potrei usare AROSICS da riga di comando?
  arosics_run_command <- base::paste(arosics_local_path, "-v -r", single_band_reference_image, "-t", single_band_image_to_coregister, "-l", base::normalizePath(path=base::paste(output_directory, "/", "AROSICS_coregistration_info.json", sep=""), winslash="/", mustWork=FALSE), "-m", base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE), "-g", base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE), sep=" ")
  # run AROSICS to get GCPs
  base::system(arosics_run_command)
  # import GCPs
  #gcp_lines <- base::readLines(con=base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE))
  # generate GCP line to be used in GDAL
  
  # if(product_type == "L0"){
  #   #remove gcp outside the boundary of the image
  #   # gcp_all <- terra::extract(x = terra::rast(coreg_proj_path_52),
  #   #                y = terra::vect(base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE)),
  #   #                ID = F,
  #   #                bind = T) 
  #   #not correct because we need to use satellite coordinates
  #   vect <- terra::vect(base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE))
  #   proj_df <- as.data.frame(terra::rast(coreg_proj_path_52), xy = F, cells = T)  
  #   library(tidytable)
  #   rowcol <- tidyterra::as_coordinates(terra::rast(coreg_proj_path_52), as.raster = FALSE) |>
  #     rename(cell = cellindex, row = rowindex, col = colindex) |>
  #   #qui bisogna invertire le righe per fare in modo che dall'angolo in alto a sinistra 
  #   #si passi all'angolo in basso a sinistra
  #     group_by(col) |>
  #     mutate(row_number = row_number()) |>
  #     mutate(max_row_number = max(row_number)) |>
  #     ungroup() |>
  #     mutate(new_row = max_row_number - row) |>
  #     select(-max_row_number,-row_number,-row) |>
  #     rename(row = new_row) 
  #     #group_by(row) |> #questo se devo invertire anche le x
  #     #mutate(col_number = row_number()) |>
  #     #mutate(max_col_number = max(col_number)) |>
  #     #ungroup() |>
  #     #mutate(new_col = max_col_number - col) |>
  #     #select(-max_col_number,-col_number,-col) |>
  #     #rename(col = new_col)
  #   
  #   proj_df_rowcol <- proj_df |> left_join(rowcol)
  #   
  #   library(tidytable)
  #   col_name <- names(proj_df |> select(-cell))
  #   gcp_spec <- proj_df_rowcol |> rename(value = any_of(col_name)) |>
  #     filter(value > 0) |>
  #     na.omit() |>
  #     select(row,col) |>
  #     rename(X_IM = col, Y_IM = row)
  #   
  #   #  filter(value > 0) 
  #   # na.omit()
  #   # select(any_of(x,y))
  #   # left_join(vect)
  #   # select(X_IM,Y_IM) => così tiro fuori gcp_spec
  #   
  #   # vect_names <- names(terra::vect(base::normalizePath(path=base::paste(output_directory, "/", "points.gpkg", sep=""), winslash="/", mustWork=FALSE)))
  #   # col_name <- names(gcp_all)[!(names(gcp_all) %in% vect_names)]
  #   # gcp_spec <- gcp_all |>
  #   #   sf::st_as_sf() |>
  #   #   tidytable::select(any_of(c("X_IM","Y_IM",col_name))) |>
  #   #   rename(value = any_of(col_name)) |>
  #   #   filter(value > 0) |>
  #   #   na.omit() |>
  #   #   select(-value)
  #   
  #   gcp_read <- read.csv(base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE), 
  #                        sep = " ", 
  #                        header = F) |>
  #     rename(X_IM = V1, Y_IM = V2, X_PROJ = V3, Y_PROJ = V4) |>
  #     left_join(gcp_spec,.) |>
  #     na.omit()
  #   
  #   gcp_to_print <- terra::vect(x = gcp_read |> as.data.frame(),
  #               geom = c("X_PROJ","Y_PROJ"),
  #               crs = "epsg:32632")
  #   
  #   terra::writeVector(gcp_to_print, 
  #                      base::normalizePath(path=base::paste(output_directory, "/", "points_sel.gpkg", sep=""), winslash="/", mustWork=FALSE),
  #                      overwrite = T)
  #   
  #   write.table(gcp_read, base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE),
  #               quote = F,
  #               row.names = F,
  #               col.names = F,
  #               sep = " ",
  #               append = F)
  # }
  
  gcp_lines <- base::readLines(con=base::normalizePath(path=base::paste(output_directory, "/", "GCP.txt", sep=""), winslash="/", mustWork=FALSE))
  
  #add z coordinates to the GCP if there is a DEM
  if(dem){
    dem_projected <- terra::rast(dem_path)
    # library(tidytable)
    gcp_x_y <- read.table(text = as.character(gcp_lines), header = FALSE) |>
      dplyr::rename(easting = V3, northing = V4) |>
      dplyr::rename(X_map = V1, Y_map = V2) |>
      tidytable::mutate(ID = row_number()) |>
      data.table::as.data.table() |>
      terra::vect(geom = c("easting", "northing"), crs = target_epsg, keepgeom = T)
    
    get_altitudes <- terra::extract(x = dem_projected, y = gcp_x_y, ID = T)
    names(get_altitudes) <- c("ID","z")
    gcp <- get_altitudes |> tidytable::left_join(gcp_x_y |> data.table::as.data.table(), by = "ID") |>
      tidytable::select(-ID) |>
      tidytable::select(X_map, Y_map, easting, northing, z) |>
      data.table::as.data.table()
    write.table(gcp,base::paste(output_directory, "/", "GCP.txt", sep=""), append = F, row.names = F, col.names = F)
    list_gcp <- lapply( 1:nrow(gcp), function(n) {paste(as.character(gcp[n,]), collapse = " ")} )
    gcp_args <- base::paste(base::paste("-gcp", list_gcp, sep=" "), collapse=" ")
  }else{
    gcp_args <- base::paste(base::paste("-gcp", gcp_lines, sep=" "), collapse=" ")
  }
  
  #orthorectification function using second-order RPF with GCPs 
  if(dem){
    gcp_prova <- gcp |> tidytable::rename(j = X_map, i = Y_map, x = easting, y = northing) |> #HERE I revert (X,Y) to (i,j) since PRISMA is taking in descending orbit #LET USER CHOOSE
      tidytable::mutate(x2 = x^2, y2 = y^2, z2 = z^2, xy = x*y, yz = y*z, xz = x*z)
    
    linear_model_i <- lm(i ~ x + y + z + xy + xz + yz + x2 + y2 + z2 + I(i*x) + I(i*y) + I(i*z) + I(i*xy) + I(i*xz) + I(i*yz) +
                           I(i*x2) + I(i*y2) + I(i*z2), data = gcp_prova)
    
    linear_model_j <- lm(j ~ x + y + z + xy + xz + yz + x2 + y2 + z2 + I(j*x) + I(j*y) + I(j*z) + I(j*xy) + I(j*xz) + I(j*yz) +
                           I(j*x2) + I(j*y2) + I(j*z2), data = gcp_prova)
    
    #####
    #####get orthoprojection of PRISMA into S2 grid
    #####
    s2_raster <- terra::rast(s2_file)
    # library(tidytable)
    s2_df <- as.data.frame(s2_raster, xy =T)[,c(1,2)] |> tidytable::mutate(flag = row_number() %% 9) |>
      tidytable::filter(flag == 0) |> tidytable::select(-flag) |> as.data.frame() 
    rm(s2_raster)
    invisible(gc())
    terra::gdalCache(10000000000000)
    terra::terraOptions(memmax = 24, memfrac = 0.9)
    s2_grid_df_ext <- terra::extract(dem_projected, s2_df, xy = T) |> tidytable::select(-ID) 
    colnames <- colnames(s2_grid_df_ext)[colnames(s2_grid_df_ext) != "x" & colnames(s2_grid_df_ext) != "y"]
    s2_grid_df <- s2_grid_df_ext |> tidytable::rename(z = all_of(colnames))
    rm(s2_grid_df_ext)
    rm(s2_df)
    invisible(gc())
    s2_grid <- s2_grid_df |> tidytable::mutate(i = predicted(linear_model_i,x,y,z), j = predicted(linear_model_j,x,y,z))
    rm(s2_grid_df)
    invisible(gc())
    
    #filter the s2_grid by extent in the PRISMA image back
    #then sample the PRISMA image back using nearest neighbor
    PRISMA_grid_x_range <- nrow(prisma_projected)
    PRISMA_grid_y_range <- ncol(prisma_projected)
    s2_grid_filtered <- s2_grid |> 
      tidytable::mutate(i = round(i), j = round(j)) |>
      tidytable::filter(i > 0 & j > 0 & i < PRISMA_grid_x_range & j < PRISMA_grid_y_range) |>
      tidytable::group_by(i,j) |>
      tidytable::mutate(z_mean = mean(z), x_mean = mean(x), y_mean = mean(y)) |>
      tidytable::ungroup() |>
      tidytable::select(-x,-y,-z) |>
      tidytable::rename(z = z_mean, x = x_mean, y = y_mean) |>
      unique()
    rm(s2_grid)
    invisible(gc())
    
    
    saveRDS(s2_grid_filtered, paste0(coreg_out_folder,"/s2_grid_filtered.rds"))
    #as.numeric(prisma_projected[raster::cellFromRowCol(prisma_projected, c(800,700,600), c(800,700,600))])
    #s2_grid_filtered$value <- as.numeric(prisma_projected[raster::cellFromRowCol(prisma_projected, s2_grid_filtered$i, s2_grid_filtered$j)])
    #s2_prova <- s2_grid_filtered |> filter(i > 800 & i < 820 & j > 800 & j < 820)
    
    #nome_colonna <- paste(colnames(s2_grid_sampled |> select(-x,-y)),"_mean", sep = "")
    #nome_vecchio_colonna <- colnames(s2_grid_sampled |> select(-x,-y))
    
    #s2_grid_filtered <- s2_grid_filtered |> filter(i > 600 & i < 900 & j > 600 & j < 900)
    ###############
    s2_grid_spectra <- prisma_projected[raster::cellFromRowCol(prisma_projected, s2_grid_filtered$i, s2_grid_filtered$j)]
    saveRDS(s2_grid_spectra,paste0(coreg_out_folder,"/s2_grid_spectra.rds"))
    invisible(gc())
    s2_grid_sampled <- cbind(s2_grid_filtered$i,s2_grid_filtered$j,s2_grid_spectra) |>
      tidytable::rename(i = any_of(c("s2_grid_filtered$i")), j = any_of(c("s2_grid_filtered$j"))) |>
      tidytable::right_join(s2_grid_filtered) |>
      tidytable::select(-i,-j,-z)
    saveRDS(s2_grid_sampled,paste0(coreg_out_folder,"/s2_grid_sampled.rds"))
    rm(s2_grid_spectra)
    rm(s2_grid_filtered)
    invisible(gc())
    #library(tidytable)
    s2_grid_sampled <- readRDS(paste0(coreg_out_folder,"/s2_grid_sampled.rds"))
    s2_grid_sampled_dt <- s2_grid_sampled |> data.table::as.data.table() 
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
    output_file <- base::paste0(output_directory,"/",gsub(".tif","_ortho.tif",basename(coreg_input_path)))
    # 
    # out_file <- set_wvl_time(output_file_path = output_file,
    #                        output_file = output_focal_again,
    #                        wvl = all_wvl)
    # 
    terra::writeRaster(x = output_focal_again,
                      filename = output_file,
                      # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                      overwrite = T
                      )
    
    rm(output_focal_again)
    invisible(gc())
    
    #s2_grid_x_range <- nrow(s2_raster)
    #s2_grid_y_range <- ncol(s2_raster)
    #s2_grid_number_of_cells <- s2_grid_x_range*s2_grid_y_range
    #s2_grid_list_of_cells <- lapply(1:s2_grid_number_of_cells, function(pixel){raster::xyFromCell(s2_raster,pixel)})
    
      #raster::xyFromCell(dem_projected, raster::cellFromRowCol(dem_projected, 800, 800))
    #file.remove(paste0(coreg_out_folder,"/raster_not_focal.tif"))
  }else{
    output_file <- base::paste0(output_directory,"/",gsub(".tif","_coreg.tif",basename(coreg_input_path)))
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
  file.remove(coreg_proj_path_52)
  file.remove(coreg_proj_path)
  
}

#_____________________________________________________________________
#regrid ----
#_____________________________________________________________________
regrid_function <- function(master_image_path, name_of_current_output_folder, regrid_input_path, resample_type){
  master <- terra::rast(master_image_path)
  slave <- terra::rast(regrid_input_path)
  
  if(terra::crs(master, describe = T)$code != terra::crs(slave, describe = T)$code){
    reproject <- terra::project(slave, y = paste0("EPSG:",terra::crs(master, describe = T)$code))
    slave <- reproject
  }
  
  output_file <- base::paste0(name_of_current_output_folder,"/",gsub(".tif","_regrid.tif",basename(regrid_input_path)))
  
  terra::extend(x = master,
                y = slave,
                fill = 0,
                filename = base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"),
                overwrite = T)
  
  terra::resample(x = slave,
                  y = terra::rast(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif")),
                  method = resample_type,
                  threads = T,
                  by_util = T,
                  filename = output_file,
                  # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                  overwrite = T)
  
  file.remove(base::paste0(name_of_current_output_folder,"/PRISMA_extend.tif"))
  
  # out_file <- set_wvl_time(output_file_path = output_file,
  #                          output_file = regridded,
  #                          wvl = all_wvl)
  # 
  # terra::writeRaster(x = out_file,
  #                    filename = output_file,
  #                    overwrite = T)
}

#_____________________________________________________________________
#crop ----
#_____________________________________________________________________
crop_function <- function(master_image_path, name_of_current_output_folder, crop_input_path){
  master <- terra::rast(master_image_path)
  slave <- terra::rast(crop_input_path)
  
  if(terra::crs(master, describe = T)$code != terra::crs(slave, describe = T)$code){
    reproject <- terra::project(slave, y = paste0("EPSG:",terra::crs(master, describe = T)$code))
    slave <- reproject
  }
  
  output_file <- base::paste0(name_of_current_output_folder,"/",gsub(".tif","_crop.tif",basename(crop_input_path)))
  
  
  terra::crop(x = slave,
              y = master,
              filename = output_file,
              # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
              overwrite = T)
  
  # out_file <- set_wvl_time(output_file_path = output_file,
  #                          output_file = cropped,
  #                          wvl = all_wvl)
  # 
  # terra::writeRaster(x = out_file,
  #                    filename = output_file,
  #                    overwrite = T)
  
}

#_____________________________________________________________________
#smooth ----
#_____________________________________________________________________
smooth_spectra <- function(smooth_file_path,
                           PRISMA_config,
                           PRISMA_bad_bands_table,
                           name_of_current_output_folder,
                           cloud_smooth,
                           full_230_bands,
                           n_threads = 1
                           ){
  output_file <- base::paste0(name_of_current_output_folder,"/",gsub(".tif","_smooth.tif",basename(smooth_file_path)))
  
  input_bad_bands <- PRISMA_bad_bands_table$band
  
  input_wvl <- PRISMA_config$center
  
  selection_vector <- 1
  
  if(full_230_bands){
    selection_vector <- c(0,1)
  }
  
  #which output bands
  output_wvl <- PRISMA_config |>
    tidytable::filter(BND_SEL %in% selection_vector) |>
    tidytable::pull(center)
  
  #print("Define spline function")
  
  
  
  #print("Read terra image")
  
  terra_image <- terra::rast(smooth_file_path)
  
  if(cloud_smooth){
    cloud <- terra::subset(terra_image,subset = 231, negate = F)
    terra_image_sub <- terra::subset(terra_image,subset = 231,negate=T)
  }else{
    terra_image_sub <- terra_image
  }
  
  #print("Apply smoothing")
  
  #terra::terraOptions(memmin = 30, print=T, progress = 1, memfrac = 0.8, verbose = T)
  
  #terra::gdalCache(1000000)
  # bb <- c()
  # for(i in 1:(terra::nrow(terra_image_sub)*terra::ncol(terra_image_sub))){
  #   df <- as.data.frame(terra_image_sub[i]) |> pivot_longer(names_to = "band", values_to = "value")
  #   df_na <- df |> filter(is.na(value))
  #   bb <- c(bb,df_na$band)
  #   bb <- unique(bb)  
  # }
  
  smoothed_file <- terra::app(
    x = terra_image_sub,
    fun = spline_fun,
    band_center_input = input_wvl, 
    bad_bands_pos = input_bad_bands, 
    band_center_output = output_wvl,
    df=40,
    
    cores = n_threads
  )
  
  if(cloud_smooth){
    out_file <- c(smoothed_file,cloud)
  }else{
    out_file <- smoothed_file
  }
  
  terra::writeRaster(x = out_file,
                     filename = output_file,
                     # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                     overwrite = T
    )
  # return(NULL)
}

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

#_____________________________________________________________________
#add_PRISMA_metadata ----
#_____________________________________________________________________
add_PRISMA_metadata <- function(name_of_current_output_folder,
                                metadata_file_path,
                                all_wvl,
                                cloud
                                ){
  output_file <- terra::rast(metadata_file_path)
  
  if(cloud){
    base::names(output_file) <- c(all_wvl,"cloud_mask")
  }else{
    base::names(output_file) <- all_wvl
  }
  
  if(grepl(pattern= "OFFL", x = basename(metadata_file_path))){
    terra::time(output_file) <- rep(as.Date(substring(stringr::str_match(basename(metadata_file_path), "OFFL_\\s*(.*?)\\s*_HCO")[2],1,8),format = "%Y%m%d"),
                                 terra::nlyr(output_file))
  }else if(grepl(pattern= "STD", x = basename(metadata_file_path))){
    terra::time(output_file) <- rep(as.Date(substring(stringr::str_match(basename(metadata_file_path), "STD_\\s*(.*?)\\s*_HCO")[2],1,8),format = "%Y%m%d"),
                                 terra::nlyr(output_file))
  }else{
    stop("ERROR. The he5 file has no standard name (it does not contain neither STD nor OFFL).")
  }
  
  output_file_path <- base::paste0(name_of_current_output_folder,"/",gsub(".tif","_addmetadata.tif",basename(metadata_file_path)))
  
  terra::writeRaster(x = output_file,
                     filename = output_file_path,
                     # wopt = base::list(gdal = c("COMPRESS=LZW", "TILED=YES")),
                     overwrite = T)
  
}
