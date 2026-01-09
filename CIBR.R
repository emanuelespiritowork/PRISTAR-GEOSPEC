#prova


CIBR <- function(rho_x,rho_y,rho_z,lambda_x,lambda_y,lambda_z){
  library(terra)
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\Experimental plan CIBR.docx"
  w_x <- (lambda_z-lambda_y)/(lambda_z-lambda_x)
  w_z <- (lambda_y-lambda_x)/(lambda_z-lambda_x)
  
  denom <- w_x*rho_x + w_z*rho_z
  
  denom[denom == 0] <- NA
  
  cibr_value <- 1 - (rho_y)/(denom)
  
  return(cibr_value)
}

RE_CIBR <- function(PRISMA_img){
  library(terra)
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\USGS_speclib_analyses\AbsBandIntervals_PRISMA.csv"
  rho_x <- PRISMA_img[[19]]
  rho_y <- PRISMA_img[[34]]
  rho_z <- PRISMA_img[[46]]
  
  lambda_x <- 542
  lambda_y <- 674
  lambda_z <- 796
  
  return(CIBR(rho_x = rho_x,
              rho_y = rho_y,
              rho_z = rho_z,
              lambda_x = lambda_x,
              lambda_y = lambda_y,
              lambda_z = lambda_z))
}

VegW_CIBR <- function(PRISMA_img){
  library(terra)
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\USGS_speclib_analyses\AbsBandIntervals_PRISMA.csv"
  rho_x <- PRISMA_img[[73]]
  rho_y <- PRISMA_img[[84]]
  rho_z <- PRISMA_img[[91]]
  
  lambda_x <- 1067
  lambda_y <- 1185
  lambda_z <- 1262
  
  return(CIBR(rho_x = rho_x,
              rho_y = rho_y,
              rho_z = rho_z,
              lambda_x = lambda_x,
              lambda_y = lambda_y,
              lambda_z = lambda_z))
}

LCb_CIBR <- function(PRISMA_img){
  library(terra)
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\USGS_speclib_analyses\AbsBandIntervals_PRISMA.csv"
  
  rho_y <- PRISMA_img[[176]]
  rho_z <- PRISMA_img[[190]]
  
  lambda_x <- 2010
  lambda_y <- 2094
  lambda_z <- 2206
  
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\Experimental plan CIBR.docx"
  lambda_x_right <- 2036
  lambda_x_left <- 1993
  rho_x <- PRISMA_img[[164]] + (PRISMA_img[[169]]-PRISMA_img[[164]])*(lambda_x-lambda_x_left)/(lambda_x_right-lambda_x_left)
    
  return(CIBR(rho_x = rho_x,
              rho_y = rho_y,
              rho_z = rho_z,
              lambda_x = lambda_x,
              lambda_y = lambda_y,
              lambda_z = lambda_z))
}

Soil_CIBR <- function(PRISMA_img){
  library(terra)
  #from "\\10.0.1.243\projects\CROP_NPV\4_PROCESSING\2024_Experiments\CIBR\USGS_speclib_analyses\AbsBandIntervals_PRISMA.csv"
  rho_x <- PRISMA_img[[176]]
  rho_y <- PRISMA_img[[190]]
  rho_z <- PRISMA_img[[194]]
  
  lambda_x <- 2094
  lambda_y <- 2206
  lambda_z <- 2237
  
  return(CIBR(rho_x = rho_x,
              rho_y = rho_y,
              rho_z = rho_z,
              lambda_x = lambda_x,
              lambda_y = lambda_y,
              lambda_z = lambda_z))
}


CIBR_stack <- c(LCb,RE,Soil,VegW)
names(CIBR_stack) <- c("LCb","RE","Soil","VegW")

terra::writeRaster(CIBR_stack,
                   "//10.0.1.243/nr_data/3_rs_data/PRISMA/Piacenza/PRS_L1_STD_OFFL_20250704102730_20250704102734_0001/CIBR/CIBR_stack.tif",
                   overwrite=T)

























