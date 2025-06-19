rm(list = ls())


library(here)
library(terra)
library(tidytable)

input_image_path <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2024/L2C/PRS_L2C_STD_20240627101417_20240627101422_0001/coreg/prs_warp__shifted_to__S2_20240626_10m_B08.bsq"

# 20240407: "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2024/L2C/PRS_L2C_STD_20240407101730_20240407101734_0001/coreg/prs_warp__shifted_to__S2_20240412_10m_B08.bsq"
# 20240512: "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2024/L2C/PRS_L2C_STD_20240512102102_20240512102107_0001/coreg/prs_warp__shifted_to__S2_20240512_10m_B08.bsq"
# 20240610: "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2024/L2C/PRS_L2C_STD_20240610102101_20240610102105_0001/coreg/prs_warp__shifted_to__S2_20240512_10m_B08.bsq"
# 20240627: "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/2024/L2C/PRS_L2C_STD_20240627101417_20240627101422_0001/coreg/prs_warp__shifted_to__S2_20240626_10m_B08.bsq"

image_date <- "20240627" #"20240407"

output_folder <- "E:/PRISMA_images_smoothed/"
output_path <-  paste0(output_folder, image_date, "_PRISMA_smoothed.tif")

PRISMA_config <- fread(here("config_data", "FWHM", "PRISMA_spectral_configuration.csv")) %>%
    mutate(band_row = row_number()) 

PRISMA_bad_bands_table <- fread(here("config_data", "PRISMA_band_selections.csv")) %>%
    filter(BB_SUPER_V3 == 1)

input_bad_bands <- PRISMA_bad_bands_table$band

input_wvl <- PRISMA_config$center

output_wvl <- PRISMA_config %>%
    filter(BND_SEL  == 1) %>%
    pull(center)


spline_fun <- function(pixel, band_center_input, bad_bands_pos, band_center_output, df = 40) {
    # togliamo le bande cattive
    ref_valid <- pixel[-bad_bands_pos]
    wvl_valid <- band_center_input[-bad_bands_pos]
    
    if (any(is.na(ref_valid))) {
        return(rep(NA_real_, length(band_center_output)))
    }
    
    sp <- stats::smooth.spline(x = wvl_valid, y = ref_valid, df = df)
    y_smooth <- predict(sp, x = band_center_output)$y
    # limitiamo a zero
    y_smooth[y_smooth < 0] <- 0
    return(y_smooth)
}


terra_image <- rast(input_image_path)


output_image <- app(
    x = terra_image,
    fun = spline_fun,
    band_center_input = input_wvl, 
    bad_bands_pos = input_bad_bands, 
    band_center_output = output_wvl,
    df=40,
    
    cores = 7,                     
    filename = output_path,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
)



output_image <- app(terra_image, function(i, band_center_input, bad_bands_pos, band_center_output){
   
    center_to_spline <- band_center_input[-bad_bands_pos]
    ref_to_spline <- i[-bad_bands_pos]
    
    if (any(is.na(ref_to_spline))) {
        return(rep(NA, length(band_center_output)))
    }
    
    spline_function <- stats::smooth.spline(center_to_spline, ref_to_spline, df = 40)
    
    res <- predict(spline_function, x = band_center_output)$y
    
    res[res < 0] <- 0
    
    return(res)
    
}, band_center_input = input_wvl, bad_bands_pos = input_bad_bands, band_center_output = output_wvl,
cores = 7, filename = paste0(output_folder, image_date, "_PRISMA_smoothed.tif")
) 


