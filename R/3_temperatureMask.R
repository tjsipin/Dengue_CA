#' Author: TJ Sipin 
#' Date: 09/27/24
#' Purpose:
    #' Instead of using point-estimates of minimum and maximum temperature suitability, 
    #' use the lower bound of the 95% CI for minimum temperature (6.9)
    #' and the upper bound of the 95% CI for maximum temperature (39.7)
    #' to give a more conservative approach
#' Input: temperature rasters found in Data/1_DataProcessing/rasters/bcm/
#' Output:

library(terra)
library(dplyr)

tmin_filenames = list.files(
    "Data/1_DataProcessing/rasters/bcm",
    pattern=paste(paste0("tmn", 2010:2020), collapse = "|")
)

tmax_filenames = list.files(
    "Data/1_DataProcessing/rasters/bcm",
    pattern=paste(paste0("tmx", 2010:2020), collapse = "|")
)

# Set month abbs
months = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
mm = str_pad(1:12, width = 2, "0", side = "left")

monthlyTempFunc = function(mon, mm){
    # Declare output tif file
    output_filename = paste0("Data/1_DataProcessing/rasters/bcm/monMean/", mm, "_", mon, "_temperature_mask.tif")
    message(output_filename)
    
    # Get tmn files
    tmn_filenames = list.files(
        "Data/1_DataProcessing/rasters/bcm",
        pattern=paste(paste0("tmn", 2010:2020, mon), collapse = "|"),
        full.names = T
    )
    # Get tmx files
    tmx_filenames = list.files(
        "Data/1_DataProcessing/rasters/bcm",
        pattern=paste(paste0("tmx", 2010:2020, mon), collapse = "|"),
        full.names = T
    )
    # Get tmn rast
    tmn_rast = rast(tmn_filenames) %>% 
        mean(na.rm = T) 
    tmn_rast = ifel(tmn_rast < 6.9, 0, 1)
    names(tmn_rast) = "tmn"
    
    # Get tmx rast
    tmx_rast = rast(tmx_filenames) %>% 
        mean(na.rm = T) 
    tmx_rast = ifel(tmx_rast > 39.7, 0, 1)
    names(tmx_rast) = "tmx"
    
    trange_rast = rast(list(tmn_rast, tmx_rast)) %>% 
        sum(na.rm = T)
    
    trange_rast = ifel(trange_rast == 2, 1, 0)
    names(trange_rast) = paste0(mon, "_env_suit")
    
    writeRaster(trange_rast, output_filename, overwrite = T)
}

map2(.x = months, .y = mm, monthlyTempFunc)