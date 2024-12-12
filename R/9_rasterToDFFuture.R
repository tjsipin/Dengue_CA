library(tidyverse)
library(tidymodels)
library(ranger)
library(terra)
library(sf)
conflicted::conflicts_prefer(terra::resample)
conflicted::conflicts_prefer(dplyr::filter)

months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)

mid_mon_years = expand.grid(
    months = ordered(months),
    years = 2040:2050
) 

end_mon_years = expand.grid(
    months = ordered(months),
    years = 2090:2099
)

# Get CA polygon for masking and cropping USGS data 
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs=crs("EPSG:3310")) # for consistency with BCM data

# Create function to take BCM and USGS LC rasters and output data frames with values per pixel
rastToDF = function(month, year, RCP){
    # Create date string
    this_date_str = lubridate::my(paste0(month, year))
    print(this_date_str)
    
    # Define output filenames
    outputrast_filename = paste0("Data/1_DataProcessing/rasters/future/", RCP, "/future_", this_date_str, "_", RCP, ".tif")
    outputdf_filename = paste0("Data/1_DataProcessing/df/future/", RCP, "/future_", this_date_str, "_", RCP, ".rds")
    
    # If output files already exist, then skip
    if(file.exists(outputrast_filename) & file.exists(outputdf_filename)) return(NULL)
    
    if(file.exists(outputrast_filename)) {
        outputrast = rast(outputrast_filename)
    } else {
        # Read in BCM data
        bcmrast_filename = paste0("Data/1_DataProcessing/rasters/bcm/", RCP, "/CA_BCM_MIROC_rcp", RCP, "_", this_date_str, ".tif")
        if(!file.exists(bcmrast_filename)) return(NULL)
        # Mask and crop BCM data to CA
        bcmrast = rast(bcmrast_filename) %>% 
            terra::mask(CA) %>%
            terra::crop(CA)
        
        # Read in usgs data
        usgsrast_filename_resampled = paste0("Data/1_DataProcessing/rasters/future/lc/CONUS_A2_y", year, "_resampled.tif")
        usgsrast_filename = paste0("Data/usgs/CONUS_A2_y", year, ".tif")
        
        # If resampled data does not yet exist, create it. 
        # Otherwise, read it in
        if(!file.exists(usgsrast_filename_resampled)){
            # Project and resample to BCM raster
            # Mask and crop to CA
            usgsrast = rast(usgsrast_filename) %>% 
                terra::project("WGS84", method="near") %>% 
                terra::project(bcmrast, method="near") %>% 
                terra::mask(CA) %>% 
                terra::crop(CA) %>% 
                terra::resample(bcmrast, method="near")
            names(usgsrast) = "lc_type"
            
            writeRaster(usgsrast, usgsrast_filename_resampled)
        } else{
            usgsrast = rast(usgsrast_filename_resampled)
        }
        
        # Combine BCM and USGS rasters
        outputrast = rast(list(bcmrast, usgsrast))
        
        # Remove USGS raster from global environment to clear up memory usage
        rm(usgsrast)
        
        writeRaster(outputrast, outputrast_filename)
    }
    
    # Get DF from raster
    outputdf = outputrast %>% 
        as.data.frame(xy=T) %>% 
        # Only want complete cases
        dplyr::filter(complete.cases(.))
    
    saveRDS(outputdf, outputdf_filename)
    
    rm(outputrast, outputdf, bcmrast)
}

month_year_grid = expand.grid(
    months = str_pad(1:12, width = 2, side = "left", pad = "0"),
    years = c(2040:2050, 2090:2099)
)
map2(
    .x = month_year_grid$months,
    .y = month_year_grid$years, 
    ~ rastToDF(month = .x, year = .y, RCP = "45")
)

map2(
    .x = month_year_grid$months,
    .y = month_year_grid$years,
    ~ rastToDF(month = .x, year = .y, RCP = "60")
)

map2(
    .x = month_year_grid$months,
    .y = month_year_grid$years,
    ~ rastToDF(month = .x, year = .y, RCP = "85")
)
