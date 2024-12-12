library(tidyverse)
library(tidymodels)
library(ranger)
library(terra)
library(sf)
library(sp)
conflicted::conflicts_prefer(terra::resample)
conflicted::conflicts_prefer(dplyr::filter)

months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)
mon_years = expand.grid(
    months = ordered(months),
    years = 2010:2024
)

# Get CA polygon for masking and cropping USGS data 
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs=CRS("EPSG:3310")) # for consistency with BCM data


# Create function to take BCM and USGS LC rasters and output data frames with values per pixel
rasterToDF = function(mon, year){
    data_filename = paste0("Data/1_DataProcessing/df/env_lc/cleaned_env_", mon, year, ".rds")
    print(paste(mon, year))
    # If the cleaned data frame does not exist, make it
    if(!file.exists(data_filename)){
        # List out BCM raster files
        bcm_files = list.files(
            path="Data/1_DataProcessing/rasters/bcm",
            pattern=paste0(year, mon, ".tif"),
            full.names=T
        )
        
        # Create raster out of files, masking/cropping within CA borders
        bcm_rast = rast(bcm_files) %>% 
            terra::mask(CA) %>% 
            terra::crop(CA)
        
        # Make names only environmental variables (remove year-month)
        names(bcm_rast) = names(bcm_rast) %>% str_split_i(pattern="2", i=1)
        
        # Read in USGS data
        usgs_filename = paste0("Data/1_DataProcessing/rasters/lc_v2/Resample_CONUS_A2_y", year, ".tif")
        
        # If resampled USGS data raster does not exist, make it
        if(!file.exists(usgs_filename)){
            # USGS file name raw/downloaded from Google Drive
            usgs_filename_v0 = paste0("Data/1_DataProcessing/rasters/lc/Reproject_CONUS_A2_y", year, "_.tif")
            
            t0=Sys.time()
            
            # Read in raw file name
            usgs_rast =  rast(usgs_filename_v0) %>% 
                # Project first to WGS84
                terra::project("EPSG:4326", method="near") %>% 
                # Then project to the BCM raster above (NAD83)
                terra::project(bcm_rast, method="near") %>%
                # Then mask and crop to CA
                terra::mask(CA) %>% 
                terra::crop(CA) %>% 
                # Finally resample to BCM raster above to get to that resolution
                terra::resample(bcm_rast, method="near")
            
            Sys.time()-t0
            
            names(usgs_rast) = "lc_type"
            writeRaster(usgs_rast, usgs_filename)
        } else{
            usgs_rast = rast(usgs_filename)
        }
        
        # Combine BCM and USGS LC rasters
        master_rast = rast(list(bcm_rast, usgs_rast))
        
        # Make data frame out of raster
        master_df = master_rast %>% 
            as.data.frame(xy=T) %>% 
            # Create environmental variables to be predicted on
            dplyr::mutate(
                mean_temp = (tmn + tmx)/2,
                tws = ifelse(pet == 0, 0, aet/pet)
            ) %>% 
            # Remove all rows with NA
            na.omit()
        
        saveRDS(master_df, data_filename)
        
        rm(master_rast, master_df, usgs_rast)
        
    } else{
        return(NULL)
    }
}

map2(
    .x = mon_years$months,
    .y = mon_years$years,
    ~ rasterToDF(.x, .y)
)
