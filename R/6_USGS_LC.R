#' Author: TJ Sipin
#' Date: April 7, 2024
#' Purpose:
#' Extract BCM data at each trap station-month-year.
#' Input: Presence-pseudoabsence points -- Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv
#' Output: Presence-pseudoabsence points with LC values -- Data/1_DataProcessing/df/USGS_p_psa_extract_2016_2023.csv

# Packages ----------------------------------------------------------------

library(readr)
library(dplyr)
library(raster)
library(sf)
library(ggplot2)
library(ggthemes)
library(viridis)
library(stringi)
library(tidyr)
library(purrr)

conflicted::conflict_prefer("select", winner = "dplyr")
conflicted::conflict_prefer("filter", winner = "dplyr")
conflicted::conflicts_prefer(raster::extract)



# Read in data ------------------------------------------------------------

p_psa = read_csv(" ")


# Create grid comprised with variable-year-month combinations
vars = c(
    "lc_type"
)
years = p_psa$Year %>% unique() %>% sort()

yv = expand.grid(var = vars, year = years, stringsAsFactors = F) 

### 2020 raster data ###
# ymv_2020 = which(ymv$year==2020)

# Extract function
extractFun = function(yyyy){
    t1 = Sys.time()
    # locate file name
    this_file_name = paste0("Data/1_DataProcessing/rasters/lc/Reproject_CONUS_A2_y", yyyy, "_.tif")
    # return NULL if file does not exist
    if(
        !file.exists(this_file_name)
    ){
        return(NULL)
    }
    
    print(this_file_name)
    # turn to raster
    this_raster = raster(this_file_name)
    
    
    # Get p_psa data at target month and year
    this_p_psa = p_psa %>% 
        filter(Year == yyyy)
    # flag for if there are no occurrence/absence points in month-year
    if(nrow(this_p_psa) == 0){
        return(NULL)
    }
    
    # Convert p_psa data to SPDF and transform to CRS of the raster data (NAD83)
    this_p_psa_spdf = this_p_psa %>% 
        SpatialPointsDataFrame(
            coords = data.frame(.$Longitude, .$Latitude),
            data = .,
            proj4string = CRS("EPSG:4326")
        ) %>% 
        spTransform(crs(this_raster))
    
    t2 = Sys.time()
    # Extract values at each point
    this_extract = extract(
        this_raster,
        this_p_psa_spdf
    )
    print(paste0(this_file_name, " extract: ", Sys.time() - t2))
    
    out = this_p_psa %>% 
        cbind(var = "lc_type", val = this_extract %>% as.factor()) 
    print(paste0(this_file_name, " total: ", Sys.time() - t1))
    
    return(out)
}

t0 = Sys.time()
extract_all = map(years, extractFun)
print(paste0("extract_all total: ", Sys.time() - t0))


# Pivot wider -------------------------------------------------------------

extract_all_df_long = extract_all %>% 
    bind_rows() 

extract_all_df = extract_all_df_long %>% 
    pivot_wider(
        names_from = "var",
        values_from = "val"
    ) %>% 
    mutate(lc_type = as.factor(lc_type))

write_csv(extract_all_df, "Data/1_DataProcessing/df/USGS_p_psa_extract_2016_2023.csv")
extract_all_df = read_csv("Data/1_DataProcessing/df/USGS_p_psa_extract_2016_2023.csv")

# t0 = Sys.time()
# extract_2020 = map(ymv_2020, extractFun)
# Sys.time() - t0
# 
# extract_2020_df_long = extract_2020 %>% 
#     bind_rows() 
# 
# extract_2020_df = extract_2020_df_long %>% 
#     pivot_wider(
#         names_from = "var",
#         values_from = "val"
#     )
# 
# dir.create("Data/1_DataProcessing/bcm/df/")
# write_csv(extract_2020_df, "Data/1_DataProcessing/bcm/df/p_psa_extract_2020.csv")
