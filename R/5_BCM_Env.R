#' Author: TJ Sipin
#' Date: April 7, 2024
#' Purpose:
#' Extract BCM data at each trap station-month-year.
#' Input: Presence-pseudoabsence points -- Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv
#' Output: Presence-pseudoabsence points with BCM values -- Data/1_DataProcessing/df/BCM_p_psa_extract_2016_2023.csv

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

p_psa = read_csv("Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv")


# Create grid comprised with variable-year-month combinations
vars = c(
    'aet', 'cwd', 'pet', 'ppt', 'tmn', 'tmx'
)
years = 2016:2023
months_names = c(
    'jan', 'feb', 'mar', 'apr', 'may', 'jun',
    'jul', 'aug', 'sep', 'oct', 'nov', 'dec'
)
months_nums = data.frame(
    mon_name = months_names,
    mon_num = 1:12
)

ymv = expand.grid(var = vars, year = years, mon_name = months_names, stringsAsFactors = F) %>% 
    full_join(months_nums)

# Extract function
extractFun = function(ymv_ind){
    t1 = Sys.time()
    # locate file name
    this_file_name = paste0("Data/1_DataProcessing/rasters/bcm/", paste0(ymv[ymv_ind, 1:3], collapse=""), ".tif")
    # return NULL if file does not exist
    if(
        !file.exists(this_file_name)
    ){
        return(NULL)
    }
    
    print(this_file_name)
    # turn to raster
    this_raster = raster(this_file_name)
    
    # Get variable, month, and year names
    var = paste0(ymv[ymv_ind, 1], recycle0 = T)
    month = ymv[ymv_ind, 4]
    year = ymv[ymv_ind, 2]
    
    # Get p_psa data at target month and year
    this_p_psa = p_psa %>% 
        filter(Year == year, Month == month)
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
        cbind(var = var, val = this_extract) 
    print(paste0(this_file_name, " total: ", Sys.time() - t1))
    
    return(out)
}

t0 = Sys.time()
extract_all = map(1:nrow(ymv), extractFun)
print(paste0("extract_all total: ", Sys.time() - t0))


# Pivot wider -------------------------------------------------------------

extract_all_df_long = extract_all %>% 
    bind_rows() 

extract_all_df = extract_all_df_long %>% 
    pivot_wider(
        names_from = "var",
        values_from = "val"
    )

write_csv(extract_all_df, "Data/1_DataProcessing/df/BCM_p_psa_extract_2016_2023.csv")
extract_all_df = read_csv("Data/1_DataProcessing/df/BCM_p_psa_extract_2016_2023.csv")