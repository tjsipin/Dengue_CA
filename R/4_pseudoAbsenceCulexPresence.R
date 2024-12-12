#' Author: TJ Sipin
#' Date: September 19, 2024
#' Purpose:
#' Get pseudo-absence points for each month
#' Input: 
    #' Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv
    #' Data/1_DataProcessing/Culex/AedesQuinq_presencestatus.csv
#' Output: Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv
rm(list=ls())
set.seed(030824)
setwd("/home/tjsipin/network-storage/Dengue_CA")

# Packages ----------------------------------------------------------------

library(terra)
library(raster)
library(purrr)
library(dplyr)
library(flexsdm)
tidymodels::tidymodels_prefer()
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyterra)
library(readr)
library(tidyr)
library(sf)
library(stringr)
library(tigris)
# Data --------------------------------------------------------------------

# Declare directory to save into
output_dir = paste0(
    "Data/1_DataProcessing/pseudoabs/geo_env/"
)
output_filename = paste0(output_dir, "thinned_pr_bg_data_years.csv")

# Create directory
dir.create(output_dir, recursive=T, showWarnings = F)

# read in presence data
aaeg_p_spdf = read.csv("Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv") %>% 
    # Filter to species we want (A. aegypti)
    filter(species == "Aedes aegypti") %>% 
    # Species column for thinning function
    mutate(species = "aaeg") %>% 
    # Filter to years after initial introduction of A. aeg.
    filter(Year >= 2016) %>% 
    # Select only lat/long and month variables
    select(longitude, latitude, Month, Year, species) %>% 
    # set observation ID and presence values
    mutate(
        pr_ab = 1,
        type_id = row_number()
    ) %>% 
    # turn to spatial points data frame
    SpatialPointsDataFrame(
        coords = data.frame(.$longitude, .$latitude),
        data = .,
        proj4string = CRS("EPSG:4326")
    ) %>%
    spTransform(CRS("EPSG:3310"))

quinq_p_spdf = read.csv("Data/1_DataProcessing/Culex/AedesQuinq_presencestatus.csv") %>% 
    filter(pseudostatus == "OnlyQuinqs") %>% 
    # Species column for thinning function
    mutate(species = "quinq") %>% 
    # Filter to years after initial introduction of A. aeg.
    filter(year >= 2016) %>% 
    # Select only lat/long and month variables
    select(longitude, latitude, Month=month, Year=year, species) %>% 
    # set observation ID and presence values
    mutate(
        pr_ab = 0,
        type_id = row_number()
    ) %>% 
    # turn to spatial points data frame
    SpatialPointsDataFrame(
        coords = data.frame(.$longitude, .$latitude),
        data = .,
        proj4string = CRS("EPSG:4326")
    ) %>%
    spTransform(CRS("EPSG:3310"))

aaeg_quinx_p = rbind(aaeg_p_spdf, quinq_p_spdf) %>% 
    as.data.frame() %>% 
    rename(x = coords.x1, y = coords.x2, Longitude = longitude, Latitude = latitude) %>% 
    mutate(master_id = row_number())
# Set month abbs
months = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec') 

# Call in california shapefile
cali = maps::map('state', plot=F, fill=T, regions="california") %>%
    st_as_sf() %>%
    st_transform(CRS("EPSG:3310"))

cali_vect = terra::vect(cali)

# path to environmental rasters
env_raster_path = "Data/1_DataProcessing/rasters/bcm/monMean"
env_rasters = map(list.files(env_raster_path, full.names = T, pattern='.tif'), rast)
# Reorder the rasters by month
env_stack = rast(env_rasters) %>% 
    mask(cali_vect) %>% 
    crop(cali_vect)

# given month, get spatially thinned presences and write pseudo-absence locations as RDS


psuedoabsenceFunc = function(mon){
    tryCatch( 
        {
            # Progress check: print month
            print(mon)
            
            mm = str_pad(mon, width=2, pad='0')
            
            t0 = Sys.time() # get time of start
            
            # Declare output filename
            output_filename = paste0(output_dir, mm, "_thinned_presence_absence.rds")
            
            # Get environmental suitability raster for month
            this_env_rast = env_stack[[mon]] 
            
            # Subset presence points to target month 
            this_aaeg_quinq_p = aaeg_quinx_p %>% 
                filter(Month == mon) %>% 
                tibble() %>% 
                select(
                    Longitude, Latitude, x, y, Month, Year, pr_ab, species, master_id, type_id
                ) %>% 
                tibble()
            
            # Subset observed data by species
            this_aaeg_p = this_aaeg_quinq_p %>% 
                filter(species == "aaeg") %>% 
                mutate(row_id = row_number())
            this_quinq_p = this_aaeg_quinq_p %>% 
                filter(species == "quinq") %>% 
                mutate(row_id = row_number())
            
            
            t0 = Sys.time()
            this_aaeg_pf = spThin::thin(
                loc.data = this_aaeg_p, # location data
                lat.col = "Latitude",
                long.col = "Longitude",
                spec.col = "species",
                thin.par = 0.275, # distance (in km) for records to be separated by
                reps = 4, # number of times to repeat the thinning process
                write.files = F,
                locs.thinned.list.return = T,
                out.dir = paste0(this_dir, "thinned"),
                out.base = paste0(mon, "_pf")
            )
            Sys.time() - t0
            
            ## Choose the repetition with the most points
            this_aaeg_pf_best = this_aaeg_pf[[1]] %>% 
                mutate(row_id = rownames(.) %>% as.integer()) %>% 
                # mutate(thinned = T) %>% 
                left_join(this_aaeg_p) %>% 
                tibble() %>% 
                select(-row_id)
            
            t0 = Sys.time()
            this_quinq_pf = spThin::thin(
                loc.data = this_quinq_p, # location data
                lat.col = "Latitude",
                long.col = "Longitude",
                spec.col = "species",
                thin.par = 0.275, # distance (in km) for records to be separated by
                reps = 4, # number of times to repeat the thinning process
                write.files = F,
                locs.thinned.list.return = T,
                out.dir = paste0(this_dir, "thinned"),
                out.base = paste0(mon, "_pf")
            )
            Sys.time() - t0
            
            ## Choose the repetition with the most points
            this_quinq_pf_best = this_quinq_pf[[1]] %>% 
                mutate(row_id = rownames(.) %>% as.integer()) %>% 
                # mutate(thinned = T) %>% 
                left_join(this_quinq_p) %>% 
                tibble() %>% 
                select(-row_id)
            
            # Create spatvect out of presences for buffering step
            this_aaeg_pf_spatvect = this_aaeg_pf_best %>% 
                terra::vect(geom = c("Longitude", "Latitude"), "EPSG:4326")
            # Create spatvect out of pseudo-absences (quinq) to mask out points that are within the 275 meter presence buffers
            this_quinq_pf_spatvect = this_quinq_pf_best %>% 
                terra::vect(geom = c("Longitude", "Latitude"), "EPSG:4326")
            
            ### Knock out quinq points that are closer than 275 meters of presence points ###
            # Make buffer around presence points
            this_aaeg_pf_spatvect_buffer = buffer(this_aaeg_pf_spatvect, width = 275) %>% 
                terra::project(crs("EPSG:3310")) 
            # Make buffer around quinq points
            this_quinq_pf_spatvect_buffer = buffer(this_quinq_pf_spatvect, width = 275) %>% 
                terra::project(crs("EPSG:3310")) 
            
            # Find quinq points that intersect with presence point buffers
            this_quinq_pf_intersects = is.related(
                this_quinq_pf_spatvect_buffer, this_aaeg_pf_spatvect_buffer,
                relation = "intersects"
            )
            
            # Remove intersected points from quinx presence points
            this_quinq_pf_filtered = this_quinq_pf_spatvect[!this_quinq_pf_intersects]
            
            # Create 270 m buffer around quinx points
            this_quinq_pf_filtered_spatvect_buffer = buffer(this_quinq_pf_filtered, width = 275) %>% 
                terra::project(crs("EPSG:3310"))
            
            # Join the buffers with the environmental rasters for getting weights to sample points
            # We don't want points inside the buffers (spatial thin)
            this_env_rast_invbuff = this_env_rast %>% 
                classify(matrix(c(NA, NA, 1, 1, 0, 0), ncol = 2, byrow = T)) %>% 
                terra::mask(this_aaeg_pf_spatvect_buffer, inverse = T, updatevalue = -1) %>% 
                terra::mask(this_quinq_pf_filtered_spatvect_buffer, inverse = T, updatevalue = -2)
            
            this_quinq_pf_filtered_aaeg_pf = this_quinq_pf_filtered_spatvect_buffer %>% 
                rbind(this_aaeg_pf_spatvect_buffer) %>% 
                as.data.frame() %>% 
                tibble() %>% 
                SpatialPointsDataFrame(
                    coords = data.frame(.$x, .$y),
                    data = .,
                    proj4string = CRS("EPSG:3310")
                ) %>% 
                # Transform to NAD 83
                spTransform(CRS("EPSG:4326")) %>% 
                as.data.frame() %>% 
                rename(Longitude = coords.x1, Latitude = coords.x2) %>% 
                tibble()
            
            # Create raster sampling weights based on environmental suitability raster for the month
            # using 0.9 for unsuitable and 0.1 for suitable (checked in with Andy)
            # and 0 for cells in the buffers around the presence points
            this_env_rast_probabilities = this_env_rast_invbuff %>% 
                classify(matrix(c(-2, 0, -1, 0, 0, 0.9, 1, 0.1), ncol = 2, byrow = T))
            
            # this_env_rast_probabilities %>% plot(main="June env. suitability")
            # ggplot() +
            #     geom_sf(data=ca) +
            #     geom_spatraster(data=this_env_rast_probabilities) #+
            # geom_spatvector(data=this_aaeg_p_spatvect_buffer, color='blue') +
            # geom_spatvector(data=this_quinq_p_spatvect[this_quinq_p_intersects], color='red', alpha = 1, size = 0.2) +
            # geom_spatvector(data=this_quinq_p_filtered_spatvect_buffer, color='orange', alpha = 1, size = 0.2) +
            # xlim(c(-120, -119.2)) +
            # ylim(c(36.2, 37))
            
            #this_env_rast_probabilities %>% plot()
            t0=Sys.time()
            
            this_aaeg_bg = spatSample(
                x = this_env_rast_probabilities,
                size = nrow(this_aaeg_pf_best) + 1000,
                method = 'weights',
                xy = T, 
                as.points = F
            ) %>% 
                mutate(
                    pr_ab = 0,
                    Month = mon,
                    species = "aaeg"
                ) %>% 
                select(x, y, pr_ab, Month, species) %>%
                SpatialPointsDataFrame(
                    coords = data.frame(.$x, .$y),
                    data = .,
                    proj4string = CRS("EPSG:3310")
                ) %>% 
                # Transform to NAD 83
                spTransform(CRS("EPSG:4326")) %>% 
                as.data.frame() %>% 
                rename(Longitude = coords.x1, Latitude = coords.x2) %>% 
                tibble() %>% 
                mutate(type_id = row_number(), master_id = NA, Year = NA)
            
            Sys.time() - t0
            
            # Bind quinx/aaeg data with bg data
            this_aaeg_quinq_prbg = this_quinq_pf_filtered_aaeg_pf %>% 
                # select(-row_number) %>%
                rbind(this_aaeg_bg) %>% 
                as.data.frame() %>% 
                # rename(Longitude = coords.x1, Latitude = coords.x2) %>% 
                mutate(row_id = rownames(.) %>% as.integer()) %>%
                tibble() %>% 
                mutate(species_standin = "na")
            
            # Start thinning process for A. aeg. presence points
            t0 = Sys.time()
            this_aaeg_quinq_prbg_pf = spThin::thin(
                loc.data = this_aaeg_quinq_prbg, # location data
                lat.col = "Latitude",
                long.col = "Longitude",
                spec.col = "species_standin",
                thin.par = 0.275, # distance (in km) for records to be separated by
                reps = 4, # number of times to repeat the thinning process
                write.files = F,
                locs.thinned.list.return = T,
                out.dir = paste0(this_dir, "thinned"),
                out.base = paste0(mon, "_pf")
            )
            Sys.time() - t0
            
            this_aaeg_quinq_prbg_pf_best = this_aaeg_quinq_prbg_pf[[1]] %>% 
                # tibble() %>% 
                mutate(row_id = rownames(.) %>% as.integer()) %>% 
                # mutate(thinned = T) %>% 
                left_join(this_aaeg_quinq_prbg) %>% 
                select(-species_standin) %>% 
                tibble() %>% 
                select(-row_id)
            
            
            
            ### Assign Year to spatially sampled PSA points
            this_year_weights = this_aaeg_p %>% 
                group_by(Year) %>% 
                summarize(prop = n()/nrow(.)) %>% 
                ungroup() %>% 
                arrange(Year) 
            
            na_size = this_aaeg_quinq_prbg_pf_best %>%
                filter(is.na(Year)) %>%
                nrow()
            
            na_years = sample(this_year_weights$Year, size = na_size, prob = this_year_weights$prop, replace = T)
            
            this_aaeg_quinq_prbg_pf_best_nonnayears = this_aaeg_quinq_prbg_pf_best %>% 
                filter(!is.na(Year))
            
            this_aaeg_quinq_prbg_pf_best_nayears = this_aaeg_quinq_prbg_pf_best %>% 
                filter(is.na(Year))
            # set Year to sampled years based on frequency
            this_aaeg_quinq_prbg_pf_best_nayears$Year = na_years
            
            
            this_aaeg_quinq_prbg_pf_best_years =  this_aaeg_quinq_prbg_pf_best_nonnayears %>% 
                rbind(this_aaeg_quinq_prbg_pf_best_nayears)
            
            saveRDS(this_aaeg_quinq_prbg_pf_best_years, output_filename)
            
            message(paste0(mon, " has finished."))
            message(paste0(mon, ": ", t0 - Sys.time()))
            
        },
        error = function(cond){
            message(paste0(mon, " has failed. See month directory for more details."))
            
            message(cond)
        }
    )
}

### To run in parallel, uncomment these lines of code
# registerDoParallel(detectCores() - 1)
# foreach(m = 1:12)%dopar%{
#     psuedoabsenceFunc(m)
# }
# stopImplicitCluster()

### To not run in parallel
for(m in 1:12){
    psuedoabsenceFunc(m)
}


# get all pseudo-absence files and reorder
# As of August 8, only able to get pseudo-absences for March to November
pseudoabs_files = list.files(
    path = output_dir,
    recursive = T, pattern = "thinned_presence_absence.rds", full.names = T
)

p_bg_records = map(pseudoabs_files, readRDS) %>%
    map(bind_rows) %>%
    bind_rows() %>% 
    mutate(pr_ab = as.factor(pr_ab)) %>% 
    mutate(
        type = case_when(
            (pr_ab==1) & (species=="aaeg") ~ "Presence",
            (pr_ab==0) & (species=="quinq") ~ "PSA (quinx)",
            (pr_ab==0) & (species=="aaeg") ~ "PSA (function)"
        )
    )

# Output combined presence/pseudo-absence data frame
write_csv(p_bg_records, output_filename)