#' Author: TJ Sipin
#' Date: March 20, 2024
#' Purpose:
#' Get absence and presence point locations by month-year,
#' such that there is an equal number of presence and absence points per month-year.
#' Include presence points of culex as pseudo-absence points
#' Steps:
#' Filter presence points through thinning process (270 m) to create data set A.
#' Create pseduo-absence points from thinned presence points
#' Input: Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv
rm(list=ls())
set.seed(032024)
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
library(ggmap)

# Data --------------------------------------------------------------------

p_bg_records = read_csv("Data/1_DataProcessing/pseudoabs/geo_env/thinned_pr_bg_data_years.csv") 

p_bg_records %>% 
    group_by(Month) %>% 
    summarize(
        Presence = sum(type=="Presence")/n(),
        PSA_function = sum(type=="PSA (function)")/n(),
        PSA_quinx = sum(type=="PSA (quinx)")/n()
    ) %>% 
    ungroup()

p_bg_records_vect = p_bg_records %>% 
    vect(geom=c('x', 'y'), crs="EPSG:3310")

# Access CA counties
ca = counties("California")

# Access CA state shapefile
cali = maps::map('state', plot=F, fill=T, regions="california") %>%
    st_as_sf() %>%
    # st_transform("+proj=longlat +datum=WGS84")
    st_transform(CRS("EPSG:3310"))
cali_vect = terra::vect(cali)

# path to environmental rasters
env_raster_path = "Data/1_DataProcessing/rasters/bcm/monMean"
env_rasters = map(list.files(env_raster_path, full.names = T, pattern='.tif'), rast)
# Reorder the rasters by month and mask/crop to CA borders
env_stack = rast(env_rasters) %>% 
    mask(cali_vect) %>% 
    crop(cali_vect)

# Plot presence/background points by month
ggplot() +
    geom_sf(data = ca) +
    geom_spatvector(
        data=p_bg_records_vect %>% 
            filter(Month %in% 1:12), 
        aes(color = type),
        alpha = 0.25, 
        cex = 0.3
    ) +
    facet_grid(Year~Month) +
    labs(
        title = "Presence/background by month-year"
    ) +
    theme_classic()

# Plot presence/background points by month
ggplot() +
    geom_sf(data = ca) +
    geom_spatvector(
        data=p_bg_records_vect %>% 
            filter(Month %in% 1:12), 
        aes(color = type),
        alpha = 0.25, 
        cex = 0.3
    ) +
    facet_wrap(~Month) +
    labs(
        title = "Presence/background by month"
    ) +
    theme_classic()

# Plot env suit rasters by month
ggplot() + 
    geom_spatraster(
        data = env_stack
    ) +
    facet_wrap(~lyr) +
    labs(
        title = "Environmental suitability by month",
        subtitle = "(lower bound min temp CI, upper bound max temp CI)"
    ) + 
    theme_classic()

p_bg_records %>% 
    mutate(month_year = 12*(Year - 2016) + Month) %>% 
    group_by(month_year) %>% 
    summarize(
        n_presence = sum(type=="Presence"),
        n_quinx = sum(type=="PSA (quinx)"),
        n_psa = sum(type=="PSA (function)")
    ) %>% 
    ungroup() %>% 
    pivot_longer(-month_year) %>% 
    rename(Observation_type = name, n = value) %>% 
    mutate(Year = 2016 + month_year/12) %>% 
    ggplot() + 
    geom_line(
        aes(
            x = Year,
            y = n, 
            color = Observation_type
        )
    )


aaeg_p_spdf %>% 
    as.data.frame() %>% 
    group_by(Year, Month) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    mutate(Year = ordered(Year)) %>% 
    ggplot() + 
    geom_line(aes(x = Month, y = n, color = Year)) +
    labs(title = "Unfiltered presence points by month and year")

p_bg_records %>% 
    as.data.frame() %>% 
    group_by(Year, Month, type) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    mutate(Year = ordered(Year)) %>% 
    ggplot() + 
    geom_line(aes(x = Month, y = n, color = Year)) +
    labs(title = "Unfiltered points by month and year") +
    facet_wrap( ~ type)


# Summarize the number of presence and pseudo-absence points (per unit area perhaps as well) 
# that fall within suitable vs. unsuitable regions, respectively? 
# Just to give us a sense of how well sampled suitable vs. unsuitable regions are?

countPointsByRegion = function(mon){
    env_rast = env_stack[[mon]] 
    env_rast_suitable = ifel(env_rast==1, 1, NA)
    env_rast_unsuitable = ifel(env_rast==0, 1, NA)
    
    pr = p_bg_records_vect %>% 
        filter(
            Month==mon,
            pr_ab==1
        )
    
    bg = p_bg_records_vect %>% 
        filter(
            Month==mon,
            pr_ab==0,
            species=="aaeg"
        )
    
    pr_extract = terra::extract(
        x=env_rast, 
        y=pr
    )
    names(pr_extract) = c("ID", "suitable")
    pr_count_1 = pr_extract %>% 
        filter(suitable==1) %>% 
        nrow()
    
    pr_count_0 = pr_extract %>% 
        filter(suitable==0) %>% 
        nrow()
    
    bg_extract = terra::extract(
        x=env_rast, 
        y=bg
    )
    names(bg_extract) = c("ID", "suitable")
    bg_count_1 = bg_extract %>% 
        filter(suitable==1) %>% 
        nrow()
    
    bg_count_0 = bg_extract %>% 
        filter(suitable==0) %>% 
        nrow()
    
    env_rast_f = freq(env_rast)
    env_rast_f$area = env_rast_f$count * 0.27*0.27 # in km^2
    
    tibble(
        month = mon,
        observation_type = c("Presence", "Presence", "Absence", "Absence"),
        count = c(pr_count_1, pr_count_0, bg_count_1, bg_count_0),
        area = c(env_rast_f$area[2], env_rast_f$area[1], env_rast_f$area[2], env_rast_f$area[1]),
        suitability = rep(c("suitable", "unsuitable"), 2)
    ) %>% 
        mutate(
            count_per_km2 = count/area
        ) #%>% 
    #pivot_wider(names_from = "suitability", values_from = c("area", "count_per_km2"))
    
}

count_points_by_region = map(1:12, countPointsByRegion) %>% 
    bind_rows()

checkCloseN = function(mon){
    print(mon)
    this_month_data = p_bg_records_vect %>%
        filter(Month == mon)
    
    dist.mat <- this_month_data %>% 
        distance()
    # SpatialPointsDataFrame(
    #     data = .,
    #     coords = data.frame(.$Longitude, .$Latitude),
    #     proj4string = CRS("EPSG:4326")
    # ) %>% 
    # st_as_sf() %>% 
    # st_distance() # Great Circle distance since in lat/lon
    # Number within 270m: Subtract 1 to exclude the point itself
    num.270 <- apply(dist.mat, 1, function(x) {
        sum(x <= 270) - 1
    })
    
    neighbors.index <- apply(dist.mat, 1, function(x){
        as.list(which(x <= 270))
    })
    
    
    # Calculate nearest distance
    nn.dist <- apply(dist.mat, 1, function(x) {
        return(sort(x, partial = 2)[2])
    })
    # Get index for nearest distance
    nn.index <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[2] })
    
    # neighbors.index <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[2:21] })
    
    
    n.data <- this_month_data %>% 
        as.data.frame() %>% 
        tibble()
    colnames(n.data)[1] <- "neighbor"
    colnames(n.data)[2:ncol(n.data)] <- 
        paste0("n.", colnames(n.data)[2:ncol(n.data)])
    mydata2 <- data.frame(this_month_data,
                          n.data[nn.index, ],
                          n.distance = nn.dist,
                          radius270 = num.270)
    rownames(mydata2) <- seq(nrow(mydata2))
    
    mydata2
    # dist_data = dist_data %>% 
    #   full_join(mydata2)
}

t0=Sys.time()
check_close_n = map(1:12, checkCloseN) %>% 
    bind_rows() %>% 
    tibble()
Sys.time()-t0

check_close_n %>% 
    arrange(n.distance) %>% 
    select(n.distance, radius270, everything())
