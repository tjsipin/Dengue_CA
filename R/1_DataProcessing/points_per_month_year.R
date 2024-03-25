#' Author: TJ Sipin
#' Date: March 20, 2024
#' Purpose:
#' Get absence and presence point locations by month-year,
#' such that there is an equal number of presence and absence points per month-year

set.seed(032024)
setwd("/home/tjsipin/Dengue_CA")

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

# Data --------------------------------------------------------------------

thinned_presence_pabsence_data = read_csv("Data/1_DataProcessing/pseudoabs/geo_env/aedes_thinned_presence_psa.csv") 

thinned_presence_data = thinned_presence_pabsence_data %>% 
  filter(pr_ab == 1) %>% 
  mutate(round_Long = round(Longitude, 4),
         round_Lat = round(Latitude, 4))

thinned_presence_data_dates = read_csv("Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv") %>% 
  # turn to spatial points data frame
  SpatialPointsDataFrame(
    coords = data.frame(.$longitude, .$latitude),
    data = .,
    proj4string = CRS("+proj=longlat +datum=WGS84")
  ) %>% 
  spTransform(CRS("EPSG:3310")) %>% 
  as.data.frame() %>% 
  mutate(round_Long = round(longitude, 4),
         round_Lat = round(latitude, 4)) %>% 
  right_join(thinned_presence_data, by=c("round_Long", "round_Lat", "Month")) %>%
  # right_join(thinned_presence_data, by = c("longitude"="Longitude", "latitude"="Latitude", "Month")) %>% 
  select(Year, Month, Longitude, Latitude, x, y) %>%
  # select(Year, Month, x, y) %>% 
  distinct() %>% 
  na.omit() %>% 
  tibble() %>% 
  mutate(pr_ab = 1)

presence_data_n = thinned_presence_data_dates %>% 
  group_by(Month, Year) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  arrange(Year, Month)

pabsence_data = thinned_presence_pabsence_data %>% 
  filter(pr_ab == 0) 

pabsence_data_dates = data.frame(
  Month = integer(0L),
  Year = integer(0L),
  Longitude = double(0L),
  Latitude = double(0L),
  x = double(0L),
  y = double(0L),
  pr_ab = integer(0L)
)

for(i in 1:nrow(presence_data_n)){
  month = presence_data_n$Month[i]
  year = presence_data_n$Year[i]
  n = presence_data_n$n[i]
  
  this_pabs_data = pabsence_data %>% 
    filter(Month == month) %>% 
    slice_sample(n = n) %>% 
    mutate(Year = year, pr_ab = 0) %>% 
    select(Month, Year, Longitude, Latitude, x, y, pr_ab)
  
  pabsence_data_dates = pabsence_data_dates %>% 
    rbind(this_pabs_data)
}

p_psa_month_year = rbind(thinned_presence_data_dates, pabsence_data_dates) 

p_psa_month_year_sf = p_psa_month_year %>% 
  SpatialPointsDataFrame(
    coords = data.frame(.$Longitude, .$Latitude),
    data = .,
    proj4string = CRS("+proj=longlat +datum=WGS84")
  ) %>% 
  st_as_sf()

write_csv(p_psa_month_year, "Data/1_DataProcessing/pseudoabs/geo_env/thinned_p_psa_month_year.csv")

write_sf(p_psa_month_year_sf, "Data/1_DataProcessing/pseudoabs/geo_env/thinned_p_psa_month_year_sf.shp")
