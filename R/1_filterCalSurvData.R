#' Author: TJ Sipin
#' Date: Feb 9, 2024
#' Purpose: Filter CalSurv data in similar ways to how Sam had, but include site_code as an ID variable
#' Input: Data/0_Raw/CalSurv_AedesSpecies.csv
#' Output: Data/1_DataProcessing/Aedes/aedes_slim_mixedfemales.csv

rm(list = ls())

library(dplyr)
library(lubridate)
library(tidyverse)


output_dir = "Data/1_DataProcessing/Aedes/"
dir.create(output_dir)

raw = read.csv("Data/0_Input/CalSurv_AedesSpecies.csv")

# filter to just aegypti
speciesFilter = raw %>% 
    filter(species == "Aedes aegypti") %>% 
    mutate(females.mixed_perTrpNt = females...mixed/trap_nights)

# remove observations with trap problems;
# those with trap problems show NA in count data
trapFilter = speciesFilter %>% 
  filter(trap_problem == "N")

# filter to observations in state
stateFilter = trapFilter %>% 
    filter(calculated_state == "California")

nonNAFilter = stateFilter %>% 
    filter(!is.na(females...mixed)) 

posMosFilter = nonNAFilter %>% 
    filter(females...mixed > 0)

# create date ID columns
dateCols = posMosFilter %>% 
    mutate(
        collection_date = mdy(collection_date), 
        Year = year(collection_date),
        Month = month(collection_date),
        woy = week(collection_date)
    ) %>% 
    # select final variables that match with Sam's aedes_slim_mixedfemales_v2.csv
    select(
        collection_date, Year, Month, woy, collection_id,
        site_code, calculated_county, calculated_city, 
        longitude, latitude, trap_type, num_trap, trap_nights, trap_problem, 
        species, females_mixed = females...mixed, coordinate_precision
    )

detectedCoordsFilter = dateCols %>% 
    filter(coordinate_precision != "Unknown")


write_csv(dateCols, str_c(output_dir, "aedes_slim_mixedfemales.csv"))
write_csv(detectedCoordsFilter, str_c(output_dir, "aedes_slim_mixedfemales_coordsfilter.csv"))

