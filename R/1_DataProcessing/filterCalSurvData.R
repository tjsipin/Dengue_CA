#' Author: TJ Sipin
#' Date: Feb 9, 2024
#' Purpose: Filter CalSurv data in similar ways to how Sam had, but include site_code as an ID variable
#' Input: Data/0_Raw/CalSurv_AedesSpecies.csv
#' Output: Data/1_DataProcessing/Aedes/aedes_slim_mixedfemales_v3.csv
library(dplyr)
library(lubridate)
library(tidyverse)

output_dir = "Data/1_DataProcessing/Aedes/"

raw = read_csv("Data/0_Raw/CalSurv_AedesSpecies.csv")

# filter to observations in state
stateFilter = raw %>% 
  filter(calculated_state == "California")

# remove observations with trap problems;
# those with trap problems show NA in count data
trapFilter = stateFilter %>% 
  filter(trap_problem == "N")

# create date ID columns
dateCols = trapFilter %>% 
  mutate(
    collection_date = mdy(collection_date), 
    Year = year(collection_date),
    Month = month(collection_date),
    woy = week(collection_date)
  )

# select final variables that match with Sam's aedes_slim_mixedfemales_v2.csv
finalOutput = dateCols %>% 
  select(
    collection_date, Year, Month, woy, collection_id,
    site_code, calculated_county, calculated_city, 
    longitude, latitude, trap_type, num_trap, trap_nights, trap_problem, 
    species, females_mixed = `females - mixed`
  )

write_csv(finalOutput, str_c(output_dir, "aedes_slim_mixedfemales_v3.csv"))
