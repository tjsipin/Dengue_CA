#' Author: TJ Sipin
#' Date: Feb 9, 2024
#' Purpose: Create occurrence points data set of mixed female Aedes across CA.
#' Input: Data/1_DataProcessing/Aedes/aedes_slim_mixedfemales.csv
#' Output: 
    #' Data/1_DataProcessing/Aedes/aedes_abundNA_mixedfemales.csv
    #' Data/1_DataProcessing/Aedes/aedes_abundNAzero_mixedfemales.csv
    #' Data/1_DataProcessing/Aedes/aedes_presenceAbsence_mixedfemales.csv
    #' Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv

library(dplyr)
library(tidyverse)

rm(list = ls())
output_dir = "Data/1_DataProcessing/Aedes/"
input_data = read_csv("Data/1_DataProcessing/Aedes/aedes_slim_mixedfemales.csv")

trapID_dict = input_data %>% 
    select(calculated_city, calculated_county, site_code) %>% 
    distinct() %>% 
    mutate(trapID = row_number())

joinTrapID = input_data %>% 
    full_join(trapID_dict)

# Obtain mosPerTrapNight for abundance analysis
abund_NA = joinTrapID %>% 
    mutate(mosPerTrapNight = females_mixed/trap_nights) %>% 
    ungroup()

# Coerce NA mosPerTrapNight to 0 
abund_NA_zero = abund_NA %>% 
    mutate(mosPerTrapNight = ifelse(is.na(mosPerTrapNight), 0, mosPerTrapNight)) 

# Create presence/absence points data frame
presence_absence_df = joinTrapID %>% 
    group_by(Year, Month, woy, longitude, latitude, trapID, calculated_city, calculated_county, site_code, collection_id, species) %>% 
    summarise(females_mixed = sum(females_mixed, na.rm = T) > 0) %>% 
    ungroup()

# Create presence points data frame
presence_df = presence_absence_df %>% 
    filter(females_mixed)

# boxplot mosPerTrapNight by county
ggplot(
    abund_NA_zero, 
    aes(
        x = log(mosPerTrapNight+1), 
        y = fct_reorder(calculated_county %>% as.factor(), log(mosPerTrapNight+ 1), median, na.rm = T, .desc = F)
    )
) + 
    geom_boxplot() + 
    facet_grid(rows = "species") +
    ggthemes::theme_solarized_2() + 
    ylab("County") + 
    xlab("log(mosquitoes per trap night)")

# plot mean mosPerTrapNight per year by county 
ggplot(
    abund_NA_zero %>% 
        group_by(Year, calculated_county)
) + 
    geom_smooth(aes(x = Year, y = log(mosPerTrapNight + 1)), method= "lm") + 
    facet_wrap(facets = "calculated_county")


write_csv(abund_NA, str_c(output_dir, "aedes_abundNA_mixedfemales.csv"))
write_csv(abund_NA_zero, str_c(output_dir, "aedes_abundNAzero_mixedfemales.csv"))
write_csv(presence_absence_df, str_c(output_dir, "aedes_presenceAbsence_mixedfemales.csv"))
write_csv(presence_df, str_c(output_dir, "aedes_presence_mixedfemales.csv"))
