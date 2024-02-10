#' Author: TJ Sipin
#' Date: Feb 9, 2024
#' Purpose: Create occurrence points data set of mixed female Aedes across CA.
#' Input: 
#' 

library(dplyr)
library(tidyverse)

output_dir = "Data/1_DataProcessing/Aedes/"
input_data = read_csv("Data/1_DataProcessing/Aedes/aedes_slim_mixedfemales_v3.csv") 

trapID_dict = input_data %>% 
  select(calculated_city, calculated_county, site_code) %>% 
  distinct() %>% 
  mutate(trapID = row_number())

joinTrapID = input_data %>% 
  full_join(trapID_dict)

# Obtain mosPerTrapNight for abundance analysis
abund = joinTrapID %>% 
  mutate(mosPerTrapNight = females_mixed/trap_nights) %>% 
  ungroup()

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
  abund, 
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


write_csv(abund, str_c(output_dir, "aedes_abund_mixedfemales.csv"))
write_csv(presence_absence_df, str_c(output_dir, "aedes_presenceAbsence_mixedfemales.csv"))
write_csv(presence_df, str_c(output_dir, "aedes_presence_mixedfemales.csv"))
