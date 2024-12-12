#' Author: TJ Sipin
#' Date: April 22, 2024
#' Input:
    #' BCM data
    #' USGS data
    #' Presence/pseudo-absence data
#' Output:
    #' Combined data set with environmental, land cover, and vector data
    

setwd("/home/tjsipin/network-storage/Dengue_CA/")
# Packages ----------------------------------------------------------------

library(readr)
library(dplyr)


# Data --------------------------------------------------------------------

bcm_env = read_csv("Data/1_DataProcessing/df/BCM_p_psa_extract_2016_2023.csv")
usgs_lc = read_csv("Data/1_DataProcessing/df/USGS_p_psa_extract_2016_2023.csv")

data = full_join(bcm_env, usgs_lc) 

write_csv(data, "Data/1_DataProcessing/df/aedesAegypti_env_lc.csv")
