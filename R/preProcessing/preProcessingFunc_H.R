library(sf) # to read in shapefiles
library(geosphere) # for distance functions
library(terra) # to read in rasters
# library(rgeos) # TODO
library(dplyr) # always
library(tidyr) # always
library(cluster) # TODO
library(lubridate) # for date functions
library(purrr) # for map function
library(sp) # for spatial points data frame
library(readr)
tidymodels::tidymodels_prefer()

# set seed
set.seed(123123)

splitData = function(data){
    data_v1 = data %>% 
        mutate(id = row_number())
    # create training and testing partitions
    data_split = rsample::initial_split(data_v1, strata = Year, prop = 0.8)
}

# create preprocess function to transform all data in training and testing by training data standardization
standardizeData = function(data, training){
    # gather median, mean, and standard deviation of each variable in training
    medians = training %>% 
        select(-pr_ab, -id, -Year, -Longitude, -Latitude, -Month, -lc_type) %>% 
        summarize_if(is.numeric, ~median(.x, na.rm = T))
    
    names_vec = base::intersect(names(data), names(medians))
    
    for(variable in names_vec){
        # print(variable)
        median_ = medians[, variable][[1]]
        data = data %>% 
            mutate_at(variable, ~replace_na(.x, median_)) 
    }
    return(data)
}


preprocessFunc = function(data_split){
    training = training(data_split) %>% 
        mutate(lc_type = as.factor(lc_type)) %>% 
        mutate(
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        ) %>% 
        # replace NA with median at month-year level
        group_by(Year, Month) %>% 
        mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE)))) %>%
        ungroup() %>% 
        # replace NA with median at year level
        group_by(Year) %>% 
        mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE)))) %>% 
        ungroup() %>% 
        # replace NA with median at full data set level
        mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))
    
    testing = testing(data_split) %>% 
        mutate(lc_type = as.factor(lc_type))
    
    # replace missing with training medians (crude, open to suggestions)
    
    testing = testing %>% 
        standardizeData(training) %>%
        mutate(
            pr_ab = as.factor(pr_ab),
            lc_type = as.factor(lc_type),
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        )
    
    return(list(
        training = training,
        testing = testing
    ))
}

