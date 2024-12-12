# library(readr)
# # Data --------------------------------------------------------------------
# 
# data = read_csv("Data/1_DataProcessing/df/aedesAegypti_env_lc_E.csv") %>% 
#     mutate(pr_ab = factor(pr_ab))
# 
# 
# source("R/1_DataProcessing/F/preProcessing/preProcessingFunc_F.R")
# source("R/1_DataProcessing/F/preProcessing/spatialFoldsFunc_F.R")
# set.seed(123123)
# t0 = Sys.time()
# split_data = splitData(data)
# training = training(split_data)
# testing = testing(split_data)
# 
# preprocessed_data = preprocessFunc(split_data)
# Sys.time() - t0
# t1 = Sys.time()
# training = preprocessed_data$training 
# Sys.time() - t1
# testing = preprocessed_data$testing 
# 
# training_L = training
# training_L$L_level = sample(1:3, size = nrow(training_L), replace = T)
# training_L1_XGB = training_L %>% 
#     filter(L_level == 1)
# training_L1_MaxENT = training_L %>% 
#     filter(L_level == 2)
# training_L2_Meta = training_L %>% 
#     filter(L_level == 3)
# 
# 
# saveRDS(training_L1_XGB, "Data/1_DataProcessing/modeling/training_L1XGB_2016_2022_F.rds")
# saveRDS(training_L1_MaxENT, "Data/1_DataProcessing/modeling/training_L1Maxent_2016_2022_F.rds")
# saveRDS(training_L2_Meta, "Data/1_DataProcessing/modeling/training_L1Meta_2016_2022_F.rds")
# saveRDS(testing, "Data/1_DataProcessing/modeling/testing_2016_2022_F.rds")
# 
# training_L1_XGB_cv = training_L1_XGB %>% 
#     getSpatialFolds_v2()
# training_L1_MaxENT_cv = training_L1_MaxENT %>% 
#     getSpatialFolds_v2()
# training_L2_Meta_cv = training_L2_Meta %>% 
#     getSpatialFolds_v2()
# 
# Sys.time() - t1
# 
# saveRDS(training_L1_XGB_cv, "Data/1_DataProcessing/modeling/training_L1XGB_cv_2016_2022_F.rds")
# saveRDS(training_L1_MaxENT_cv, "Data/1_DataProcessing/modeling/training_L1_MaxENT_cv_2016_2022_F.rds")
# saveRDS(training_L2_Meta_cv, "Data/1_DataProcessing/modeling/training_L2_Meta_cv_2016_2022_F.rds")
# Sys.time() - t1
# 
# 
# ##################
library(dplyr)
library(readr)
# Data --------------------------------------------------------------------

data = read_csv("Data/1_DataProcessing/df/aedesAegypti_env_lc_H.csv") %>% 
    mutate(pr_ab = factor(pr_ab))


source("R/1_DataProcessing/H/preProcessing/preProcessingFunc_H.R")
source("R/1_DataProcessing/H/preProcessing/spatialFoldsFunc_H.R")
set.seed(123123)
t0 = Sys.time()
split_data = splitData(data)
training = training(split_data)
testing = testing(split_data)

preprocessed_data = preprocessFunc(split_data)
Sys.time() - t0
t1 = Sys.time()
training = preprocessed_data$training 
Sys.time() - t1
testing = preprocessed_data$testing 

training_L = training
training_L$L_level = sample(1:2, size = nrow(training_L), replace = T)
training_L1 = training_L %>% 
    filter(L_level == 1)
training_L2 = training_L %>% 
    filter(L_level == 2)


saveRDS(training_L1, "Data/1_DataProcessing/modeling/training_L1_2016_2023_H.rds")
saveRDS(training_L2, "Data/1_DataProcessing/modeling/training_L2_2016_2023_H.rds")


saveRDS(testing, "Data/1_DataProcessing/modeling/testing_2016_2023_H.rds")

training_L1_cv = training_L1 %>% 
    getSpatialFolds_v2()
training_L2_cv = training_L2 %>% 
    getSpatialFolds_v2()

saveRDS(training_L1_cv, "Data/1_DataProcessing/modeling/training_L1_cv_2016_2023_H.rds")
saveRDS(training_L2_cv, "Data/1_DataProcessing/modeling/training_L2_cv_2016_2023_H.rds")
