library(tidyverse)
library(stacks)
library(maxnet)
library(terra)
library(ggthemes)
library(dplyr)
library(purrr)
library(tidyterra)
tidymodels::tidymodels_prefer()

#Read in data
training_L1 = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

training_L2 = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

testing = readRDS("Data/1_DataProcessing/modeling/testing_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

raster.master.path = "Data/1_DataProcessing/modeling/stacking/rasters/monthly/binary/"
png.master.path = "Data/1_DataProcessing/modeling/stacking/png/monthly/binary/"
dir.create(raster.master.path, recursive=T, showWarnings=T)
dir.create(png.master.path, recursive=T, showWarnings=T)


#Read in models
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgboost/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking/meta_learner_rf_mask_randomForest.rds")

#Get thresholds for each model
training_L1_preds = training_L1 %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
        maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) 

training_L2_preds = training_L2 %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
        maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) %>% 
    mutate(
        meta_pred = predict(meta_model, ., 'prob')$.pred_1
    ) %>% 
    mutate(
        meta_pred_mask = ifelse(lc_type %in% c(2, 13, 14), meta_pred, 0)
    )

testing_preds = testing %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
        maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) %>% 
    mutate(
        meta_pred = predict(meta_model, ., 'prob')$.pred_1
    ) %>% 
    mutate(
        meta_pred_mask = ifelse(lc_type %in% c(2, 13, 14), meta_pred, 0)
    )

thresh_xgb_min = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(xgb_pred_mask > 0) %>% 
    pull(xgb_pred_mask) %>% 
    min() 

thresh_maxent_min = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(maxent_pred_mask > 0) %>% 
    pull(maxent_pred_mask) %>% 
    min() 

thresh_meta_min = training_L2_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(meta_pred_mask > 0) %>% 
    pull(meta_pred_mask) %>% 
    min() 

thresh_xgb_roc = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = training_L1_preds$pr_ab, 
            predictor= training_L1_preds$xgb_pred_mask, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

thresh_maxent_roc = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = training_L1_preds$pr_ab, 
            predictor= training_L1_preds$maxent_pred_mask, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

thresh_meta_roc = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = training_L2_preds$pr_ab, 
            predictor= training_L2_preds$meta_pred_mask, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

thresh_xgb_p05 = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(xgb_pred_mask > 0) %>% 
    pull(xgb_pred_mask) %>% 
    quantile(0.05)

thresh_maxent_p05 = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(maxent_pred_mask > 0) %>%
    pull(maxent_pred_mask) %>% 
    quantile(0.05)

thresh_meta_p05 = training_L2_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(meta_pred_mask > 0) %>% 
    pull(meta_pred_mask) %>% 
    quantile(0.05)

thresh_xgb_p10 = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(xgb_pred_mask > 0) %>% 
    pull(xgb_pred_mask) %>% 
    quantile(0.1)

thresh_maxent_p10 = training_L1_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(maxent_pred_mask > 0) %>% 
    pull(maxent_pred_mask) %>% 
    quantile(0.1)

thresh_meta_p10 = training_L2_preds %>% 
    filter(pr_ab == 1) %>% 
    filter(meta_pred_mask > 0) %>% 
    pull(meta_pred_mask) %>% 
    quantile(0.1)


thresholds = tibble(
    model = c("XGBoost", "MaxENT", "Meta"),
    min_threshold = c(thresh_xgb_min, thresh_maxent_min, thresh_meta_min),
    roc_threshold = c(thresh_xgb_roc, thresh_maxent_roc, thresh_meta_roc),
    p05_threshold = c(thresh_xgb_p05, thresh_maxent_p05, thresh_meta_p05),
    p10_threshold = c(thresh_xgb_p10, thresh_maxent_p10, thresh_meta_p10)
)

write_csv(thresholds, "Data/1_DataProcessing/modeling/stacking/binary_thresholds.csv")

kableExtra::kbl(
    thresholds, 
    col.names = c("Model", "Minimum thresh.", "Thresh. that maximizes ROC", "P05 thresh.", "P10 thresh."), 
    caption = "Thresholds by Model",
    digits = 3
) %>%
    kableExtra::kable_classic(full_width = F, html_font = "Cambria")

mons = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
mon_year = expand_grid(
    mon = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'),
    year = 2010:2023
)


makeBinaryRaster = function(this_m, this_y, type=c("min", "roc", "p05", "p10")){
    print(paste(this_m, this_y, type))
    
    # Collect month digits
    m = tibble(
        mons = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'),
        mm = str_pad(1:12, 2, 'left', '0')
    ) %>% 
        filter(mons==this_m) %>%
        pull(mm)
    
    # Call in input raster given month and year
    input_rast_filename = paste0("Data/1_DataProcessing/modeling/stacking/rasters/monthly/pred_", this_m, this_y, ".tif")
    if(!file.exists(input_rast_filename)) return(NULL)
    input_rast = rast(input_rast_filename)
    
    # Declare output raster filename
    output_rast_filename = paste0(raster.master.path, type, "/", m, "_", this_m, this_y, "_", type, "_thresh.tif")
    dir.create(paste0(raster.master.path, type), recursive = T, showWarnings = F)

    # If output raster exists, skip
    if(file.exists(output_rast_filename)){
        print(paste0("File exists: ", output_rast_filename))
        return(NULL)
    }

    # Declare threshold type as string
    threshold_type = paste0(type, "_threshold")
    
    # Call input raster
    input_rast = rast(input_rast_filename)
    # XGBoost prediction of input raster
    input_rast_xgb = input_rast$xgb_pred_mask
    xgb_thresh = thresholds %>% 
        filter(model == "XGBoost") %>% 
        pull(threshold_type)
    # If predicted probability is less than the threshold, set to 0. Otherwise, set to 1
    output_rast_xgb = ifel(input_rast_xgb < xgb_thresh, 0, 1)
    
    # MaxENT prediction of input raster
    input_rast_maxent = input_rast$maxent_pred_mask
    maxent_thresh = thresholds %>% 
        filter(model == "MaxENT") %>% 
        pull(threshold_type)
    # If predicted probability is less than the threshold, set to 0. Otherwise, set to 1
    output_rast_maxent = ifel(input_rast_maxent < maxent_thresh, 0, 1)
    
    # Meta prediction of input raster
    input_rast_meta = input_rast$meta_pred_mask
    meta_thresh = thresholds %>% 
        filter(model == "Meta") %>% 
        pull(threshold_type)
    # If predicted probability is less than the threshold, set to 0. Otherwise, set to 1
    output_rast_meta = ifel(input_rast_meta < meta_thresh, 0, 1)
    
    # Combine all binary outputs
    output_rast = rast(list(output_rast_xgb, output_rast_maxent, output_rast_meta))
    
    writeRaster(output_rast, output_rast_filename, overwrite=T)
}

for(year in 2010:2023){
    for(mon in mons){
        for(type in c("min", "roc", "p05", "p10")){
            makeBinaryRaster(this_m = mon, this_y = year, type = type)
        }
    }
}


ggplotMonthlyRasters = function(this_m, this_y){
    print(paste(this_m, this_y))
    
    m = tibble(
        mons = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'),
        mm = str_pad(1:12, 2, 'left', '0')
    ) %>% 
        filter(mons==this_m) %>%
        pull(mm)
    
    # input_min_filenames = list.files(
    #     path = paste0(raster.master.path, "min/"),
    #     pattern = this_m,
    #     full.names = T
    # )
    
    input_min_filename = paste0(raster.master.path, "min/", m, "_", this_m, this_y, "_min_thresh.tif")
    input_roc_filename = paste0(raster.master.path, "roc/", m, "_", this_m, this_y, "_roc_thresh.tif")
    input_p05_filename = paste0(raster.master.path, "p05/", m, "_", this_m, this_y, "_p05_thresh.tif")
    input_p10_filename = paste0(raster.master.path, "p10/", m, "_", this_m, this_y, "_p10_thresh.tif")
    
    input_min_rast = rast(input_min_filename)
    input_roc_rast = rast(input_roc_filename)
    input_p05_rast = rast(input_p05_filename)
    input_p10_rast = rast(input_p10_filename)
    
    xgb_rast = rast(list(
        input_min_rast$xgb_pred_mask, 
        input_roc_rast$xgb_pred_mask, 
        input_p05_rast$xgb_pred_mask, 
        input_p10_rast$xgb_pred_mask
    ))
    names(xgb_rast) = c("min", "roc", "p05", "p10")
    
    maxent_rast = rast(list(
        input_min_rast$maxent_pred_mask, 
        input_roc_rast$maxent_pred_mask, 
        input_p05_rast$maxent_pred_mask, 
        input_p10_rast$maxent_pred_mask
    ))
    names(maxent_rast) = c("min", "roc", "p05", "p10")
    
    meta_rast = rast(list(
        input_min_rast$meta_pred_mask, 
        input_roc_rast$meta_pred_mask, 
        input_p05_rast$meta_pred_mask, 
        input_p10_rast$meta_pred_mask
    ))
    
    names(meta_rast) = c("min", "roc", "p05", "p10")
    dir.create(paste0(png.master.path, "xgb"), showWarnings = F)
    dir.create(paste0(png.master.path, "maxent"), showWarnings = F)
    dir.create(paste0(png.master.path, "meta"), showWarnings = F)
    
    #XGBoost plot
    ggplot() +
        geom_spatraster(
            data=as.factor(xgb_rast)
        ) +
        facet_wrap(~lyr) + 
        labs(
            title = "XGBoost output on 4 thresholds",
            subtitle = paste(str_to_sentence(this_m), this_y)
        ) + 
        # viridis::scale_fill_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        scale_fill_discrete(limits = factor(c(0, 1))) +
        theme_void()
    ggsave(paste0(png.master.path, "xgb/", m, "_xgb_binary_", this_m, this_y, ".png"),  height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)
    
    #MaxENT plot
    ggplot() +
        geom_spatraster(
            data=as.factor(maxent_rast)
        ) +
        facet_wrap(~lyr) + 
        labs(
            title = "MaxENT output on 4 thresholds",
            subtitle = paste(str_to_sentence(this_m), this_y)
        ) + 
        scale_fill_discrete(limits = factor(c(0, 1))) +
        theme_void()
    ggsave(paste0(png.master.path, "maxent/", m, "_maxent_binary_", this_m, this_y, ".png"),  height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)
    
    #Meta plot
    ggplot() +
        geom_spatraster(
            data=as.factor(meta_rast)
        ) +
        facet_wrap(~lyr) + 
        labs(
            title = "Meta model output on 4 thresholds",
            subtitle = paste(str_to_sentence(this_m), this_y)
        ) + 
        scale_fill_discrete(limits = factor(c(0, 1))) +
        theme_void()
    ggsave(paste0(png.master.path, "meta/", m, "_meta_binary_", this_m, this_y, ".png"),  height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)
    
}


for(year in 2010:2023){
    for(mon in mons){
        ggplotMonthlyRasters(this_m = mon, this_y = year)
    }
}
dir.create(paste0(png.master.path, "sum"))

sumBinaryRaster = function(this_y, type=c('min', 'roc', 'p05', 'p10')){
    output_filename = paste0(png.master.path, "sum/", "meta_binary_sum_year=", this_y, "_type=", type, ".png")
    if(file.exists(output_filename)) return(NULL)
    input_rasts = list.files(
        path = paste0(raster.master.path, type),
        pattern = paste0(this_y),
        all.files = T,
        full.names = T
    )
    
    xgb_rast_sum = map(
        input_rasts,
        function(fi){
            r = rast(fi)
            return(r$xgb_pred_mask)
        }
    ) %>% 
        rast() %>% 
        sum()
    names(xgb_rast_sum) = "xgb"
    maxent_rast_sum = map(
        input_rasts,
        function(fi){
            r = rast(fi)
            return(r$maxent_pred_mask)
        }
    ) %>% 
        rast() %>% 
        sum()
    names(maxent_rast_sum) = "maxent"
    meta_rast_sum = map(
        input_rasts,
        function(fi){
            r = rast(fi)
            return(r$meta_pred_mask)
        }
    ) %>% 
        rast() %>% 
        sum()
    names(meta_rast_sum) = "meta"
    
    master_rast = rast(list(xgb_rast_sum, maxent_rast_sum, meta_rast_sum))
    
    ggplot() +
        geom_spatraster(
            data=as.factor(master_rast)
        ) +
        facet_wrap(~lyr) + 
        labs(
            title = paste0("Sum of binary output (", type, ")"),
            subtitle = paste0("For year ", this_y)
        ) + 
        viridis::scale_fill_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        theme_void()
    ggsave(output_filename,  height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)
}

for(year in 2010:2022){
    for(type in c("min", "roc", "p05", "p10")){
        print(paste(year, type))
        sumBinaryRaster(this_y = year, type = type)
    }
}

