
library(DALEX)
library(terra)
library(tidyterra)

#' Steps:
#' 1) Read in training data
#' 2) Get min/max temperatures
#' 3) Determine n observations for each month-year raster
# Read in training data
training_L1 = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
training_L2 = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
training = rbind(training_L1, training_L2)

# Read in models
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")

# Define min/max temperatures
max_temp = max(training$mean_temp)
min_temp = min(training$mean_temp)

### Define function to get number of pixels that are outside the extreme temperatures
nPixelsExtremeTemps = function(mon, year){
    print(paste0(mon, year))
    # Read in input data
    input_data_filename = paste0("Data/1_DataProcessing/df/env_lc/cleaned_env_", mon, year, ".rds")
    if(!file.exists(input_data_filename)) return(NULL)
    input_data = readRDS(input_data_filename)
    
    # Get filtered data
    filtered_data = input_data %>% 
        filter(
            (mean_temp > max_temp) | (mean_temp < min_temp)
        )
    # Define number of observations for month-year
    n = nrow(filtered_data)
    
    # Create output data frame
    output = tibble(
        mon = mon,
        year = year,
        n = n
    )
    return(output)
}

### Current conditions
months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
months_years = expand.grid(
    months = months,
    years = as.character(2010:2022)
)
extreme_current_n = map2(
    months_years$months,
    months_years$years,
    ~ nPixelsExtremeTemps(mon = .x, year = .y)
) %>% bind_rows()

extreme_current_mon_years = extreme_current_n %>% 
    filter(n > 0) %>% 
    select(mon, year)


getExtremePoints = function(mon, year){
    print(paste0(mon, year))
    # Read in input data
    input_data_filename = paste0("Data/1_DataProcessing/df/env_lc/cleaned_env_", mon, year, ".rds")
    if(!file.exists(input_data_filename)) return(NULL)
    input_data = readRDS(input_data_filename)
    
    
    # Get filtered data
    filtered_data = input_data %>% 
        filter(
            (mean_temp > max_temp) | (mean_temp < min_temp)
        ) %>% 
        mutate(tws = aet/pet) %>% 
        select(x, y, mean_temp, ppt, tws, lc_type)
    
    return(filtered_data)
    
}

extreme_points_current = map2(
    months_years$mon,
    months_years$year,
    ~ getExtremePoints(mon = .x, year = .y)
) %>% bind_rows()

## Get extreme point predictions
extreme_points_current_pred = extreme_points_current %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')$.pred_1,
        maxent_pred = predict(maxent_model, ., 'prob')$.pred_1
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

extreme_points_current_pred_long = extreme_points_current_pred %>% 
    select(mean_temp, ppt, tws, xgb_pred, maxent_pred) %>% 
    slice_sample(n = 10000) %>% 
    pivot_longer(-c(xgb_pred, maxent_pred), names_to = "env_name", values_to = "env_val") %>% 
    pivot_longer(-c(env_name, env_val), names_to = "pred_name", values_to = "pred_val")

extreme_points_current_pred_long %>% 
    ggplot() +
    geom_boxplot(
        aes(
            y = env_val
        )
    ) +
    facet_wrap(~ env_name, scales = "free") +
    labs(
        title = "Boxplots of environmental predictors for extreme temperature (current)"
    ) +
    theme_classic()
    

extreme_points_current_pred %>% 
    select(mean_temp, ppt, tws, xgb_pred, maxent_pred) %>% 
    mutate(extreme_temp = case_when(
        mean_temp < min_temp ~ "Low",
        mean_temp > max_temp ~ "High"
    )) %>% 
    slice_sample(n = 10000) %>% 
    pivot_longer(-c(xgb_pred, maxent_pred, extreme_temp), names_to = "env_name", values_to = "env_val") %>% 
    pivot_longer(-c(env_name, env_val, extreme_temp), names_to = "pred_name", values_to = "pred_val") %>% 
    ggplot() +
    geom_boxplot(
        aes(
            y = pred_val,
            color = extreme_temp
        )
    ) +
    labs(
        title = "Boxplots of predictions per extreme temperature type (current)"
    ) +
    theme_classic()

#### DALEX explain models
xgb_explainer = explain(
    model = xgb_model,
    data = extreme_points_current_pred %>% 
        select(mean_temp, ppt, tws)
)
xgb_profile = model_profile(
    explainer = xgb_explainer
)

xgb_profile %>% plot(title = "XGBoost")

maxent_explainer = explain(
    model = maxent_model,
    data = extreme_points_current_pred %>% 
        select(mean_temp, ppt, tws)
)
maxent_profile = model_profile(
    explainer = maxent_explainer
)

maxent_profile %>% plot(title = "MaxENT")

extreme_points_current_pred %>% 
    select(mean_temp, ppt, tws, xgb_pred, maxent_pred) %>% 
    slice_sample(n = 10000) %>% 
    pivot_longer(-c(xgb_pred, maxent_pred), names_to = "env_name", values_to = "env_val") %>% 
    pivot_longer(-c(env_name, env_val), names_to = "pred_name", values_to = "pred_val") %>% 
    ggplot() +
    geom_point(
        aes(
            x = env_val,
            y = pred_val,
            color = pred_name
        )
    ) +
    facet_wrap(~ env_name, scales = "free")


currentExtremeRast = function(mon){
    print(paste0(mon))
    output_rast_filename = paste0("Data/1_DataProcessing/modeling/stacking_G/sensitivity/extreme/current/", mon, "_extreme.tif")
    dir.create("Data/1_DataProcessing/modeling/stacking_G/sensitivity/extreme/current")
    # if(file.exists(output_rast_filename)) return(NULL)
    # Read in input data
    input_rast_filename = paste0("Data/1_DataProcessing/rasters/bcm/means/current/", mon, "_env_mean.tif")
    input_rast = rast(input_rast_filename)
    
    msk = ifel((input_rast$mean_temp < min_temp) | (input_rast$mean_temp > max_temp), 1, NA) 
    input_rast_mask = input_rast %>% 
        terra::mask(msk) 
    input_rast_mask %>% 
        plot()
    
    writeRaster(output_rast, output_rast_filename, overwrite=T)
}
map(months, currentExtremeRast)

feb_current_rast = rast(paste0("Data/1_DataProcessing/modeling/stacking_G/sensitivity/extreme/current/", "feb", "_extreme.tif"))

### Define function to get number of pixels that are outside the extreme temperatures
nPixelsExtremeTemps_v2 = function(mon, decade, RCP){
    print(paste0(mon, decade, RCP))
    # Read in input raster
    input_rast_filename = paste0("Data/1_DataProcessing/rasters/bcm/means/", decade, "/", RCP, "/", mon, "_env_mean.tif")
    if(!file.exists(input_rast_filename)) return(NULL)
    input_rast = rast(input_rast_filename)
    
    input_df = input_rast %>% 
        as.data.frame() 
    
    temp_filter_df = input_df %>% 
        filter((mean_temp < min_temp) | (mean_temp > max_temp))
    
    # Define number of observations for month-year
    n = nrow(temp_filter_df)
    
    # Create output data frame
    output = tibble(
        mon = mon,
        decade = decade,
        RCP = RCP,
        n = n
    )
    return(output)
}



current_combs_df = tibble(
    mon = ordered(months),
    decade = "current",
    RCP = ""
)

mid_combs_df = expand_grid(
    mon = ordered(months),
    decade = "mid",
    RCP = c("45", "60", "85")
)

end_combs_df = expand_grid(
    mon = ordered(months),
    decade = "end",
    RCP = c("45", "60", "85")
)

combs_df = current_combs_df %>% 
    rbind(mid_combs_df) %>% 
    rbind(end_combs_df) 

n_pixels_extreme_temps_v2 = map(
    1:nrow(combs_df),
    ~ nPixelsExtremeTemps_v2(combs_df$mon[.x], combs_df$decade[.x], combs_df$RCP[.x])
) %>% 
    bind_rows()


### Summer for end decade of RCP 8.5 looks like a good place to visualize
#### Get mean temperature for these months
summer_end_8.5_mean_temp = rast(
    list(
        rast(paste0("Data/1_DataProcessing/rasters/bcm/means/", "end", "/", "85", "/", "jun", "_env_mean.tif"))$mean_temp,
        rast(paste0("Data/1_DataProcessing/rasters/bcm/means/", "end", "/", "85", "/", "jul", "_env_mean.tif"))$mean_temp,
        rast(paste0("Data/1_DataProcessing/rasters/bcm/means/", "end", "/", "85", "/", "aug", "_env_mean.tif"))$mean_temp
    )
)
names(summer_end_8.5_mean_temp) = c("Mean temperature for June", "Mean temperature for July", "Mean temperature for August")

summer_end_8.5_mean_temp_mask = ifel((summer_end_8.5_mean_temp < min_temp) | (summer_end_8.5_mean_temp > max_temp), 1, NA)
summer_end_8.5_mean_temp_masked = summer_end_8.5_mean_temp %>% 
    mask(summer_end_8.5_mean_temp_mask) 

ggplot() + 
    geom_spatvector(data=CA) +
    geom_spatraster(data = summer_end_8.5_mean_temp_masked, na.rm = T) +
    facet_wrap(~lyr) +
    scale_fill_continuous(na.value = "transparent") +
    theme_classic() +
    labs(
        title = "Mean temperatures for summer in the 2090s (RCP 8.5)"
    )
    
#### Get predictions for these months
summer_end_8.5_preds = map(
    c("jun", "jul", "aug"),
    function(m){
        data_rast = rast(paste0("Data/1_DataProcessing/rasters/bcm/means/", "end", "/", "85", "/", m, "_env_mean.tif"))
        
        temp_mask = ifel((data_rast$mean_temp < min_temp) | (data_rast$mean_temp > max_temp), 1, NA)
        
        data_rast_masked = data_rast %>% 
            mask(temp_mask)
        
        indivPredFunc = function(...){
            predict(..., type='prob')$.pred_1
            
        }
        t0 = Sys.time()
        xgb_pred_rast = terra::predict(object=data_rast_masked, model=xgb_model, fun=indivPredFunc, na.rm=T)
        names(xgb_pred_rast) = "xgb_pred_mask"
        xgb_pred_rast = rast(list(data_rast_masked, xgb_pred_rast))
        # xgb_pred_rast$xgb_pred_mask = ifel(xgb_pred_rast$lc_type %in% c(2, 13, 14), xgb_pred_rast$xgb_pred, 0)
        Sys.time() - t0
        maxent_pred_rast = terra::predict(object=data_rast_masked, model=maxent_model, fun=indivPredFunc, na.rm=T)
        names(maxent_pred_rast) = "maxent_pred_mask"
        maxent_pred_rast = rast(list(data_rast_masked, maxent_pred_rast))
        # maxent_pred_rast$maxent_pred_mask = ifel(maxent_pred_rast$lc_type %in% c(2, 13, 14), maxent_pred_rast$maxent_pred, 0)
        Sys.time() - t0
        
        
        indiv_pred_rast = rast(
            list(xgb_pred_rast$xgb_pred_mask, maxent_pred_rast$maxent_pred_mask) # not really mask but need to name it for meta model
        ) %>% 
            mask(CA)
        metaPredFunction = function(...){
            p = predict(..., 'prob')$.pred_1
            p
        }
        t0 = Sys.time()
        meta_pred_rast = terra::predict(indiv_pred_rast, meta_model, fun=metaPredFunction, na.rm=T) %>% 
            mask(CA)
        names(meta_pred_rast) = "meta_pred"
        meta_pred_rast_combine = rast(list(
            indiv_pred_rast,
            meta_pred_rast
        ))
        
        meta_pred_rast_combine
    }
)
    
summer_end_8.5_preds = rast(list(
    summer_end_8.5_preds[[1]]$xgb_pred_mask,
    summer_end_8.5_preds[[2]]$xgb_pred_mask,
    summer_end_8.5_preds[[3]]$xgb_pred_mask,
    summer_end_8.5_preds[[1]]$maxent_pred_mask,
    summer_end_8.5_preds[[2]]$maxent_pred_mask,
    summer_end_8.5_preds[[3]]$maxent_pred_mask,
    summer_end_8.5_preds[[1]]$meta_pred,
    summer_end_8.5_preds[[2]]$meta_pred,
    summer_end_8.5_preds[[3]]$meta_pred
))

names(summer_end_8.5_preds) = c(
    "Mean prediction for June (XGBoost)", "Mean prediction for July (XGBoost)", "Mean prediction for August (XGBoost)",
    "Mean prediction for June (MaxENT)", "Mean prediction for July (MaxENT)", "Mean prediction for August (MaxENT)",
    "Mean prediction for June (Meta)", "Mean prediction for July (Meta)", "Mean prediction for August (Meta)"
)

ggplot() + 
    geom_spatvector(data=CA) +
    geom_spatraster(data = summer_end_8.5_preds, na.rm = T) +
    facet_wrap(~lyr, nrow = 3) +
    scale_fill_continuous(na.value = "transparent") +
    theme_classic() +
    labs(
        title = "Mean predictions for summer in the 2090s (RCP 8.5)"
    )

ggplot() +
    geom_spatvector(data=CA) +
    geom_spatraster(data = summer_end_8.5_preds, na.rm = T) +
    facet_wrap(~lyr, nrow = 3) +
    xlim(c(150000, 539724.2)) + 
    ylim(c(-450000, 0)) +
    scale_fill_continuous(na.value = "transparent") +
    theme_classic() +
    labs(
        title = "Mean predictions for summer in the 2090s (RCP 8.5)"
    )