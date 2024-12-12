library(raster)
library(terra)
library(tidyterra)
tidymodels_prefer()


xgb_moderasterxgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")

env_var_filenames = list.files(
    "Data/1_DataProcessing/df/future/45/",
    pattern = paste0(paste0(2090:2099, "-04-01"), collapse = "|"),
    full.names = T
)

env_var_data_list = map(
    env_var_filenames,
    readRDS
)

pred_list = map(
    env_var_data_list,
    function(env_data){
        env_data_clean = env_data %>% 
            mutate(
                mean_temp = (tmn + tmx)/2,
                tws = aet/pet
            ) %>% 
            select(x, y, mean_temp, ppt, tws, lc_type)
        
        pred_data = env_data_clean %>% 
            mutate(
                xgb_pred = predict(xgb_model, ., 'prob')$.pred_1,
                maxent_pred = predict(maxent_model, ., 'prob')$.pred_1
            )
        
        pred_data
    }
)

map(
    env_var_data_list,
    ~ function(env_data){
        env_data_clean = env_data %>% 
            mutate(
                mean_temp = (tmn + tmx)/2,
                tws = aet/pet
            ) %>% 
            select(mean_temp, ppt, tws)
        
        xgb_explainer = DALEX::explain(
            xgb_model,
            data = env_data_clean
        )
        xgb_profile = DALEX::model_profile(
            xgb_explainer,
            N = 10000
        )
        
        maxent_explainer = DALEX::explain(
            maxent_model,
            data = env_data_clean
        )
        maxent_profile = DALEX::model_profile(
            maxent_explainer,
            N = 10000
        )
        
        
    }
)


t0=Sys.time()
env_var_rasts = map(
    1:length(env_var_data_list),
    function(i){
        data_i = env_var_data_list[[i]] %>% 
            mutate(
                mean_temp = (tmn + tmx)/2,
                tws = aet/pet
            ) %>% 
            select(x, y, mean_temp, ppt, tws)
        
        rast_i = rasterFromXYZ(
            data_i,
            crs = "EPSG:3310"
        ) %>% 
            rast() 
        
        rast_i
    }
) 
Sys.time()-t0

env_var_plots = map(
    1:length(env_var_rasts),
    function(i){
        rast_i = env_var_rasts[[i]]
        
        g = ggplot() +
            geom_spatraster(data = rast_i) +
            facet_wrap(~lyr) +
            scale_fill_continuous(na.value = "transparent")
        
        return(g)
    }
)

env_pred_rasts = map(
    1:length(pred_list),
    function(i){
        date_i = env_var_filenames[[i]] %>% 
            str_split_i(pattern = "future_", 2) %>% 
            str_split_i(pattern = "-01_45", 1)
        
        data_i = pred_list[[i]]
        rast_i = rasterFromXYZ(
            data_i,
            crs = "EPSG:3310"
        ) %>% 
            rast() 
        
        # gg_i = ggplot() +
        #     geom_spatraster(data = rast_i) +
        #     facet_wrap(~lyr, scales="free") +
        #     scale_fill_continuous(na.value = "transparent")
        
        return(
            list(
                raster = rast_i,
                date = date_i
            )
        )
    }
)

env_pred_rast_mean = map(
    c("mean_temp", "ppt", "tws", "xgb_pred", "maxent_pred"),
    function(var){
        mean_var = map(
            1:length(env_pred_rasts),
            function(i){
                rast_i = env_pred_rasts[[i]]$raster[[var]]
                rast_i
            }
        ) %>% 
            rast() %>% 
            terra::mean(na.rm=T)
        names(mean_var) = var
        return(mean_var)
    }
) %>% 
    rast()

env_pred_rast_mean %>% plot()

names(env_pred_rast_mean) = c("mean_temp", "ppt", "tws", "xgb_pred_mask", "maxent_pred_mask")

env_pred_rast_mean_pred = predict(env_pred_rast_mean, meta_model, type = 'prob', na.rm = T)

env_pred_rast_mean_pred_comb = rast(list(
    env_pred_rast_mean, env_pred_rast_mean_pred$.pred_1
))

example_bcm_data = rast("Data/bcm/rasters/45/CA_BCM_MIROC_rcp45_2039-10-01.tif")
CA = tigris::counties(state="CA") %>% 
    sf::st_transform(., crs = sf::st_crs(example_bcm_data))

ggplot() + 
    geom_sf(data=CA, color = 'red') +
    geom_spatraster(data=env_pred_rast_mean_pred_comb$.pred_1, alpha = 0.75)

agg_env_pred_rast_mean_pred_comb_meta_pred = aggregate(
    env_pred_rast_mean_pred_comb$.pred_1,
    fact = 3
) %>% 
    terra::project("EPSG:3857")

library(leaflet)

pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(agg_env_pred_rast_mean_pred_comb_meta_pred),
                    na.color = "transparent")
leaflet() %>% 
    addTiles() %>%
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred, opacity = 0.8) %>% 
    addLegend(pal = pal, values = values(agg_env_pred_rast_mean_pred_comb_meta_pred), title = "Meta model predictions")

#' lisa couper meeting:
#' walkthrough of troubleshooting for SDM
#' aka biasing for non-analog conditions
#' 