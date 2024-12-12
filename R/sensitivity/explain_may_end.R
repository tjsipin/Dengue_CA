library(raster)
library(terra)
library(tidyterra)
tidymodels_prefer()


xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")

env_var_filenames = list.files(
    "Data/1_DataProcessing/df/future/45/",
    pattern = paste0(paste0(2090:2099, "-05-01"), collapse = "|"),
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
Sys.time() - t0

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
    env_pred_rast_mean_pred_comb,
    fact = 3
) %>% 
    terra::project("EPSG:3857")

library(leaflet)


temp_vals = values(agg_env_pred_rast_mean_pred_comb_meta_pred$mean_temp)
ppt_vals = values(agg_env_pred_rast_mean_pred_comb_meta_pred$ppt)
tws_vals = values(agg_env_pred_rast_mean_pred_comb_meta_pred$tws)
pred_vals = seq(0, 1, 0.000001)

temp_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), temp_vals,
                    na.color = "transparent")
ppt_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), ppt_vals,
                         na.color = "transparent")
tws_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), tws_vals,
                         na.color = "transparent")
pred_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), pred_vals,
                         na.color = "transparent")



may_end_4.5_leaflet = leaflet() %>% 
    addTiles() %>%
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$mean_temp, opacity = 0.5, group = 'Mean temperature') %>% 
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$ppt, opacity = 0.5, group = 'Precipitation') %>% 
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$tws, opacity = 0.5, group = 'TWS') %>% 
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$xgb_pred_mask, opacity = 0.5, group = 'XGBoost') %>% 
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$maxent_pred_mask, opacity = 0.5, group = 'MaxENT') %>% 
    addRasterImage(agg_env_pred_rast_mean_pred_comb_meta_pred$.pred_1, opacity = 0.5, group = 'Meta') %>% 
    addLegend("bottomright", pal = temp_pal, title="Mean temperature",
              values = temp_vals, group="Mean temperature") %>%
    addLegend("bottomright", pal = ppt_pal, title="Precipitation",
              values = ppt_vals, group="Precipitation") %>% 
    addLegend("bottomright", pal = temp_pal, title="TWS",
              values = tws_vals, group="TWS") %>%
    addLegend("bottomright", pal = pred_pal, title="Prediction",
              values = pred_vals) %>% 
    addLayersControl(
        baseGroups = c("OSM (default)"),
        overlayGroups = c("Mean temperature", "Precipitation", "TWS", "XGBoost", "MaxENT", "Meta"),
        options = layersControlOptions(collapsed = FALSE)
    ) 

htmlwidgets::saveWidget(may_end_4.5_leaflet, file="Data/1_DataProcessing/modeling/stacking_G/sensitivity/leaflet/may_end_4.5_leaflet.html")

#' lisa couper meeting:
#' walkthrough of troubleshooting for SDM
#' aka biasing for non-analog conditions
#' 