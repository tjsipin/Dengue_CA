library(tidyverse)
library(tidymodels)
library(randomForest)
library(terra)
library(tidyterra)

library(tidyverse)
library(tidymodels)
library(rsample)
library(caret)
library(ranger)
library(mlrMBO)
library(lhs)
library(leaflet)
library(xgboost)
library(maxnet)
library(sf)

conflicted::conflicts_prefer(terra::resample)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# set working directory
# create directory to output rasters 
setwd("/home/tjsipin/network-storage/Dengue_CA/")
dir.create("Data/1_DataProcessing/modeling/stacking/rasters/monthly", recursive = T)
months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgboost/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking/meta_learner_rf_mask_randomForest.rds")

mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)

mon_years = expand.grid(
    months = ordered(months),
    years = 2000:2023
) 

# Get CA polygon for masking and cropping USGS data 
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs = st_crs("EPSG:3310")) # for consistency with BCM data

predRaster = function(mon, year){
    out.dir = paste0("Data/1_DataProcessing/modeling/stacking/rasters/monthly/")
    output_file = paste0(out.dir, "pred_", mon, year, ".tif")
    
    if(file.exists(output_file)){
        return(NULL)
    }
    month = mon_month %>% 
        filter(mon_abb == mon) %>% 
        pull(Month)
    
    data_filename = paste0("Data/1_DataProcessing/df/env_lc/cleaned_env_", mon, year, ".rds")
    
    
    print(paste0(mon, year))
    if(!file.exists(data_filename)) return(NULL)
    if(file.exists(output_file)) return(NULL)
    data = readRDS(data_filename)
    
    data_rast = raster::rasterFromXYZ(
        data,
        crs = crs(example_bcm_data)
    ) %>%
        mask(CA) %>% 
        rast()
    
    
    indivPredFunc = function(...){
        predict(..., type='prob')$.pred_1
        
    }
    t0 = Sys.time()
    xgb_pred_rast = terra::predict(object=data_rast, model=xgb_model, fun=indivPredFunc, na.rm=T)
    names(xgb_pred_rast) = "xgb_pred"
    xgb_pred_rast = rast(list(data_rast, xgb_pred_rast))
    xgb_pred_rast$xgb_pred_mask = ifel(xgb_pred_rast$lc_type %in% c(2, 13, 14), xgb_pred_rast$xgb_pred, 0)
    Sys.time() - t0
    maxent_pred_rast = terra::predict(object=data_rast, model=maxent_model, fun=indivPredFunc, na.rm=T)
    names(maxent_pred_rast) = "maxent_pred"
    maxent_pred_rast = rast(list(data_rast, maxent_pred_rast))
    maxent_pred_rast$maxent_pred_mask = ifel(maxent_pred_rast$lc_type %in% c(2, 13, 14), maxent_pred_rast$maxent_pred, 0)
    Sys.time() - t0
    
    indiv_pred_rast = rast(
        list(xgb_pred_rast$xgb_pred_mask, maxent_pred_rast$maxent_pred_mask)
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
        data_rast$lc_type,
        indiv_pred_rast,
        meta_pred_rast
    ))
    meta_pred_rast_combine$meta_pred_mask = ifel(meta_pred_rast_combine$lc_type %in% c(2, 13, 14), meta_pred_rast_combine$meta_pred, 0) %>% 
        mask(CA)
    meta_pred_rast_final = rast(list(
        meta_pred_rast_combine$meta_pred_mask,
        meta_pred_rast_combine$xgb_pred_mask,
        meta_pred_rast_combine$maxent_pred_mask
    ))
    
    terra::writeRaster(meta_pred_rast_final, output_file, overwrite=T)
    Sys.time()-t0
}

months_years = expand.grid(
    months = months,
    years = as.character(2010:2023)
)

map2(
    .x = months_years$months,
    .y = months_years$years,
    ~ predRaster(mon = .x, year = .y)
)


CA_counties = tigris::counties(state = "CA") %>% 
    sf::st_transform(crs="EPSG:3310") %>%
    vect()
# Get means of each month of year -----------------------------------------

monthMeans = function(mon, m){
    output_rast_filename = paste0("Data/1_DataProcessing/modeling/stacking/rasters/means/", mon, "_means.tif")
    output_png_filename = paste0("Data/1_DataProcessing/modeling/stacking/png/means/", m, '_', mon, "_means.png")
    
    files = list.files(
        path=paste0("Data/1_DataProcessing/modeling/stacking/rasters/monthly/"),
        pattern=paste0(mon),
        full.names=T
    )

    xgbRasts = map(
        files,
        ~ rast(.x)[["xgb_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(xgbRasts) = "xgb_pred_mask_lc"

    maxentRasts = map(
        files,
        ~ rast(.x)[["maxent_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(maxentRasts) = "maxent_pred_mask_lc"

    metaRasts = map(
        files,
        ~ rast(.x)[["meta_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(metaRasts) = "meta_pred_mask_lc"

    masterRast = rast(list(xgbRasts, maxentRasts, metaRasts))
    writeRaster(masterRast, output_rast_filename, overwrite=T)
    masterRast = rast(output_rast_filename)
    
    ggplot() +
        tidyterra::geom_spatraster(
            data=masterRast
        ) +
        geom_spatvector(data=CA_counties, color = "red3", fill=NA) +
        scale_fill_gradientn(colors = viridis::viridis(10), name = expression("Predicted probability\nof occurrence"), limits=c(0, 1)) +
        facet_wrap(~lyr) +
        ggtitle(str_to_title(paste0(mon, ": Aedes Aegypti in California (current average)"))) +
        ggthemes::theme_tufte(ticks = F)

    ggsave(output_png_filename, scale=5, bg = 'white')
}

dir.create("Data/1_DataProcessing/modeling/stacking/rasters/means/", recursive=T)
dir.create("Data/1_DataProcessing/modeling/stacking/png/means/", recursive=T)

map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y)
)
