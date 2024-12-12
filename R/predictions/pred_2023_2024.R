library(tidyverse)
library(tidymodels)
library(ranger)
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

output_dir = "Data/1_DataProcessing/modeling/stacking_H/rasters/monthly/calibration/"
dir.create(output_dir, recursive = T)

xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_H/meta_learner_rf_mask_randomForest.rds")

months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)

mon_years = expand.grid(
    months = ordered(months, levels = months),
    years = 2023:2025
) 
# Get CA polygon for masking and cropping USGS data 
example_bcm_data = rast("Data/bcm/rasters/45/CA_BCM_MIROC_rcp45_2039-10-01.tif")
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs = st_crs(example_bcm_data)) # for consistency with BCM data

predRaster = function(mon, year, RCP){
    out.dir = paste0(output_dir, "/", RCP, "/")
    dir.create(out.dir, recursive=T, showWarnings = F)
    output_file = paste0(out.dir, "pred_", mon, year, ".tif")
    
    month = mon_month %>% 
        filter(mon_abb == mon) %>% 
        pull(Month)
    
    ymd = ymd(paste0(year, "-", month, "-01"))
    
    data_filename = paste0("Data/1_DataProcessing/df/future/", RCP, "/future_", ymd, "_", RCP, ".rds")
    
    
    
    
    print(paste0(RCP, mon, year))
    if(!file.exists(data_filename)) return(NULL)
    if(file.exists(output_file)) return(NULL)
    data = readRDS(data_filename) %>% 
        mutate(
            mean_temp = (tmx + tmn)/2,
            tws = aet/pet
        ) %>% 
        select(x, y, mean_temp, ppt, tws, lc_type)
    if(nrow(data)==0) return(NULL)
    
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
    
    rm(data)
    
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


# RCP 4.5
map2(
    .x = mon_years$months,
    .y = mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '45')
)
# RCP 6.0
map2(
    .x = mon_years$months,
    .y = mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '60')
)
# RCP 8.5
map2(
    .x = mon_years$months,
    .y = mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '85')
)