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

output_dir = "Data/1_DataProcessing/modeling/stacking/rasters/monthly/future"
dir.create(output_dir, recursive = T)

xgb_model = readRDS("Data/1_DataProcessing/modeling/xgboost/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking/meta_learner_rf_mask_randomForest.rds")


months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)

mid_mon_years = expand.grid(
    months = ordered(months, levels = months),
    years = 2040:2050
) 

end_mon_years = expand.grid(
    months = ordered(months),
    years = 2090:2099
)

# Get CA polygon for masking and cropping USGS data 
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs = st_crs("EPSG:3310")) # for consistency with BCM data

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
    .x = mid_mon_years$months,
    .y = mid_mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '45')
)

map2(
    .x = end_mon_years$months,
    .y = end_mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '45')
)

# RCP 6.0
map2(
    .x = mid_mon_years$months,
    .y = mid_mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '60')
)

map2(
    .x = end_mon_years$months,
    .y = end_mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '60')
)

# RCP 8.5
map2(
    .x = mid_mon_years$months,
    .y = mid_mon_years$years,
    ~ predRaster(mon = .x, year = .y, RCP = '85')
)

map2(
    .x = end_mon_years$months,
    .y = end_mon_years$years,
    function(.x, .y){
        predRaster(mon = .x, year = .y, RCP = '85')
    }
)

CA = tigris::counties(state = "CA") %>%
    sf::st_transform(crs="EPSG:3310") %>%
    vect()

monthMeans = function(mon, m, RCP, decade){
    if(decade == "mid"){
        years = 2040:2050
        plot_years = "2040-2050"
    } else if(decade == "end"){
        years = 2090:2099
        plot_years = "2090-2099"
    }
    print(paste(mon, m, RCP, decade))

    dir.create(paste0("Data/1_DataProcessing/modeling/stacking/rasters/means/future/", RCP), showWarnings = F, recursive = T)
    dir.create(paste0("Data/1_DataProcessing/modeling/stacking/png/means/future/", RCP), showWarnings = F, recursive = T)

    output_rast_filename = paste0("Data/1_DataProcessing/modeling/stacking/rasters/means/future/", RCP, "/", decade, "_", m, mon, "_means.tif")
    output_png_filename = paste0("Data/1_DataProcessing/modeling/stacking/png/means/future/", RCP, "/", decade, "_", m, mon, "_means.png")
    if(file.exists(output_png_filename)){
        return(NULL)
    }

    files = list.files(
        path=paste0("Data/1_DataProcessing/modeling/stacking/rasters/monthly/future/", RCP, "/"),
        pattern=paste(paste0(mon, years), collapse = "|"),
        full.names=T
    )

    ## Special note for RCP 8.5 and end of century:
    ## 8.5 Jan 2092 bad extent, skip
    if((RCP=="85") & (mon=='jan') & (decade=='end')){
        mod_years = c(2090, 2091, 2093:2099)
        files = list.files(
            path=paste0("Data/1_DataProcessing/modeling/stacking/rasters/monthly/future/", RCP, "/"),
            pattern=paste(paste0(mon, mod_years), collapse = "|"),
            full.names=T
        )
    }
    
    # Collect XGBoost predictions and take the mean
    xgbRasts = map(
        files,
        ~ rast(.x)[["xgb_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(xgbRasts) = "xgb_pred_mask_lc"

    # Collect MaxENT predictions and take the mean
    maxentRasts = map(
        files,
        ~ rast(.x)[["maxent_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(maxentRasts) = "maxent_pred_mask_lc"

    # Collect Meta model predictions and take the mean
    metaRasts = map(
        files,
        ~ rast(.x)[["meta_pred_mask"]]
    ) %>%
        rast() %>%
        terra::mean(na.rm=T)

    names(metaRasts) = "meta_pred_mask_lc"

    # Combine all individual model outputs
    masterRast = rast(list(xgbRasts, maxentRasts, metaRasts))
    writeRaster(masterRast, output_rast_filename, overwrite=T)

    # Plot the raster as ggplot and save
    ggplot() +
        tidyterra::geom_spatraster(
            data=masterRast
        ) +
        geom_spatvector(data=CA, color = "red3", fill=NA) +
        scale_fill_gradientn(colors = viridis::viridis(10), name = expression("Predicted probability\nof occurrence"), limits=c(0, 1)) +
        facet_wrap(~lyr) +
        labs(
            title = str_to_title(paste0(mon, ": Aedes Aegypti in California (", plot_years, " average)")),
            subtitle = paste0(substr(RCP, 1, 1), ".", substr(RCP, 2, 2))
        ) +
        ggthemes::theme_tufte(ticks = F)
    ggsave(output_png_filename, scale=1, bg = 'white')
}


map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '45', decade = 'mid')
)
map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '45', decade = 'end')
)

map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '60', decade = 'mid')
)
map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '60', decade = 'end')
)

map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '85', decade = 'mid')
)
map2(
    .x = months,
    .y = 1:12,
    ~ monthMeans(mon = .x, m = .y, RCP = '85', decade = 'end')
)
