library(raster)
library(terra)
library(tidyterra)
library(leaflet)
library(sf)
tidymodels_prefer()


xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")

# Get CA polygon for masking and cropping USGS data 
example_bcm_data = rast("Data/bcm/rasters/45/CA_BCM_MIROC_rcp45_2039-10-01.tif")
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs = "EPSG:3310") # for consistency with BCM data

months = c(
    'jan', 'feb', 'mar', 'apr', 'may', 'jun',
    'jul', 'aug', 'sep', 'oct', 'nov', 'dec'
)
outputHTML = function(month, decade=c('current', 'mid', 'end'), RCP=NULL){
    print(paste0(month, decade, RCP))
    t0 = Sys.time()
    mm = tibble(
        mon = months,
        mm = str_pad(1:12, width = 2, "0", side="left")
    ) %>% 
        filter(mon==month) %>% 
        pull(mm)
    
    output_filename = paste0(
        "Data/1_DataProcessing/modeling/stacking_G/sensitivity/leaflet/", 
        decade, "/", 
        RCP, "/", 
        "month=", month, 
        "+decade=", decade, 
        "+RCP=", RCP, "_leaflet.html"
    )
    
    # if(file.exists(output_filename)){
    #     return(NULL)
    # }
    
    if(decade=="current"){
        env_var_path = paste0(
            "Data/1_DataProcessing/df/env_lc"
        )
    } else{
        env_var_path = paste0(
            "Data/1_DataProcessing/df/", RCP
        )
    }
    
    if(decade=='current'){
        env_var_pattern = paste0(month)
    } else if(decade=='mid'){
        env_var_pattern = paste0(paste0(2040:2050, "-", mm, "-01"), collapse = "|")
    } else if(decade=='end'){
        env_var_pattern = paste0(paste0(2090:2099, "-", mm, "-01"), collapse = "|")
    }
    
    env_var_filenames = list.files(
        env_var_path,
        pattern = env_var_pattern,
        full.names = T
    )
    
    env_var_data_list = map(
        env_var_filenames,
        readRDS
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
                select(x, y, mean_temp, ppt, tws, lc_type)
            
            rast_i = rasterFromXYZ(
                data_i,
                crs = "EPSG:3310"
            ) %>% 
                rast() 
            
            rast_i
        }
    ) 
    Sys.time()-t0
    
    tmean_rast = map(
        env_var_rasts,
        function(r){
            r$mean_temp
        }
    ) %>% 
        rast() %>% 
        mean(na.rm = T)
    names(tmean_rast) = "mean_temp"
    ppt_rast = map(
        env_var_rasts,
        function(r){
            r$ppt
        }
    ) %>% 
        rast() %>% 
        mean(na.rm = T)
    names(ppt_rast) = "ppt"
    tws_rast = map(
        env_var_rasts,
        function(r){
            r$tws
        }
    ) %>% 
        rast() %>% 
        mean(na.rm = T)
    names(tws_rast) = "tws"
    lc_rast = map(
        env_var_rasts,
        function(r){
            r$lc_type
        }
    ) %>% 
        rast() %>% 
        modal(na.rm = T)
    names(lc_rast) = "lc_type"
    
    if(decade=="current"){
        pred_path = "Data/1_DataProcessing/modeling/stacking_G/rasters/means"
        pred_pattern = month
    } else{
        pred_path = paste0("Data/1_DataProcessing/modeling/stacking_G/rasters/means/future/", RCP)
        pred_pattern = paste0(decade, ".*", mon)
    }
    
    pred_file = list.files(
        pred_path,
        pred_pattern,
        full.names = T
    )
    
    pred_rast = rast(pred_file) %>% 
        crop(CA, mask=T)
    
    env_rast = rast(list(tmean_rast, ppt_rast, tws_rast, lc_rast)) %>% 
        crop(pred_rast, mask=T)
    
    comb_rast = rast(list(env_rast, pred_rast))
    
    
    
    CA_counties = tigris::counties(state="CA") %>% 
        sf::st_transform(., crs = "EPSG:3310")
    
    comb_rast_agg = aggregate(
        comb_rast,
        fact = 3
    ) %>% 
        terra::project("EPSG:3857")
    
    library(leaflet)
    
    
    temp_vals = values(comb_rast_agg$mean_temp)
    ppt_vals = values(comb_rast_agg$ppt)
    tws_vals = values(comb_rast_agg$tws)
    pred_vals = seq(0, 1, 0.000001)
    
    # temp_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), temp_vals,
    #                          na.color = "transparent")
    # ppt_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), ppt_vals,
    #                         na.color = "transparent")
    # tws_pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), tws_vals,
    #                         na.color = "transparent")
    # pred_pal <- colorNumeric(c("black", "#0C2C84", "#41B6C4", "#FFFFCC"), pred_vals,
    #                          na.color = "transparent")
    
    temp_pal <- colorNumeric("viridis", temp_vals,
                             na.color = "transparent")
    ppt_pal <- colorNumeric("viridis", ppt_vals,
                            na.color = "transparent")
    tws_pal <- colorNumeric("viridis", tws_vals,
                            na.color = "transparent")
    pred_pal <- colorNumeric("viridis", pred_vals,
                             na.color = "transparent")
    
    
    
    leaflet_output = leaflet() %>% 
        addTiles() %>%
        addRasterImage(comb_rast_agg$mean_temp, opacity = 0.8, colors = temp_pal, group = 'Mean temperature') %>% 
        addRasterImage(comb_rast_agg$ppt, opacity = 0.8, colors = ppt_pal, group = 'Precipitation') %>% 
        addRasterImage(comb_rast_agg$tws, opacity = 0.8, colors = tws_pal, group = 'TWS') %>% 
        addRasterImage(comb_rast_agg$xgb_pred_mask_lc, opacity = 0.8, colors = pred_pal, group = 'XGBoost') %>% 
        addRasterImage(comb_rast_agg$maxent_pred_mask_lc, opacity = 0.8, colors = pred_pal, group = 'MaxENT') %>% 
        addRasterImage(comb_rast_agg$meta_pred_mask_lc, opacity = 0.8, colors = pred_pal, group = 'Meta') %>% 
        addLegend("topleft", pal = temp_pal, title="Mean temperature",
                  values = temp_vals, group="Mean temperature") %>%
        addLegend("topleft", pal = ppt_pal, title="Precipitation",
                  values = ppt_vals, group="Precipitation") %>% 
        addLegend("bottomleft", pal = temp_pal, title="TWS",
                  values = tws_vals, group="TWS") %>%
        addLegend("bottomright", pal = pred_pal, title="Prediction",
                  values = pred_vals) %>% 
        addLayersControl(
            baseGroups = c("OSM (default)"),
            overlayGroups = c("Mean temperature", "Precipitation", "TWS", "XGBoost", "MaxENT", "Meta"),
            options = layersControlOptions(collapsed = FALSE)
        ) 
    
    htmlwidgets::saveWidget(leaflet_output, file=output_filename)
    print(Sys.time() - t0)
}


map(
    months,
    ~outputHTML(.x, decade="current")
)
map(
    months,
    ~outputHTML(.x, decade="mid", RCP="45")
)


