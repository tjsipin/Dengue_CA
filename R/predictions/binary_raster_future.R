#' Use MLR for stacking, not tidyverse

library(tidyverse)
library(stacks)
library(maxnet)
library(terra)
library(ggthemes)
tidymodels::tidymodels_prefer()

#Read in data

raster.master.path = "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/"
png.master.path = "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/"


training_L1 = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

training_L2 = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

testing = readRDS("Data/1_DataProcessing/modeling/testing_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

raster.master.path = "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/"
png.master.path = "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/"
dir.create(raster.master.path, recursive=T, showWarnings=T)
dir.create(png.master.path, recursive=T, showWarnings=T)


#Read in models
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgboost/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking/meta_learner_rf_mask_randomForest.rds")

#Read in thresholds
thresholds = read_csv("Data/1_DataProcessing/modeling/stacking/binary_thresholds.csv")

kableExtra::kbl(
    thresholds, 
    col.names = c("Model", "Minimum thresh.", "Thresh. that maximizes ROC", "P05 thresh.", "P10 thresh."), 
    caption = "Thresholds by Model",
    digits = 3
) %>%
    kableExtra::kable_classic(full_width = F, html_font = "Cambria")

months = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')



# Create binary rasters per month per model per threshold type ------------


makeBinaryRaster = function(mon, RCP=c("45", "60", "85"), decade=c("mid", "end"), type=c("min", "roc", "p05", "p10")){
    t0 = Sys.time()
    print(paste(mon, RCP, decade, type))
    
    m = tibble(
        months = months,
        mm = str_pad(1:12, 2, 'left', '0')
    ) %>% 
        filter(months == mon) %>% 
        pull(mm)
    
    input_rast_filename = list.files(
        path = paste0("Data/1_DataProcessing/modeling/stacking/rasters/means/future/", RCP),
        pattern = paste0(decade, "_\\d*", mon),
        full.names = T
    )
    
    if(length(input_rast_filename)==0) return(NULL)
    
    output_rast_filename = paste0(raster.master.path, RCP, "/", decade, "/", type, "/", m, "_", mon, "_", decade, "_", type, "_thresh.tif")
    dir.create(paste0(raster.master.path, RCP, "/", decade, "/", type), recursive = T, showWarnings = F)
    
    if(file.exists(output_rast_filename)){
        print(paste0("File exists: ", output_rast_filename))
        return(NULL)
    }
    
    if(length(input_rast_filename) > 1){
        print(input_rast_filename)
        return(NULL)
    }
    if(file.exists(output_rast_filename)){
        warning(paste0("File exists: ", output_rast_filename))
        return(NULL)
    }
    
    threshold_type = paste0(type, "_threshold")
    
    
    input_rast = rast(input_rast_filename)
    input_rast_xgb = input_rast$xgb_pred_mask_lc
    xgb_thresh = thresholds %>% 
        filter(model == "XGBoost") %>% 
        pull(threshold_type)
    output_rast_xgb = ifel(input_rast_xgb < xgb_thresh, 0, 1)
    
    input_rast_maxent = input_rast$maxent_pred_mask_lc
    maxent_thresh = thresholds %>% 
        filter(model == "MaxENT") %>% 
        pull(threshold_type)
    output_rast_maxent = ifel(input_rast_maxent < maxent_thresh, 0, 1)
    
    input_rast_meta = input_rast$meta_pred_mask_lc
    meta_thresh = thresholds %>% 
        filter(model == "Meta") %>% 
        pull(threshold_type)
    output_rast_meta = ifel(input_rast_meta < meta_thresh, 0, 1)
    
    output_rast = rast(list(output_rast_xgb, output_rast_maxent, output_rast_meta))
    
    writeRaster(output_rast, output_rast_filename, overwrite=T)
    print(Sys.time() - t0)
}

dir.create(paste0(raster.master.path, "min"), recursive=T)
dir.create(paste0(raster.master.path, "roc"), recursive=T)
dir.create(paste0(raster.master.path, "p05"), recursive=T)
dir.create(paste0(raster.master.path, "p10"), recursive=T)

#Min threshold
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='mid', type='min')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='end', type='min')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='mid', type='min')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='end', type='min')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='mid', type='min')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='end', type='min')
)

#Threshold that maximizes ROC
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='mid', type='roc')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='end', type='roc')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='mid', type='roc')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='end', type='roc')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='mid', type='roc')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='end', type='roc')
)

#p05 threshold
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='mid', type='p05')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='end', type='p05')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='mid', type='p05')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='end', type='p05')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='mid', type='p05')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='end', type='p05')
)

#p10 threshold
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='mid', type='p10')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='45', decade='end', type='p10')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='mid', type='p10')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='60', decade='end', type='p10')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='mid', type='p10')
)
map(
    months, 
    ~ makeBinaryRaster(mon=.x, RCP='85', decade='end', type='p10')
)

# Sum up binary values ----------------------------------------------------

sumBinary = function(type=c("min", "roc", "p05", "p10"), RCP=c("45", "60", "85"), decade=c("mid", "end")){
    
    t0 = Sys.time()
    output_filename = paste0("Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/", RCP, "_", decade, "_", type, "_thresh_sum.tif")
    
    # Get all binary XGBoost raster outputs for the RCP and decade, then sum 
    rasts_xgb = map(
        months,
        function(mon){
            m = tibble(
                months = months,
                mm = str_pad(1:12, 2, 'left', '0')
            ) %>% 
                filter(months == mon) %>% 
                pull(mm)
            
            input_rast_filename = paste0(
                "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", 
                RCP, "/",
                decade, "/", 
                type, "/", 
                m, "_", mon, "_", decade, "_", type, "_thresh.tif"
            )
            input_rast = rast(input_rast_filename)$xgb_pred_mask_lc
            names(input_rast) = paste0(mon, "_min")
            return(input_rast)
        }
    ) %>% 
        rast() %>% 
        sum()
    
    # Get all binary MaxENT raster outputs for the RCP and decade, then sum 
    rasts_maxent = map(
        months,
        function(mon){
            m = tibble(
                months = months,
                mm = str_pad(1:12, 2, 'left', '0')
            ) %>% 
                filter(months == mon) %>% 
                pull(mm)
            
            input_rast_filename = paste0(
                "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", 
                RCP, "/",
                decade, "/", 
                type, "/", 
                m, "_", mon, "_", decade, "_", type, "_thresh.tif"
            )
            input_rast = rast(input_rast_filename)$maxent_pred_mask_lc
            names(input_rast) = paste0(mon, "_min")
            return(input_rast)
        }
    ) %>% 
        rast() %>% 
        sum()
    
    # Get all binary Meta RF raster outputs for the RCP and decade, then sum 
    rasts_meta = map(
        months,
        function(mon){
            m = tibble(
                months = months,
                mm = str_pad(1:12, 2, 'left', '0')
            ) %>% 
                filter(months == mon) %>% 
                pull(mm)
            
            input_rast_filename = paste0(
                "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", 
                RCP, "/",
                decade, "/", 
                type, "/", 
                m, "_", mon, "_", decade, "_", type, "_thresh.tif"
            )
            input_rast = rast(input_rast_filename)$meta_pred_mask_lc
            names(input_rast) = paste0(mon, "_min")
            return(input_rast)
        }
    ) %>% 
        rast() %>% 
        sum()
    
    rasts_master = rast(list(rasts_xgb, rasts_maxent, rasts_meta))
    names(rasts_master) = c(paste0("xgb_", type), paste0("maxent_", type), paste0("meta_", type))
    
    writeRaster(rasts_master, output_filename, overwrite=T)
    print(Sys.time() - t0)
}

sumBinary(type='min', RCP='45', decade='mid')
sumBinary(type='roc', RCP='45', decade='mid')
sumBinary(type='p05', RCP='45', decade='mid')
sumBinary(type='p10', RCP='45', decade='mid')
sumBinary(type='min', RCP='45', decade='end')
sumBinary(type='roc', RCP='45', decade='end')
sumBinary(type='p05', RCP='45', decade='end')
sumBinary(type='p10', RCP='45', decade='end')

sumBinary(type='min', RCP='60', decade='mid')
sumBinary(type='roc', RCP='60', decade='mid')
sumBinary(type='p05', RCP='60', decade='mid')
sumBinary(type='p10', RCP='60', decade='mid')
sumBinary(type='min', RCP='60', decade='end')
sumBinary(type='roc', RCP='60', decade='end')
sumBinary(type='p05', RCP='60', decade='end')
sumBinary(type='p10', RCP='60', decade='end')

sumBinary(type='min', RCP='85', decade='mid')
sumBinary(type='roc', RCP='85', decade='mid')
sumBinary(type='p05', RCP='85', decade='mid')
sumBinary(type='p10', RCP='85', decade='mid')
sumBinary(type='min', RCP='85', decade='end')
sumBinary(type='roc', RCP='85', decade='end')
sumBinary(type='p05', RCP='85', decade='end')
sumBinary(type='p10', RCP='85', decade='end')


ggplotRasters = function(RCP=c('45', '60', '85'), decade=c('mid', 'end')){
    #Gather XGBoost rasters of differing thresholds
    min_xgb_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_min_thresh_sum.tif"
    ))$xgb_min
    roc_xgb_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_roc_thresh_sum.tif"
    ))$xgb_roc
    p05_xgb_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p05_thresh_sum.tif"
    ))$xgb_p05
    p10_xgb_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p10_thresh_sum.tif"
    ))$xgb_p10
    master_xgb_rast = rast(list(min_xgb_thresh, roc_xgb_thresh, p05_xgb_thresh, p10_xgb_thresh))


    #Gather MaxENT rasters of differing thresholds
    min_maxent_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_min_thresh_sum.tif"
    ))$maxent_min
    roc_maxent_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_roc_thresh_sum.tif"
    ))$maxent_roc
    p05_maxent_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p05_thresh_sum.tif"
    ))$maxent_p05
    p10_maxent_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p10_thresh_sum.tif"
    ))$maxent_p10
    master_maxent_rast = rast(list(min_maxent_thresh, roc_maxent_thresh, p05_maxent_thresh, p10_maxent_thresh))


    #Gather meta rasters of differing thresholds
    min_meta_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_min_thresh_sum.tif"
    ))$meta_min
    roc_meta_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_roc_thresh_sum.tif"
    ))$meta_roc
    p05_meta_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p05_thresh_sum.tif"
    ))$meta_p05
    p10_meta_thresh = rast(paste0(
        "Data/1_DataProcessing/modeling/stacking/rasters/means/future/binary/", RCP, "/",
        RCP, "_", decade, "_p10_thresh_sum.tif"
    ))$meta_p10
    master_meta_rast = rast(list(min_meta_thresh, roc_meta_thresh, p05_meta_thresh, p10_meta_thresh))

    if(decade=='mid'){
        caption_decade = "2040-2050"
    } else{
        caption_decade = "2090-2099"
    }
    if(RCP=='45'){
        caption_RCP = '4.5'
    } else if(RCP=='60'){
        caption_RCP = '6.0'
    } else if(RCP=='85'){
        caption_RCP = '8.5'
    }
    gg_caption = paste0("Years ", caption_decade, "with RCP ", caption_RCP)

    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(master_xgb_rast)) +
        facet_wrap(~lyr) +
        ggtitle("Months where Aedes Aegypti is present", "(thresholds using XGBoost v6)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        theme_void()
    ggsave(paste0("Data/1_DataProcessing/modeling/stacking/png/means/future/binary/", RCP, "_", decade, "_xgb_thresh.png"), height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)

    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(master_maxent_rast)) +
        facet_wrap(~lyr) +
        ggtitle("Months where Aedes Aegypti is present", "(thresholds using MaxENT v6)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        theme_void()
    ggsave(paste0("Data/1_DataProcessing/modeling/stacking/png/means/future/binary/", RCP, "_", decade, "_maxent_thresh.png"),  height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)

    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(master_meta_rast)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        facet_wrap(~lyr) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ggtitle("Months where Aedes Aegypti is present", "(thresholds using Meta v6)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:12)), discrete = T, option = "A") +
        theme_void()
    ggsave(paste0("Data/1_DataProcessing/modeling/stacking/png/means/future/binary/", RCP, "_", decade, "_meta_thresh.png"), height = 1080, width = 1080, units='px', bg = 'white', scale = 1.5)
}

ggplotRasters(RCP='45', decade='mid')
ggplotRasters(RCP='45', decade='end')
ggplotRasters(RCP='60', decade='mid')
ggplotRasters(RCP='60', decade='end')
ggplotRasters(RCP='85', decade='mid')
ggplotRasters(RCP='85', decade='end')


# Plot binary plots for each mean month -----------------------------------

ggplotRasterMeanMonth = function(mon, RCP=c('45', '60', '85'), decade=c('mid', 'end')){
    m = tibble(
        months = months,
        mm = str_pad(1:12, 2, 'left', '0')
    ) %>% 
        filter(months == mon) %>% 
        pull(mm)
    
    input_filename_min = paste0(raster.master.path, RCP, "/", decade, "/min/", m, "_", mon, "_", decade, "_min_thresh.tif")
    input_filename_roc = paste0(raster.master.path, RCP, "/", decade, "/roc/", m, "_", mon, "_", decade, "_roc_thresh.tif")
    input_filename_p05 = paste0(raster.master.path, RCP, "/", decade, "/p05/", m, "_", mon, "_", decade, "_p05_thresh.tif")
    input_filename_p10 = paste0(raster.master.path, RCP, "/", decade, "/p10/", m, "_", mon, "_", decade, "_p10_thresh.tif")
    
    input_rast_min = rast(input_filename_min)
    input_rast_roc = rast(input_filename_roc)
    input_rast_p05 = rast(input_filename_p05)
    input_rast_p10 = rast(input_filename_p10)
    
    dir.create(paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/min"
    ), recursive = T, showWarnings = F)
    dir.create(paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/roc"
    ), recursive = T, showWarnings = F)
    dir.create(paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/p05"
    ), recursive = T, showWarnings = F)
    dir.create(paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/p10"
    ), recursive = T, showWarnings = F)
    
    output_filename_min = paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/min/",
        m, "_", mon, "_", decade, "_RCP", RCP, "_min_thresh.png"
    )
    output_filename_roc = paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/roc/",
        m, "_", mon, "_", decade, "_RCP", RCP, "_roc_thresh.png"
    )
    output_filename_p05 = paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/p05/",
        m, "_", mon, "_", decade, "_RCP", RCP, "_p05_thresh.png"
    )
    output_filename_p10 = paste0(
        "Data/1_DataProcessing/modeling/stacking/png/means/future/binary/mean_month/", 
        RCP, "/", decade, "/p10/",
        m, "_", mon, "_", decade, "_RCP", RCP, "_p10_thresh.png"
    )
    
    if(decade=='mid'){
        caption_decade = '2040-2050'
    } else{
        caption_decade = '2090-2099'
    }
    gg_caption = paste0(str_to_sentence(mon), " ", caption_decade, " — RCP: ", RCP)
    
    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(input_rast_min)) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        facet_wrap(~lyr) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ggtitle("Binary classification of suitability for A. aegypti", "(Minimum threshold)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        theme_void()
    ggsave(output_filename_min, bg = 'white', scale = 5)
    
    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(input_rast_roc)) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        facet_wrap(~lyr) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ggtitle("Binary classification of suitability for A. aegypti", "(ROC threshold)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        theme_void()
    ggsave(output_filename_roc, bg = 'white', scale = 5)
    
    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(input_rast_p05)) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        facet_wrap(~lyr) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ggtitle("Binary classification of suitability for A. aegypti", "(P05 threshold)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        theme_void()
    ggsave(output_filename_p05, bg = 'white', scale = 5)
    
    ggplot() +
        tidyterra::geom_spatraster(data = as.factor(input_rast_p10)) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        facet_wrap(~lyr) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ggtitle("Binary classification of suitability for A. aegypti", "(P10 threshold)") +
        labs(caption = gg_caption) +
        viridis::scale_color_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        viridis::scale_fill_viridis(limits = factor(c(0:1)), discrete = T, option = "A") +
        theme_void()
    ggsave(output_filename_p10, bg = 'white', scale = 5)
}

map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='45', decade='mid')
)
map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='45', decade='end')
)
map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='60', decade='mid')
)
map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='60', decade='end')
)
map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='85', decade='mid')
)
map(
    .x = months,
    ~ ggplotRasterMeanMonth(.x, RCP='85', decade='end')
)