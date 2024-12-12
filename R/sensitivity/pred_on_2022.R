# XGBoost model variables: tws, ppt, mean temp
# MaxENT model variables: tws, ppt, mean temp

rm(list = ls())
library(tidyverse)
library(stacks)
library(maxnet)
library(tidyverse)
library(tidymodels)
library(ranger)
library(terra)
library(tidyverse)
library(stacks)
library(maxnet)
library(tidysdm)
library(ggpubr)

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
library(tidyterra)
tidymodels::tidymodels_prefer()
theme_set(theme_pubr())

source("R/1_DataProcessing/G/modeling/randomForest_stack_modeling_functions_G.R")

months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

#Read in data
training_L1 = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) %>% 
    filter(Year != 2022)
training_L2 = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) %>% 
    filter(Year != 2022)
## training_L2 = predict(base models)
testing = readRDS("Data/1_DataProcessing/modeling/testing_2016_2022_F.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab))

training_clusters = readRDS("Data/1_DataProcessing/modeling/training_L2_cv_2016_2022_F_v2.rds")$spatial_clusters
training_blocks = readRDS("Data/1_DataProcessing/modeling/training_L2_cv_2016_2022_F_v2.rds")$spatial_blocks

training_clusters_full = map(
    1:length(training_clusters$id),
    function(s0){
        s = training_clusters$splits[[s0]]
        as = rsample::assessment(s) %>%
            mutate(spat_fold = s0)
        as
    }
) %>% bind_rows() %>%
    sf::st_drop_geometry() %>% 
    mutate(L_level = 2) %>% 
    filter(Year != 2022)

training_blocks_full = map(
    1:length(training_blocks$id),
    function(s0){
        s = training_blocks$splits[[s0]]
        as = rsample::assessment(s) %>%
            mutate(spat_fold = s0)
        as
    }
) %>% bind_rows() %>%
    sf::st_drop_geometry() %>% 
    mutate(L_level = 2) %>% 
    filter(Year != 2022)

output.path = "Data/1_DataProcessing/modeling/stacking_G/sensitivity/2022/"
dir.create(output.path)

#Monthly model output training data
##XGBoost
xgb_model_original = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")

### XGBoost best params
    xgb_metrics_per_month = map(
        .x = 1:12,
        .f = function(m){
            metrics_filename = paste0("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/clusters/", m, "/xgb_all_folds_metrics_hurdle.rds")
            metrics = readRDS(metrics_filename) %>% 
                mutate(m=m)
            
            metrics
        }
    ) %>% 
        bind_rows() 
    
    xgb_metrics_per_month_wider = xgb_metrics_per_month %>% 
        group_by(param_set, test_fold) %>% 
        summarize(
            mean_auc = mean(auc)
        ) %>% 
        pivot_wider(names_from = "test_fold", values_from = "mean_auc") %>% 
        ungroup() %>% 
        mutate(mean_auc = rowMeans(select(., `1`:`4`)))
    
    xgb_best_parameters_setnumber = xgb_metrics_per_month_wider %>% 
        arrange(desc(mean_auc)) %>% 
        filter(row_number()==1) %>% 
        pull(param_set)
    
    xgb_best_parameters = xgb_metrics_per_month %>% 
        filter(param_set == xgb_best_parameters_setnumber) %>% 
        select(tree_depth, learn_rate, mtry, min_n, loss_reduction) %>%
        distinct()
#####

xgb_model_params = xgb_model_original$fit$params
xgb_model = boost_tree(
    mode = "classification",
    tree_depth = xgb_best_parameters$tree_depth,
    learn_rate = xgb_best_parameters$learn_rate,
    mtry = xgb_best_parameters$mtry,
    min_n = xgb_best_parameters$min_n,
    loss_reduction =  xgb_best_parameters$loss_reduction
) %>% 
    fit(
        formula = formula(
            "pr_ab ~ mean_temp + ppt + tws"
        ),
        data = training_L1
    )
# xgb_predictions = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_predictions_clusters.rds") %>% 
#     rename(xgb_pred = pred_prob) %>% 
#     mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))

xgb_predictions = training_L2 %>% 
    mutate(xgb_pred = predict(xgb_model, ., 'prob')$.pred_1) %>% 
    mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))
##MaxEnt
# maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
maxent_model = maxent(
    mode = "classification",
    engine = "maxnet", 
    feature_classes = "lqph",
    regularization_multiplier = 0.4256475 # taken from maxent CV
) %>% 
    fit(pr_ab ~ tws + ppt + mean_temp, data = training_L1)
# maxent_predictions = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_predictions_clusters.rds") %>% 
#     rename(maxent_pred = pred_prob) %>% 
#     mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))
maxent_predictions = training_L2 %>% 
    mutate(maxent_pred = predict(maxent_model, ., 'prob')$.pred_1) %>% 
    mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))

##Join together both individual models
indiv_predictions = full_join(xgb_predictions, maxent_predictions)

meta_learner_rf_mask = rand_forest(
    mode = "classification",
    engine = "randomForest",
    mtry = 1,
    trees = 2084,
    min_n = 16
) %>% fit(
    formula = formula(
        "pr_ab ~ xgb_pred_mask + maxent_pred_mask"
    ),
    data = indiv_predictions
)

testing_pred = testing %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., "prob")$.pred_1,
        maxent_pred = predict(maxent_model, ., "prob")$.pred_1
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) %>% 
    mutate(
        meta_pred = predict(meta_learner_rf_mask, ., "prob")$.pred_1
    ) %>% 
    mutate(
        meta_pred_mask = ifelse(lc_type %in% c(2, 13, 14), meta_pred, 0)
    )

saveRDS(xgb_model, "Data/1_DataProcessing/modeling/stacking_G/2022/xgb_model.rds")
saveRDS(maxent_model, "Data/1_DataProcessing/modeling/stacking_G/2022/maxent_model.rds")
saveRDS(meta_learner_rf_mask, "Data/1_DataProcessing/modeling/stacking_G/2022/meta_model.rds")
saveRDS(testing_pred, "Data/1_DataProcessing/modeling/stacking_G/2022/testing_pred.rds")
testing_pred = readRDS("Data/1_DataProcessing/modeling/stacking_G/2022/testing_pred.rds")

mon_month = tibble(
    mon_abb = months,
    Month = as.factor(1:12)
)

# Get CA polygon for masking and cropping USGS data 
example_bcm_data = rast("Data/1_DataProcessing/rasters/bcm/aet2014jan.tif")
CA = tigris::counties(state = "CA") %>% 
    sf::st_transform(crs="EPSG:3310") %>% 
    vect()

predRaster = function(mon){
    
    output_rast_filename = paste0(output.path, "pred_", mon, ".tif")
    output_png_filename = paste0(output.path, "pred_", mon, ".png")
    
    # if(file.exists(output_rast_filename)){
    #     return(NULL)
    # }
    month = mon_month %>% 
        filter(mon_abb == mon) %>% 
        pull(Month)
    
    pred_rast = rast(output_rast_filename)
    
    # data_filename = paste0("Data/1_DataProcessing/df/env_lc/cleaned_env_", mon, "2022.rds")
    # 
    # 
    # t0 = Sys.time()
    # print(paste0(mon, "2022"))
    # if(!file.exists(data_filename)) return(NULL)
    # if(file.exists(output_rast_filename)) return(NULL)
    # data = readRDS(data_filename)
    # 
    # 
    # indiv_pred_df = data %>% 
    #     mutate(
    #         lc_type = as.factor(as.character(lc_type))
    #     ) %>% 
    #     mutate(
    #         xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
    #         maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    #     ) %>% 
    #     mutate(
    #         xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
    #         maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    #     ) %>% 
    #     mutate(
    #         Month = month,
    #         Mon = mon %>% as.factor()
    #     )
    # Sys.time() - t0
    # 
    # rm(data)
    # 
    # t0 = Sys.time()
    # meta_pred_df = indiv_pred_df %>%
    #     select(x, y, Month, xgb_pred_mask, maxent_pred_mask, lc_type) %>% 
        # mutate(
        #     meta_pred = predict(meta_learner_rf_mask, ., 'prob')$.pred_1
        # ) %>%
        # mutate(
        #     meta_pred_mask = ifelse(lc_type %in% c(2, 13, 14), meta_pred, 0)
        # ) %>%
        # select(x, y, xgb_pred_mask, maxent_pred_mask, meta_pred_mask)
    # Sys.time() - t0
    # 
    # pred_rast = raster::rasterFromXYZ(meta_pred_df, crs = "EPSG:3310") %>% 
    #     terra::rast()
    # rm(indiv_pred_df)
    # rm(meta_pred_df)
    
    testing_m = testing_pred %>% 
        filter(Year == 2022) %>% 
        filter(
            Month == month
        ) %>% 
        vect(geom = c("x", "y"), crs="EPSG:3310")
            
    
    raster_plot = ggplot() +
        geom_spatvector(
            data = CA
        ) +
        geom_spatraster(
            data = pred_rast,
            alpha = 0.8
        ) +
        scale_fill_gradientn(colors = alpha(viridis::viridis(10), 0.8), name = expression("Predicted probability\nof occurrence"), limits=c(0, 1)) +
        facet_wrap(~lyr) +
        labs(
            title = "2022 predictions trained on pre-2022 data", 
            subtitle = str_to_sentence(mon)
        ) +
        theme(axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank(),
              legend.text = element_text(size=7))
    
    points_plot = ggplot() +
        geom_spatvector(
            data = CA
        ) +
        geom_spatvector(
            data = testing_m,
            aes(color = pr_ab),
            size = 0.7
        ) + 
        labs(
            caption = "Recorded presences and pseudo-absences"
        ) +
        theme(axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank())
    
    testing_m_preds = testing_m %>% 
        select(pr_ab, contains("mask")) %>% 
        rename(
            XGBoost = xgb_pred_mask,
            MaxENT = maxent_pred_mask, 
            Meta = meta_pred_mask
        ) %>% 
        pivot_longer(-pr_ab, names_to = "Model", values_to = "Prediction") %>% 
        as.data.frame() %>% 
        tibble()
    density_plot = ggplot() +
        geom_density(
            data = testing_m_preds,
            aes(
                fill = pr_ab,
                x = Prediction
            ),
            alpha = 0.7
        ) +
        facet_wrap(~ Model)
    
    # Arrange the plots on the same page
    ggarrange(
        raster_plot, 
        ggarrange(
            points_plot, density_plot,
            ncol = 2,
            widths = c(0.7, 1)
        ), 
        ncol = 1, nrow = 2,
        heights = c(1, 1),
        widths = c(0.8, 1)
    )
    
    ggsave(output_png_filename, scale = 4)
    # terra::writeRaster(pred_rast, output_rast_filename, overwrite=T)
    
    # rm(pred_rast)
}

map(
    months,
    predRaster
)


# Code appendix -----------------------------------------------------------


# # change the default setting of scales::colour_ramp:
# assignInNamespace("colour_ramp", function(colors, na.color = NA, alpha = TRUE){
#     # if (length(colors) == 0) {
#     #     stop("Must provide at least one color to create a color ramp")
#     # }
#     # colorMatrix <- grDevices::col2rgb(colors, alpha = alpha)
#     # structure(function(x) {
#     #     scales:::colour_ramp(colorMatrix, x, alpha)
#     # }, safe_palette_func = TRUE)
#     
#     # colour_ramp <- function(colors, na.color = NA, alpha = TRUE) {
#     if (length(colors) == 0) {
#         cli::cli_abort("Must provide at least one colour to create a colour ramp")
#     }
#     
#     if (length(colors) == 1) {
#         return(structure(
#             function(x) {
#                 ifelse(is.na(x), na.color, colors)
#             },
#             safe_palette_func = TRUE
#         ))
#     }
#     
#     # farver is not currently case insensitive, but col2rgb() is
#     colors <- tolower(colors)
#     lab_in <- farver::decode_colour(
#         colour = colors,
#         alpha = TRUE,
#         to = "lab",
#         na_value = "transparent"
#     )
#     
#     x_in <- seq(0, 1, length.out = length(colors))
#     l_interp <- stats::approxfun(x_in, lab_in[, 1])
#     u_interp <- stats::approxfun(x_in, lab_in[, 2])
#     v_interp <- stats::approxfun(x_in, lab_in[, 3])
#     
#     if (!alpha || all(lab_in[, 4] == 1)) {
#         alpha_interp <- function(x) NULL
#     } else {
#         alpha_interp <- stats::approxfun(x_in, lab_in[, 4])
#     }
#     
#     structure(
#         function(x) {
#             lab_out <- cbind(l_interp(x), u_interp(x), v_interp(x))
#             out <- farver::encode_colour(lab_out, alpha = alpha_interp(x), from = "lab")
#             out[is.na(out)] <- na.color
#             out
#         },
#         safe_palette_func = TRUE
#     )
#     # }
#     
# }, "scales")
