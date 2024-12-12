


# Read in -----------------------------------------------------------------

## Data
testing = readRDS("Data/1_DataProcessing/modeling/testing_2016_2022_F.rds")
## Models
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")
## Thresholds
thresholds = read_csv("Data/1_DataProcessing/modeling/stacking_G/binary_thresholds.csv")

testing_pred = testing %>% 
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
    ) %>% 
    select(
        master_id, pr_ab, mean_temp, ppt, tws, lc_type, contains("pred_mask")
    )

testing_pred_thresholds = testing_pred %>% 
    mutate(
        xgb_thresh_min = ifelse(xgb_pred_mask < thresholds$min_threshold[1], 0, 1),
        xgb_thresh_roc = ifelse(xgb_pred_mask < thresholds$roc_threshold[1], 0, 1),
        xgb_thresh_p05 = ifelse(xgb_pred_mask < thresholds$p05_threshold[1], 0, 1),
        xgb_thresh_p10 = ifelse(xgb_pred_mask < thresholds$p10_threshold[1], 0, 1),
        
        maxent_thresh_min = ifelse(maxent_pred_mask < thresholds$min_threshold[2], 0, 1),
        maxent_thresh_roc = ifelse(maxent_pred_mask < thresholds$roc_threshold[2], 0, 1),
        maxent_thresh_p05 = ifelse(maxent_pred_mask < thresholds$p05_threshold[2], 0, 1),
        maxent_thresh_p10 = ifelse(maxent_pred_mask < thresholds$p10_threshold[2], 0, 1),
        
        meta_thresh_min = ifelse(meta_pred_mask < thresholds$min_threshold[3], 0, 1),
        meta_thresh_roc = ifelse(meta_pred_mask < thresholds$roc_threshold[3], 0, 1),
        meta_thresh_p05 = ifelse(meta_pred_mask < thresholds$p05_threshold[3], 0, 1),
        meta_thresh_p10 = ifelse(meta_pred_mask < thresholds$p10_threshold[3], 0, 1)
    ) %>% 
    pivot_longer(contains("thresh"), names_to = "model_thresh_type", values_to = "model_thresh_pred") %>% 
    mutate(model_thresh_pred = as.factor(model_thresh_pred)) %>% 
    mutate(error_type = case_when(
        pr_ab == model_thresh_pred ~ "Correct",
        (pr_ab == 1) & (model_thresh_pred == 0) ~ "False negative",
        (pr_ab == 0) & (model_thresh_pred == 1) ~ "False positive"
    ) %>% as.factor())
    
ggplot(data=testing_pred_thresholds) + 
    geom_density(
        aes(
            x = mean_temp,
            fill = model_thresh_pred
        ),
        alpha = 0.7
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

tmean_error_type = ggplot(data=testing_pred_thresholds) + 
    geom_density(
        aes(
            x = mean_temp,
            fill = error_type
        ),
        alpha = 0.7
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)


ppt_error_type = ggplot(data=testing_pred_thresholds) + 
    geom_density(
        aes(
            x = log(ppt + 1),
            fill = error_type
        ),
        alpha = 0.5
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

tws_error_type = ggplot(data=testing_pred_thresholds) + 
    geom_density(
        aes(
            x = log(tws + 1),
            fill = error_type
        ),
        alpha = 0.5
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

ggarrange(
    tmean_error_type,
    ppt_error_type,
    tws_error_type,
    common.legend = T
)

tmean_error_type_boxplot = ggplot(data=testing_pred_thresholds) + 
    geom_boxplot(
        aes(
            y = mean_temp,
            fill = error_type
        ),
        alpha = 0.7
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

ppt_error_type_boxplot = ggplot(data=testing_pred_thresholds) + 
    geom_boxplot(
        aes(
            y = log(ppt + 1),
            fill = error_type
        ),
        alpha = 0.7
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

tws_error_type_boxplot = ggplot(data=testing_pred_thresholds) + 
    geom_boxplot(
        aes(
            y = log(tws + 1),
            fill = error_type
        ),
        alpha = 0.7
    ) +
    facet_wrap(~model_thresh_type, nrow = 3)

ggarrange(
    tmean_error_type_boxplot,
    ppt_error_type_boxplot,
    tws_error_type_boxplot,
    common.legend = T
)

tmean_error_type_density_2d = ggplot(
    data = testing_pred_thresholds %>% 
        filter(model_thresh_type %in% c(
            "meta_thresh_min", "meta_thresh_roc",
            "meta_thresh_p05", "meta_thresh_p10"
        ))
) +
    geom_density_2d_filled(
        aes(
            x = mean_temp,
            y = meta_pred_mask,
            # fill = error_type
        )
    ) +
    facet_wrap(~ error_type + model_thresh_type)

ppt_error_type_density_2d = ggplot(data = testing_pred_thresholds) +
    geom_density_2d_filled(
        aes(
            x = log(ppt + 1),
            y = meta_pred_mask,
            # fill = error_type
        )
    ) +
    facet_wrap(model_thresh_type ~ error_type)

tws_error_type_density_2d = ggplot(data = testing_pred_thresholds) +
    geom_density_2d_filled(
        aes(
            x = tws,
            y = meta_pred_mask,
            # fill = error_type
        )
    ) +
    facet_wrap(model_thresh_type ~ error_type)
