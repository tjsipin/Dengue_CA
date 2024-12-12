#' Use MLR for stacking, not tidyverse

library(tidyverse)
library(stacks)
library(maxnet)
library(tidysdm)
tidymodels::tidymodels_prefer()
setwd("/home/tjsipin/network-storage/Dengue_CA")
#Read in data
source("R/1_DataProcessing/H/modeling/maxent_modeling_functions_H.R")
master.path = "Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/"
dir.create(master.path, recursive=T)

des = data.frame(
    regularization_multiplier = c(1.1106301, 0.3527918, 0.8924429, 1.4073909, 0.4823738)
)

formula = formula("pr_ab ~ ppt + mean_temp + tws")

model_per_month = map(
    .x = 1,
    .f = ~ modelPerMonth(
        m = .x, 
        this_formula = formula, 
        this_training = training_clusters_full,
        master.path = master.path,
        training_type = 'clusters'
    )
)

all_metrics_per_m = map(
    .x = 1:12,
    .f = function(m){
        metrics_filename = paste0("Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/clusters/", m, "/maxent_all_folds_metrics_hurdle.rds")
        metrics = readRDS(metrics_filename)
        metrics %>%
            mutate(m = m)
    }
) %>%
    bind_rows()

all_metrics_per_m %>% 
    mutate(
        param_set = as.factor(param_set),
        test_fold = as.factor(test_fold)
    ) %>% 
    ggplot() + 
    geom_boxplot(
        aes(
            y = auc,
            color = test_fold,
            group = test_fold
        )
    ) + 
    facet_wrap(~ param_set, nrow = 1) +
    labs(
        title = "Performance of MaxENT models (n=5) faceted by spatial folds"
    ) +
    theme_classic()

all_metrics_per_m %>% 
    mutate(
        regularization_multiplier = as.factor(round(regularization_multiplier, 4)),
        test_fold = as.factor(test_fold)
    ) %>% 
    ggplot() + 
    geom_line(
        aes(
            x = m,
            y = auc,
            color = test_fold
        )
    ) +
    facet_wrap(~regularization_multiplier, nrow=1) +
    labs(
        title = "Performance of MaxENT models (n=5) faceted by spatial folds"
    ) +
    theme_classic()

all_metrics_per_m %>% 
    group_by(param_set, regularization_multiplier) %>% 
    summarize(
        mean_auc = mean(auc),
        median_auc = median(auc),
        sd_auc = sd(auc)
    ) %>% 
    ungroup() 
## param set 5 looks best
best_regularization_multiplier = all_metrics_per_m %>% 
    filter(param_set==5) %>% 
    pull(regularization_multiplier) %>% 
    unique()

best_model = maxent(
    mode = "classification",
    engine = "maxnet", 
    feature_classes = "lqph",
    regularization_multiplier = best_regularization_multiplier
) %>% 
    fit(pr_ab ~ tws + ppt + mean_temp, data = training_clusters_full)

saveRDS(best_model, "Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_best_model_clusters.rds")
best_model = readRDS("Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_best_model_clusters.rds")

zoe_model = maxent(
    mode = "classification",
    engine = "maxnet", 
    feature_classes = "lqph",
    regularization_multiplier = 1
) %>% 
    fit(pr_ab ~ tws + ppt + mean_temp, data = training_clusters_full)
saveRDS(zoe_model, "Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_zoe_model.rds")

testing_pred = testing %>%
    mutate(pred_prob = predict(best_model, ., 'prob')[, 2][[1]])
saveRDS(testing_pred, "Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_predictions_clusters.rds")
testing_pred_zoe = testing %>%
    mutate(pred_prob = predict(zoe_model, ., 'prob')[, 2][[1]])
saveRDS(testing_pred_zoe, "Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_predictions_clusters_zoe.rds")

measureAUC(
    probabilities = testing_pred_zoe$pred_prob,
    truth = as.integer(as.character(testing_pred_zoe$pr_ab)),
    positive = "1"
) # 0.8490551

maxent_explain = DALEX::explain(
    best_model,
    data = testing_pred %>%
        select(pr_ab, tws, ppt, mean_temp, lc_type),
    y = as.numeric(testing_pred$pr_ab),
    label = "v8"
)

DALEX::model_profile(maxent_explain) %>%
    plot()
