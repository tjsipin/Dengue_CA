#' Author: TJ Sipin
#' Date: Oct 23, 2024



library(tidyverse)
library(tidymodels)
library(rsample)
library(caret)
library(ranger)
library(mlrMBO)
library(lhs)

tidymodels_prefer()
# Parameters --------------------------------------------------------------

formula = formula("pr_ab ~ ppt + mean_temp + tws")
master.path = "Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/"
dir.create(master.path, recursive = T)

# read in functions
source("R/1_DataProcessing/H/modeling/xgb_modeling_functions_H.R")


# -------------------------------------------------------------------------

model_per_month = map(
    .x = 1:12,
    .f = ~ modelPerMonth(
        m = .x, 
        this_formula = formula, 
        this_training = training_clusters_full,
        master.path = master.path,
        training_type = 'clusters'
    )
)


metrics_per_month = map(
    .x = 1:12,
    .f = function(m){
        metrics_filename = paste0("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/clusters/", m, "/xgb_all_folds_metrics_hurdle.rds")
        metrics = readRDS(metrics_filename) %>% 
            mutate(m=m)
        
        metrics
    }
) %>% 
    bind_rows() 

metrics_per_month_wider = metrics_per_month %>% 
    group_by(param_set, test_fold) %>% 
    summarize(
        mean_auc = mean(auc)
    ) %>% 
    pivot_wider(names_from = "test_fold", values_from = "mean_auc") %>% 
    ungroup() %>% 
    mutate(mean_auc = rowMeans(select(., `1`:`4`)))

best_parameters_setnumber = metrics_per_month_wider %>% 
    arrange(desc(mean_auc)) %>% 
    filter(row_number()==1) %>% 
    pull(param_set)

best_parameters = metrics_per_month %>% 
    filter(param_set == best_parameters_setnumber) %>% 
    select(tree_depth, learn_rate, mtry, min_n, loss_reduction) %>%
    distinct()


best_model = boost_tree(
    mode = "classification",
    engine = "xgboost",
    tree_depth = best_parameters$tree_depth,
    learn_rate = best_parameters$learn_rate,
    mtry = best_parameters$mtry,
    min_n = best_parameters$min_n,
    loss_reduction = best_parameters$loss_reduction
) %>%
    fit(formula, data = training_clusters_full)
saveRDS(best_model, "Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_best_model_clusters.rds")
best_model = readRDS("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_best_model_clusters.rds")

testing_pred = testing %>%
    mutate(pred_prob = predict(best_model, ., 'prob')[, 2][[1]])

saveRDS(testing_pred, "Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_predictions_clusters.rds")
testing_pred = readRDS("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_predictions_clusters.rds")

measureAUC(
    probabilities = testing_pred$pred_prob,
    truth = as.integer(as.character(testing_pred$pr_ab)),
    positive = "1"
) # 0.8610141

xgb_explain = DALEX::explain(
    best_model,
    data = testing_pred %>%
        select(pr_ab, tws, ppt, mean_temp, tmean_ppt),
    y = as.numeric(testing_pred$pr_ab),
    label = "XGBoost"
)

DALEX::model_profile(xgb_explain) %>%
    plot()
DALEX::variable_importance(xgb_explain) %>% 
    plot()
