# XGBoost model variables: tws, ppt, mean temp
# MaxENT model variables: tws, ppt, mean temp

rm(list = ls())
library(tidyverse)
library(stacks)
library(maxnet)
tidymodels::tidymodels_prefer()

source("R/1_DataProcessing/H/modeling/randomForest_stack_modeling_functions_H.R")

#Read in data
training_L2 = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2023_H.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
## training_L2 = predict(base models)

training_clusters = readRDS("Data/1_DataProcessing/modeling/training_L2_cv_2016_2023_H.rds")$spatial_clusters
training_blocks = readRDS("Data/1_DataProcessing/modeling/training_L2_cv_2016_2023_H.rds")$spatial_blocks

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
    mutate(L_level = 2)

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
    mutate(L_level = 2)

testing = readRDS("Data/1_DataProcessing/modeling/testing_2016_2023_H.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

output.path = "Data/1_DataProcessing/modeling/stacking_H/"
dir.create(output.path)

#Monthly model output training data
##XGBoost
xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_best_model_clusters.rds")
# xgb_predictions = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_predictions_clusters.rds") %>% 
#     rename(xgb_pred = pred_prob) %>% 
#     mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))

xgb_predictions = training_clusters_full %>% 
    mutate(xgb_pred = predict(xgb_model, ., 'prob')$.pred_1) %>% 
    mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))
##MaxEnt
maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_zoe_model.rds")
# maxent_predictions = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_predictions_clusters.rds") %>% 
#     rename(maxent_pred = pred_prob) %>% 
#     mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))
maxent_predictions = training_clusters_full %>% 
    mutate(maxent_pred = predict(maxent_model, ., 'prob')$.pred_1) %>% 
    mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))

##Join together both individual models
indiv_predictions = full_join(xgb_predictions, maxent_predictions)

##Start CV process
map(
    1:12,
    ~ modelPerMonth(
        m = .x, this_formula = "pr_ab ~ xgb_pred_mask + maxent_pred_mask",
        master.path = output.path, this_training = indiv_predictions, training_type = "clusters"
    )
)
##See best hyperparameter set
metrics_per_month = map(
    1:12,
    function(m){
        metrics_filename = paste0("Data/1_DataProcessing/modeling/stacking_H/clusters/", m, "/rf_all_folds_metrics_hurdle.rds")
        metrics = readRDS(metrics_filename)
        metrics %>%
            mutate(
                Month = m,
                test_fold = as.factor(test_fold)
            ) %>%
            select(-param_set)
    }
) %>%
    bind_rows()

paramset_key = metrics_per_month %>%
    select(mtry, trees, min_n) %>%
    distinct() %>%
    mutate(param_set = row_number() %>% as.factor())

metrics_per_month %>%
    full_join(paramset_key) %>%
    ggplot() +
    geom_boxplot(
        aes(
            y = auc,
            group = param_set,
            color = param_set
        )
    ) +
    facet_wrap(~test_fold) +
    labs(
        title = "AUC by sptial fold and parameter set"
    )

metrics_per_month_wider = metrics_per_month %>%
    full_join(paramset_key) %>%
    group_by(param_set, test_fold) %>%
    summarize(
        mean_auc = mean(auc, na.rm = T),
        median_auc = median(auc, na.rm = T),
        sd_auc = sd(auc, na.rm = T)
    ) %>%
    pivot_wider(names_from = "test_fold", values_from = c("mean_auc", "median_auc", "sd_auc")) %>%
    ungroup() %>%
    mutate(
        mean_auc = rowMeans(select(., mean_auc_1:mean_auc_4)),
        median_auc = rowMeans(select(., median_auc_1:median_auc_4)),
        sd_auc = rowMeans(select(., sd_auc_1:sd_auc_4))
    )
metrics_per_month_wider %>%
    select(param_set, mean_auc, median_auc, sd_auc, everything()) %>%
    arrange(desc(mean_auc))

## Choosing param_set==7
best_params = metrics_per_month %>%
    full_join(paramset_key) %>%
    filter(param_set==7) %>%
    select(mtry, trees, min_n) %>%
    distinct()

meta_learner_rf = rand_forest(
    mode = "classification",
    engine = "ranger",
    mtry = best_params$mtry,
    trees = best_params$trees,
    min_n = best_params$min_n
) %>% fit(
    formula = formula(
        "pr_ab ~ xgb_pred + maxent_pred"
    ),
    data = indiv_predictions
)

meta_learner_rf_mask = rand_forest(
    mode = "classification",
    engine = "randomForest",
    mtry = best_params$mtry,
    trees = best_params$trees,
    min_n = best_params$min_n
) %>% fit(
    formula = formula(
        "pr_ab ~ xgb_pred_mask + maxent_pred_mask"
    ),
    data = indiv_predictions
)

saveRDS(meta_learner_rf_mask, paste0(output.path, "meta_learner_rf_mask_randomForest.rds"))
saveRDS(meta_learner_rf, paste0(output.path, "meta_learner_rf.rds"))
meta_learner_rf = readRDS(paste0(output.path, "meta_learner_rf.rds"))
meta_learner_rf_mask = readRDS(paste0(output.path, "meta_learner_rf_mask_randomForest.rds"))
#
meta_testing_preds = testing %>%
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')$.pred_1,
        maxent_pred = predict(maxent_model, ., 'prob')$.pred_1
    ) %>%
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) %>%
    mutate(
        rf_pred_meta_mask = predict(meta_learner_rf_mask, ., 'prob')$.pred_1,
    ) %>%
    mutate(
        rf_pred_meta_mask =ifelse(lc_type %in% c(2, 13, 14), rf_pred_meta_mask, 0),
    )
#
#
meta_testing_preds %>%
    mutate(pr_ab = as.integer(as.character(pr_ab))) %>% 
    group_by(Month) %>%
    summarize(
        # `Probability of Presence` = mean(pr_ab),
        `Meta Model` = measureAUC(
            rf_pred_meta_mask,
            pr_ab,
            positive = 1
        ),
        XGBoost = measureAUC(
            xgb_pred_mask,
            pr_ab,
            positive = 1
        ),
        MaxENT = measureAUC(
            maxent_pred_mask,
            pr_ab,
            positive = 1
        )
    ) %>%
    ungroup() %>%
    # pivot_longer(-c(Month, `Probability of Presence`), names_to = "Model", values_to = "AUC") %>%
    pivot_longer(-c(Month), names_to = "Model", values_to = "AUC") %>%
    group_by(Model) %>% 
    mutate(`Mean AUC` = mean(AUC)) %>% 
    ungroup() %>% 
    ggplot() +
    geom_line(
        aes(
            x = Month,
            y = AUC,
            color = Model
        )
    ) +
    geom_abline(aes(intercept = `Mean AUC`, slope = 0, color = Model), linetype="dashed") +
    # geom_col(
    #     aes(
    #         x = Month,
    #         y = `Probability of Presence`/3
    #     )
    # ) +
    scale_y_continuous(
        # Add a second axis and specify its features
        sec.axis = sec_axis(~., name="Probability of Presence")
    ) +
    theme_classic() +
    labs(
        title = "AUCs of Models per Month (with average across all months)"
    )

####Default parameter set:
# # A tibble: 12 × 2
# # Month   auc
# # <dbl> <dbl>
# # 1     1 0.992
# # 2     2 0.917
# # 3     3 0.952
# # 4     4 0.891
# # 5     5 0.921
# # 6     6 0.949
# # 7     7 0.934
# # 8     8 0.925
# # 9     9 0.913
# # 10    10 0.937
# # 11    11 0.903
# # 12    12 0.962

####Tuned parameter set:
# A tibble: 12 × 4
# Month `Meta Model` XGBoost MaxENT
# <dbl>        <dbl>   <dbl>  <dbl>
#     1        0.858   0.863  0.827
#     2        0.879   0.862  0.856
#     3        0.972   0.950  0.951
#     4        0.874   0.872  0.867
#     5        0.912   0.895  0.883
#     6        0.918   0.917  0.920
#     7        0.926   0.931  0.921
#     8        0.918   0.909  0.913
#     9        0.912   0.905  0.912
#    10        0.911   0.887  0.886
#    11        0.872   0.850  0.849
#    12        0.921   0.913  0.913
mlr::measureAUC(
    probabilities = meta_testing_preds$xgb_pred_mask,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9186929
#
mlr::measureAUC(
    probabilities = meta_testing_preds$maxent_pred_mask,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9151032
#
mlr::measureAUC(
    probabilities = meta_testing_preds$rf_pred_meta,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.8799412
mlr::measureAUC(
    probabilities = meta_testing_preds$rf_pred_meta_mask_0,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9298302
mlr::measureAUC(
    probabilities = meta_testing_preds$rf_pred_meta_mask_1,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.925885


# DALEX -------------------------------------------------------------------

rf_pred_meta_mask_1_explainer = DALEX::explain(
    meta_learner_rf_mask,
    data =  meta_testing_preds %>%
        select(xgb_pred_mask, maxent_pred_mask),
    y = as.numeric(meta_testing_preds$pr_ab)
)
# rf_pred_meta_mask_0_explainer = DALEX::explain(
#     meta_learner_rf_mask,
#     data =  meta_testing_preds %>%
#         select(xgb_pred_mask, maxent_pred_mask),
#     y = as.numeric(meta_testing_preds$pr_ab)
# )

rf_pred_meta_explainer = DALEX::explain(
    meta_learner_rf,
    data =  meta_testing_preds %>%
        select(xgb_pred, maxent_pred),
    y = as.numeric(meta_testing_preds$pr_ab)
)

# rf_pred_meta_mask_1_explainer %>% DALEX::variable_importance() %>% plot()
# rf_pred_meta_mask_0_explainer %>% DALEX::model_profile() %>% plot()
rf_pred_meta_mask_1_explainer %>% DALEX::model_profile() %>% plot()
rf_pred_meta_explainer %>% DALEX::model_profile() %>% plot()
