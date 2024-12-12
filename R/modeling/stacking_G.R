# XGBoost model variables: tws, ppt, mean temp
# MaxENT model variables: tws, ppt, mean temp

rm(list = ls())
library(tidyverse)
library(stacks)
library(maxnet)
tidymodels::tidymodels_prefer()

#Read in data
source("R/1_DataProcessing/E/modeling/xgb_modeling_functions_E.R")

output.path = "Data/1_DataProcessing/modeling/stacking_E/"
dir.create(output.path)

xgb_model = readRDS("Data/1_DataProcessing/modeling/xgb_E/tws_ppt_meantemp/best_model.rds")

maxent_model = readRDS("Data/1_DataProcessing/modeling/maxent_E/maxent_model.rds")

#Create random forest meta-learner
training_preds = training %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
        maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    )

ggplot(training_preds) + 
    geom_boxplot(
        aes(
            y = xgb_pred
        )
    ) + 
    facet_wrap(~lc_type)

xgb_mask_thresh = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = training_preds$pr_ab, 
            predictor= training_preds$xgb_pred_mask, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

maxent_mask_thresh = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = training_preds$pr_ab, 
            predictor= training_preds$maxent_pred_mask, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

testing_preds = testing %>% 
    mutate(
        xgb_pred = predict(xgb_model, ., 'prob')[,2][[1]],
        maxent_pred = predict(maxent_model, ., type = "prob")[,2][[1]]
    ) %>% 
    mutate(
        xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0),
        maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0)
    ) %>% 
    mutate(
        xgb_pred_mask_thresh = ifelse(xgb_pred_mask >= xgb_mask_thresh, 1, 0) %>% as.factor(),
        maxent_pred_mask_thresh = ifelse(maxent_pred_mask >= maxent_mask_thresh, 1, 0) %>% as.factor()
    )



meta_learner_rf_both_mask = rand_forest(
    mode = "classification",
    engine = "ranger"
) %>% fit(
    formula = formula("pr_ab ~ xgb_pred_mask + maxent_pred_mask"),
    data = training_preds
)

saveRDS(meta_learner_rf_both_mask, paste0(output.path, "/meta_learner_rf_both_mask.rds"))
meta_learner_rf_both_mask = readRDS(paste0(output.path, "/meta_learner_rf_both_mask.rds"))

meta_training_preds = training_preds %>% 
    mutate(
        rf_pred_meta = predict(meta_learner_rf_both_mask, ., 'prob')$.pred_1
        # glm_pred_meta = predict(meta_learner_glm, .),
        # xgb_pred_meta = predict(meta_learner_xgb, ., 'prob')$.pred_1,
        # svm_pred_meta = predict(meta_learner_svm, ., 'prob')$.pred_1,
        # svmrbf_pred_meta = predict(meta_learner_svmrbf, ., 'prob')$.pred_1
    )

meta_thresh = # create auc threshold for testing
    pROC::coords(
        pROC::roc(
            response = meta_training_preds$pr_ab, 
            predictor= meta_training_preds$rf_pred_meta, 
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]

meta_testing_preds = testing_preds %>% 
    mutate(
        rf_pred_meta = predict(meta_learner_rf_both_mask, ., 'prob')$.pred_1,
        # glm_pred_meta = predict(meta_learner_glm, .),
        # xgb_pred_meta = predict(meta_learner_xgb, ., 'prob')$.pred_1,
        # svm_pred_meta = predict(meta_learner_svm, ., 'prob')$.pred_1,
        # svmrbf_pred_meta = predict(meta_learner_svmrbf, ., 'prob')$.pred_1
    ) %>% 
    mutate(rf_pred_meta_thresh = ifelse(rf_pred_meta >= meta_thresh, 1, 0) %>% as.factor())

mlr::measureAUC(
    probabilities = meta_testing_preds$xgb_pred_mask,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9394003

mlr::measureAUC(
    probabilities = meta_testing_preds$maxent_pred_mask,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9144676

mlr::measureAUC(
    probabilities = meta_testing_preds$rf_pred_meta,
    truth = meta_testing_preds$pr_ab,
    positive = '1'
) # 0.9332317


# Confusion matrices ------------------------------------------------------

confusionMatrix(
    data = testing_preds$xgb_pred_mask_thresh,
    reference = testing_preds$pr_ab,
    positive = "1"
)

confusionMatrix(
    data = testing_preds$maxent_pred_mask_thresh,
    reference = testing_preds$pr_ab,
    positive = "1"
)

confusionMatrix(
    data = meta_testing_preds$rf_pred_meta_thresh,
    reference = meta_testing_preds$pr_ab,
    positive = "1"
)
