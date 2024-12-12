xgb_model_G = readRDS("Data/1_DataProcessing/modeling/xgb_G/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model_G = readRDS("Data/1_DataProcessing/modeling/maxent_G/tws_ppt_meantemp/maxent_best_model_clusters.rds")
meta_model_G = readRDS("Data/1_DataProcessing/modeling/stacking_G/meta_learner_rf_mask_randomForest.rds")
testing_G = readRDS("Data/1_DataProcessing/modeling/testing_2016_2022_F.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
xgb_predictions_G = testing_G %>% 
    mutate(xgb_pred = predict(xgb_model_G, ., 'prob')$.pred_1) %>% 
    mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))

maxent_predictions_G = testing_G %>% 
    mutate(maxent_pred = predict(maxent_model_G, ., 'prob')$.pred_1) %>% 
    mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))
##Join together both individual models
indiv_predictions_G = full_join(xgb_predictions_G, maxent_predictions_G)
meta_predictions_G = indiv_predictions_G %>% 
    mutate(meta_pred = predict(meta_model_G, ., 'prob')$.pred_1)


# H -----------------------------------------------------------------------

xgb_model_H = readRDS("Data/1_DataProcessing/modeling/xgb_H/tws_ppt_meantemp/xgb_best_model_clusters.rds")
maxent_model_H = readRDS("Data/1_DataProcessing/modeling/maxent_H/tws_ppt_meantemp/maxent_zoe_model.rds")
meta_model_H = readRDS("Data/1_DataProcessing/modeling/stacking_H/meta_learner_rf_mask_randomForest.rds")

testing_H = readRDS("Data/1_DataProcessing/modeling/testing_2016_2023_H.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

xgb_predictions_H = testing_H %>% 
    mutate(xgb_pred = predict(xgb_model_H, ., 'prob')$.pred_1) %>% 
    mutate(xgb_pred_mask = ifelse(lc_type %in% c(2, 13, 14), xgb_pred, 0))
maxent_predictions_H = testing_H %>% 
    mutate(maxent_pred = predict(maxent_model_H, ., 'prob')$.pred_1) %>% 
    mutate(maxent_pred_mask = ifelse(lc_type %in% c(2, 13, 14), maxent_pred, 0))
##Join together both individual models
indiv_predictions_H = full_join(xgb_predictions_H, maxent_predictions_H)
meta_predictions_H = indiv_predictions_H %>% 
    mutate(meta_pred = predict(meta_model_H, ., 'prob')$.pred_1)

meta_predictions_joined = meta_predictions_G %>% 
    select(Month, meta_pred, xgb_pred, maxent_pred) %>% 
    mutate(model = "G") %>% 
    rbind(
        meta_predictions_H %>% 
            select(Month, meta_pred, xgb_pred, maxent_pred) %>% 
            mutate(model = "H")
    )

ggplot() +
    geom_density(
        data = meta_predictions_joined,
        aes(x = meta_pred, color = model)
    ) +
    facet_wrap(
        ~Month
    ) +
    theme_classic()


ggplot() +
    geom_density(
        data = meta_predictions_joined %>% 
            filter(meta_pred > 0),
        aes(x = meta_pred, color = model)
    ) +
    facet_wrap(
        ~Month
    ) +
    theme_classic() +
    labs(title = "Meta predictions per month", subtitle = "(current conditions)", x = "")

ggplot() +
    geom_density(
        data = meta_predictions_joined %>% 
            filter(meta_pred > 0),
        aes(x = xgb_pred, color = model)
    ) +
    facet_wrap(
        ~Month
    ) +
    theme_classic() +
    labs(title = "XGBoost predictions per month", subtitle = "(current conditions)", x = "")

ggplot(
    data = meta_predictions_joined %>% 
        filter(meta_pred > 0)
) +
    geom_density(
        aes(x = maxent_pred, color = model)
    ) +
    facet_wrap(
        ~Month
    ) +
    theme_classic() +
    labs(title = "MaxENT predictions per month", subtitle = "(current conditions)", x = "")


ggplot(
    data = meta_predictions_joined %>% 
        filter(meta_pred > 0)
) +
    geom_boxplot(
        aes(y = maxent_pred, color = model)
    ) +
    facet_wrap(
        ~Month
    ) +
    theme_classic()