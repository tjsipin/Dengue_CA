library(tidyverse)
library(tidymodels)
library(rsample)
library(caret)
library(ranger)
library(mlrMBO)
library(lhs)

## Set parameter set boundaries
param_set = makeParamSet(
    makeIntegerParam("mtry", lower = 1, upper = 2),
    makeIntegerParam("trees", lower = 250L, upper = 2500L),
    makeIntegerParam("min_n", lower = 5L, upper = 20L)
)

## Create design matrix with n sets
# DO NOT RUN
# des <- generateDesign(n = 20L, param_set, fun = randomLHS)
# saveRDS(des, "Data/1_DataProcessing/modeling/stacking_H/des.rds")
des = readRDS("Data/1_DataProcessing/modeling/stacking_H/des.rds")


modelPerMonth = function(m, this_formula, master.path, this_training, training_type){
    dir.path = paste0(master.path, training_type, "/", m, "/")
    dir.create(dir.path, recursive=T)
    
    if(file.exists(paste0(dir.path, "rf_all_folds_metrics_hurdle.rds"))){
        return(NULL)
    }
    
    ## Step 1: Assign training_m and testing_m
    ## m in 1:12
    ## Train model_m on all training data except for those whose Month==m
    ## Predict (and validate) on training data whose Month==m
    training_m = this_training %>%
        filter(Month != m) %>%
        mutate(lc_type = as.factor(lc_type)) %>%
        mutate(pr_ab = as.factor(pr_ab)) %>%
        mutate(
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        )
    
    testing_m = this_training %>%
        filter(Month == m) %>%
        mutate(lc_type = as.factor(lc_type)) %>%
        mutate(pr_ab = as.factor(pr_ab)) %>%
        mutate(
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        )
    
    spat_folds = training_m %>% 
        pull(spat_fold) %>% 
        unique() %>% 
        sort()
    
    
    ## Step 2: Get metrics across spatial folds for month m per hyperparameter set
    ## Make empty parent list to populate all spatial CV folds
    all_folds = list()
    ## Initiate for loop to populate all_folds
    # For each fold in spat_folds do
    t_folds = Sys.time()
    for(fold in spat_folds){ 
        ## Make empty list to populate in fold
        in_fold = list()
        # For each row in design matrix do
        for(des_row in 1:nrow(des)) { 
            message(paste0("Month: ", m))
            message(paste0("fold: ", fold, "/", max(spat_folds)))
            message(paste0("des_row: ", des_row, "/", nrow(des)))
            # Get parameter set params
            params = des[des_row, ] 
            # Execute XGBoost function with each parameter in the set params
            this_rf = rfFunc(
                mtry = params$mtry,
                trees = params$trees,
                min_n = params$min_n,
                fold = fold, m = m, this_formula = this_formula, param_set = des_row, 
                training_m = training_m, testing_m = testing_m
            )
            ## Populate in_fold with this hyperparameter set selection
            in_fold[[des_row]] = this_rf
            message(Sys.time() - t_folds)
        }
        
        # List to DF for all metrics with hyperparameters in fold
        in_fold_df = in_fold %>%
            bind_rows() 
        
        # Populate list index with finished metric-hyperparameter set DF
        all_folds[[fold]] = in_fold_df
    }
    
    # List to DF for all metrics with all hyperparamters for all folds
    all_folds = all_folds %>% bind_rows()
    
    if(nrow(all_folds) == 0){
        return(NULL)
    }
    Sys.time() - t_folds
    
    saveRDS(all_folds, paste0(dir.path, "rf_all_folds_metrics_hurdle.rds"))
}
rfFunc <- function(mtry, trees, min_n, fold, m, this_formula, param_set, training_m, testing_m){ # returns vector
    
    cv_analysis = training_m %>% 
        filter(spat_fold != fold) %>% 
        mutate(lc_type = as.factor(lc_type)) %>% 
        mutate(pr_ab = as.factor(pr_ab)) %>% 
        mutate(
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        ) %>% 
        na.omit()
    cv_assessment = training_m %>% 
        filter(spat_fold == fold) %>% 
        mutate(lc_type = as.factor(lc_type)) %>% 
        mutate(pr_ab = as.factor(pr_ab)) %>% 
        mutate(
            mean_temp = (tmn + tmx)/2,
            tws = aet/pet
        ) %>% 
        na.omit()
    
    if(sum(cv_analysis$pr_ab==1)==0){
        return(NULL)
    } 
    if(length(unique(cv_analysis$pr_ab)) == 0){
        return(NULL)
    }
    
    if(sum(cv_assessment$pr_ab==1)==0){
        return(NULL)
    } 
    if(length(unique(cv_assessment$pr_ab)) == 0){
        return(NULL)
    }
    
    rf_model = rand_forest(
        mode = "classification",
        engine = "randomForest",
        mtry = mtry,
        trees = trees,
        min_n = min_n
    ) %>% fit(
        formula = formula(
            "pr_ab ~ xgb_pred + maxent_pred"
        ),
        data = cv_analysis
    )
    
    
    rf.pred.train = predict(rf_model, cv_analysis, 'prob')[,2][[1]]
    rf.pred.test = predict(rf_model, cv_assessment, 'prob')[,2][[1]]
    
    best_auc.thresh = pROC::coords(
        pROC::roc(
            response = cv_analysis$pr_ab,
            predictor= rf.pred.train,
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]
    
    if(length(best_auc.thresh) > 1){
        return(NULL)
    }
    
    rf.pred.test.thresh = ifelse(rf.pred.test > best_auc.thresh, 1, 0) %>% factor(levels = c(0, 1))
    
    truth = factor(cv_assessment$pr_ab, levels = c(0, 1))
    response = rf.pred.test.thresh
    if(length(truth) != length(response)){
        return(NULL)
    }
    if(length(unique(truth)) < 2){
        return(NULL)
    }
    
    test.bal = mlr::measureBAC(
        truth = truth,
        response = response
    )
    
    test.auc = mlr::measureAUC(
        truth = truth,
        probabilities = response,
        positive = '1'
    )
    
    list(
        test_fold = fold,
        thresh = best_auc.thresh,
        auc = test.auc,
        bal = test.bal,
        param_set = param_set,
        mtry = mtry,
        trees = trees,
        min_n = min_n
    )
}