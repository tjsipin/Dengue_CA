library(tidyverse)
library(tidymodels)
library(rsample)
library(caret)
library(ranger)
library(mlrMBO)
library(lhs)


training = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
training_clusters = readRDS("Data/1_DataProcessing/modeling/training_L1_cv_2016_2022_F_v2.rds")$spatial_clusters
training_blocks = readRDS("Data/1_DataProcessing/modeling/training_L1_cv_2016_2022_F_v2.rds")$spatial_blocks

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
    mutate(L_level = 1)

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
    mutate(L_level = 1)

testing = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2022_F_v2.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 


modelPerMonth = function(m, this_formula, master.path, this_training, training_type){
    dir.path = paste0(master.path, training_type, "/", m, "/")
    dir.create(dir.path, recursive=T)
    
    ## Step 1: Assign training_m and testing_m
        ## m in 1:12
        ## Train model_m on all training data except for those whose Month==m
        ## Predict (and validate) on training data whose Month==m
        ## testing_m goes on to L2 model as its training data (training_L2)
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
    
    ## Step 2: Define hyperparameter design matrix
    ## Set parameter set boundaries
    param_set = makeParamSet(
        makeIntegerParam("tree_depth", lower = 1, upper = 15),
        makeNumericParam("learn_rate", lower = 0.01, upper = 0.3),
        makeIntegerParam("mtry", lower = 1, upper = 6), # number of covariates = 6
        makeIntegerParam("min_n", lower = 1, upper = 40),
        makeNumericParam("loss_reduction", lower = 0, upper = 15)
    )
    
    ## Create design matrix with n sets
    des <- generateDesign(n = 10L * 4L, param_set, fun = randomLHS)
    
    ## Step 3: Get metrics across spatial folds for month m per hyperparameter set
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
            print(des_row)
            # Get parameter set params
            params = des[des_row, ] 
            # Execute XGBoost function with each parameter in the set params
            this_xgb = xgbFunc(
                tree_depth = params$tree_depth, learn_rate = params$learn_rate,
                mtry = params$mtry, min_n = params$min_n,
                loss_reduction = params$loss_reduction,
                fold = fold, m = m, this_formula = this_formula, param_set = des_row,
                training_m = training_m, testing_m = testing_m
            )
            ## Populate in_fold with this hyperparameter set selection
            in_fold[[des_row]] = this_xgb
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
    
    saveRDS(all_folds, paste0(dir.path, "xgb_all_folds_metrics_hurdle.rds"))
    all_folds = readRDS(paste0(dir.path, "xgb_all_folds_metrics_hurdle.rds"))
    
    # With all hyperparameter sets in des, choose the one that maximizes performance for the entire month aggregated by spatial fold
    best_hyperparam_setnumber = all_folds %>% 
        bind_rows() %>% 
        distinct() %>% 
        group_by(param_set) %>% 
        summarize(auc = mean(auc, na.rm = T)) %>% 
        # filter(abs(mbe) == min(abs(mbe))) %>%
        ungroup() %>% 
        arrange(desc(auc)) %>% 
        filter(row_number() == 1)
    
    best_hyperparams = des %>% 
        filter(row_number() == best_hyperparam_setnumber$param_set)
    
    # Define final XGBoost model for month m trained on training_m
    xgb_final_model = boost_tree(
        mode = "classification",
        engine = "xgboost",
        tree_depth = best_hyperparams$tree_depth, 
        learn_rate = best_hyperparams$learn_rate,
        mtry = best_hyperparams$mtry,
        min_n = best_hyperparams$min_n,
        loss_reduction = best_hyperparams$loss_reduction
    ) %>% 
        fit(this_formula, data = training_m)
    
    saveRDS(xgb_final_model, paste0(dir.path, "xgb_final_model.rds"))
    
    # Get prediction set to feed into L2 model
    testing_m_pred = testing_m %>% 
        mutate(
            pred_prob = predict(xgb_final_model, ., 'prob')[, 2][[1]]
        )
    measureAUC(
        truth = testing_m_pred$pr_ab,
        probabilities = testing_m_pred$pred_prob,
        positive = '1'
    ) 
    saveRDS(testing_m_pred, paste0(dir.path, "testing_m_pred.rds"))
}
xgbFunc <- function(tree_depth, learn_rate, mtry, min_n, loss_reduction, fold, m, this_formula, param_set, training_m, testing_m){ # returns vector
    
    
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
    
    xgb.model = boost_tree(
        mode = "classification",
        engine = "xgboost",
        tree_depth = tree_depth,
        learn_rate = learn_rate,
        mtry = mtry,
        min_n = min_n,
        loss_reduction = loss_reduction
    ) %>%
        fit(this_formula, data = cv_analysis)
    
    
    xgb.pred.train = predict(xgb.model, cv_analysis, 'prob')[,2][[1]]
    xgb.pred.test = predict(xgb.model, cv_assessment, 'prob')[,2][[1]]
    
    best_auc.thresh = pROC::coords(
        pROC::roc(
            response = cv_analysis$pr_ab,
            predictor= xgb.pred.train,
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]
    
    if(length(best_auc.thresh) > 1){
        return(NULL)
    }
    
    xgb.pred.test.thresh = ifelse(xgb.pred.test > best_auc.thresh, 1, 0) %>% factor(levels = c(0, 1))
    
    truth = factor(cv_assessment$pr_ab, levels = c(0, 1))
    response = xgb.pred.test.thresh
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
        tree_depth = tree_depth,
        learn_rate = learn_rate,
        mtry = mtry,
        min_n = min_n,
        loss_reduction = loss_reduction
    )
}
