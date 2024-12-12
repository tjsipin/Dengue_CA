library(tidyverse)
library(tidymodels)
library(rsample)
library(caret)
library(ranger)
library(mlrMBO)
library(lhs)

training = readRDS("Data/1_DataProcessing/modeling/training_L1_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 
training_clusters = readRDS("Data/1_DataProcessing/modeling/training_L1_cv_2016_2023.rds")$spatial_clusters
training_blocks = readRDS("Data/1_DataProcessing/modeling/training_L1_cv_2016_2023.rds")$spatial_blocks

training_clusters_full = map(
    1:length(training_clusters$id),
    function(s0){
        s = training_clusters$splits[[s0]]
        as = rsample::assessment(s) %>%
            mutate(spat_fold = s0)
        as
    }
) %>% bind_rows() %>%
    sf::st_drop_geometry()

training_blocks_full = map(
    1:length(training_blocks$id),
    function(s0){
        s = training_blocks$splits[[s0]]
        as = rsample::assessment(s) %>%
            mutate(spat_fold = s0)
        as
    }
) %>% bind_rows() %>%
    sf::st_drop_geometry()


testing = readRDS("Data/1_DataProcessing/modeling/training_L2_2016_2023.rds") %>% 
    mutate(lc_type = as.factor(lc_type)) %>% 
    mutate(pr_ab = as.factor(pr_ab)) 

## Set parameter set boundaries
param_set = makeParamSet(
    makeNumericParam("regularization_multiplier", lower = 0.1, upper = 1.5)
)


## Create design matrix with n sets
des <- generateDesign(n = 5L, param_set, fun = randomLHS)


modelPerMonth = function(m, this_formula, master.path, this_training, training_type){
    dir.path = paste0(master.path, training_type, "/", m, "/")
    dir.create(dir.path, recursive=T)
    
    if(file.exists(paste0(dir.path, "testing_m_pred.rds"))){
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
            this_maxent = maxentFunc(
                regularization_multiplier = params,
                fold = fold, m = m, this_formula = this_formula, param_set = des_row, 
                training_m = training_m, testing_m = testing_m
            )
            ## Populate in_fold with this hyperparameter set selection
            in_fold[[des_row]] = this_maxent
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
    
    saveRDS(all_folds, paste0(dir.path, "maxent_all_folds_metrics_hurdle.rds"))
    all_folds = readRDS(paste0(dir.path, "maxent_all_folds_metrics_hurdle.rds"))
    
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
        filter(row_number() == best_hyperparam_setnumber$param_set) %>% 
        pull(regularization_multiplier)
    
    # Define final XGBoost model for month m trained on training_m
    t0 = Sys.time()
    maxent_final_model = maxent(
        mode = "classification",
        engine = "maxnet", 
        feature_classes = "lqph",
        regularization_multiplier = best_hyperparams
    ) %>% 
        fit(this_formula, data = training_m)
    
    saveRDS(maxent_final_model, paste0(dir.path, "maxent_final_model.rds"))
    Sys.time() - t0
    
    # Get prediction set to feed into L2 model
    testing_m_pred = testing_m %>% 
        mutate(
            pred_prob = predict(maxent_final_model, ., 'prob')[, 2][[1]]
        )
    measureAUC(
        truth = testing_m_pred$pr_ab,
        probabilities = testing_m_pred$pred_prob,
        positive = '1'
    ) 
    saveRDS(testing_m_pred, paste0(dir.path, "testing_m_pred.rds"))
}

## Base function to create a MaxENT model
## Input: 
### Regularization multiplier (hyperparameter)
### Spatial fold
### Month
### Formula (universal)
### Parameter set number ID
### Training given month (all training data minus month)
### Testing given month (training data)
## Output:
### Returns vector of performance
maxentFunc <- function(regularization_multiplier, fold, m, this_formula, param_set, training_m, testing_m){ 
    
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
    t0=Sys.time()
    maxent_model = maxent(
        mode = "classification",
        engine = "maxnet", 
        feature_classes = "lqph",
        regularization_multiplier = regularization_multiplier
    ) %>% 
        fit(formula, data = cv_analysis)
    Sys.time()-t0
    
    maxent.pred.train = predict(maxent_model, cv_analysis, 'prob')[,2][[1]]
    maxent.pred.test = predict(maxent_model, cv_assessment, 'prob')[,2][[1]]
    
    best_auc.thresh = pROC::coords(
        pROC::roc(
            response = cv_analysis$pr_ab,
            predictor= maxent.pred.train,
            levels=c(0,1), auc = TRUE
        ), "best", ret = "threshold"
    )[1][[1]]
    
    if(length(best_auc.thresh) > 1){
        return(NULL)
    }
    
    maxent.pred.test.thresh = ifelse(maxent.pred.test > best_auc.thresh, 1, 0) %>% factor(levels = c(0, 1))
    
    truth = factor(cv_assessment$pr_ab, levels = c(0, 1))
    response = maxent.pred.test.thresh
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
        regularization_multiplier = regularization_multiplier
    )
}
