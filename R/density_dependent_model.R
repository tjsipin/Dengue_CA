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
# paramter `bw`
bw = 500

formula = formula("pr_ab ~ ppt + mean_temp + tws")

ppt_density = training_L1 %>% 
    select(pr_ab, ppt) %>% 
    mutate(pr_ab = as.integer(as.character(pr_ab))) %>% 
    mutate(
        bucket = cut(ppt, breaks = quantile(unique(training_L1$ppt), probs = seq(0, 1, length.out = bw + 1)), include.lowest = T) # parameter: bw
    ) %>% 
    group_by(bucket) %>% 
    summarize(pr_ab = mean(pr_ab)) %>% 
    ungroup() %>% 
    mutate(
        bucket = as.character(bucket) %>% 
            str_replace("\\(", "") %>% 
            str_replace("\\)", "") %>% 
            str_replace("\\[", "") %>% 
            str_replace("\\]", "") 
    ) %>% 
    mutate(
        min = str_split_i(bucket, pattern = ",",1) %>% as.numeric(),
        max = str_split_i(bucket, pattern = ",",2) %>% as.numeric()
    ) %>% 
    select(pr_ab, min, max) %>% 
    arrange(min) %>% 
    as.data.table()

mean_temp_density = training_L1 %>% 
    select(pr_ab, mean_temp) %>% 
    mutate(pr_ab = as.integer(as.character(pr_ab))) %>% 
    mutate(
        bucket = cut(mean_temp, breaks = quantile(unique(training_L1$mean_temp), probs = seq(0, 1, length.out = bw + 1)), include.lowest = T) # parameter: bw
    ) %>% 
    group_by(bucket) %>% 
    summarize(pr_ab = mean(pr_ab)) %>% 
    ungroup() %>% 
    mutate(
        bucket = as.character(bucket) %>% 
            str_replace("\\(", "") %>% 
            str_replace("\\)", "") %>% 
            str_replace("\\[", "") %>% 
            str_replace("\\]", "") 
    ) %>% 
    mutate(
        min = str_split_i(bucket, pattern = ",",1) %>% as.numeric(),
        max = str_split_i(bucket, pattern = ",",2) %>% as.numeric()
    ) %>% 
    select(pr_ab, min, max) %>% 
    arrange(min) %>% 
    as.data.table()

# Convert data frames to data tables
eval_dt <- as.data.table(training_L2)

# Perform a non-equi join using data.table syntax
ppt_pred <- eval_dt[ppt_density, on = .(ppt >= min, ppt <= max), 
                    ppt_pred := i.pr_ab] %>% 
    tibble() %>% 
    mutate(id = row_number()) %>% 
    select(id, pr_ab, ppt_pred)

tmean_pred <- eval_dt[mean_temp_density, on = .(mean_temp >= min, mean_temp <= max), 
                    tmean_pred := i.pr_ab] %>% 
    tibble() %>% 
    mutate(id = row_number()) %>% 
    select(id, pr_ab, tmean_pred)

res = ppt_pred %>% 
    full_join(tmean_pred) 

# View the result
print(result)
measureAUC(result$presence_rate, result$pr_ab, positive = 1)

## With parameter `bw`


kd = MASS::kde2d(training_L1$mean_temp, training_L1$ppt)
