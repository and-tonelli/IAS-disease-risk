#### XGB ####
library(tidymodels)
library(yardstick)
require(parsnip)
require(xgboost)

DatasetMainModel

Prediction_df_final_ready <- read_csv("Data/Prediction_df_final_ready.csv")
Prediction_df_final_ready2 <- Prediction_df_final_ready %>% 
  mutate(log_sum_cit = median(Prediction_df_final_ready$log_sum_cit))

# Ten rounds of nested cross validation
for (s in round(runif(10, 1, 999), 0)){
  
  # 1. Split into training and test sets
  set.seed(s)
  DatasetMainModel$vir_bin <- fct_relevel(DatasetMainModel$vir_bin, "1")
  data_split <- initial_split(DatasetMainModel, prop = 0.8, strata = vir_bin)
  train_data <- training(data_split)
  test_data <- testing(data_split)
  
  # 2. 10-fold cross-validation on training set
  folds <- vfold_cv(train_data, v = 10, strata = vir_bin)
  
  # 3. Define recipe
  xgb_recipe <- recipe(vir_bin ~ trait_sim_gow + foraging_sim + phylo_sim + log_sum_cit + Overlapping_Cells_log,
                       data = train_data)
  
  # 4. Define model and hyperparameter ranges
  xgb_spec <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),  # gamma
    sample_size = tune(),
    mtry = tune()
  ) %>%
    set_mode("classification") %>%
    set_engine("xgboost", scale_pos_weight = tune())
  
  # 5. Define workflow
  xgb_workflow <- workflow() %>%
    add_model(xgb_spec) %>%
    add_recipe(xgb_recipe)
  
  # 6. Define tuning grid
  set.seed(s)
  
  xgb_param <- parameters(xgb_spec) %>%
    
    update(
      mtry = finalize(mtry(), train_data[c(5:7, 12, 15)]),
      scale_pos_weight = scale_pos_weight(range = c(0.8, 3))
    )
  
  # For parallel running
  library(doParallel)
  cl <- makePSOCKcluster(parallel::detectCores() - 2)
  registerDoParallel(cl)
  
  # 7. Bayesian optimization
  set.seed(123)
  xgb_bayes <- tune_bayes(
    xgb_workflow,
    resamples = folds,
    param_info = xgb_param,
    initial = 10,
    iter = 20,      # <- iterations of Bayesian search
    metrics = metric_set(recall, roc_auc, precision, specificity),
    control = control_bayes(verbose = TRUE, no_improve = 5, save_pred = TRUE)
  )
  
  stopCluster(cl)
  
  xgb_bayes %>% pluck(".metrics")
  
  # 8. Select best hparameters by recall
  best_params <- select_best(xgb_bayes, metric = c("recall"))
  
  # 9. Finalize workflow with best params
  final_workflow <- finalize_workflow(xgb_workflow, best_params)
  
  # 10. Fit final model on full training data
  final_fit <- final_workflow %>%
    last_fit(split = data_split, metrics = metric_set(roc_auc, precision, recall, specificity))
  
  # 11. Estimate performance on test set
  test_metrics <- collect_metrics(final_fit)
  assign(paste0("test_metrics_", s), test_metrics)
  
  # 12. Fit model on full native network
  final_fit <- final_workflow %>% 
    fit(data = DatasetMainModel)
  
  # Save individual models
  assign(paste0("final_fit_", s), final_fit)
  
  # 13a. Predict alien network
  preds <- predict(final_fit, Prediction_df_final_ready, type = "prob")
  
  Prediction_df <- Prediction_df_final_ready
  
  Prediction_df$preds <- preds$.pred_1
  Prediction_df %<>% 
    mutate(pred_bin = ifelse(preds > 0.5, 1, 0))
  
  # Save predictions
  assign(paste0("PredictionRegular_", s), Prediction_df)
  
  # 13b. Predict alien network with median citations
  preds2 <- predict(final_fit, Prediction_df_final_ready2, type = "prob")
  
  Prediction_df2 <- Prediction_df_final_ready2
  
  Prediction_df2$preds <- preds2$.pred_1
  Prediction_df2 %<>% 
    mutate(pred_bin = ifelse(preds > 0.5, 1, 0))
  
  # Save predictions
  assign(paste0("PredictionMedianCit_", s), Prediction_df2)
}

# Check test set metrics
grep("test_metrics_", names(.GlobalEnv), value=TRUE) %>% 
  lapply(., get) %>% 
  # lapply(., mutate(TSS = .$.estimate[.$.metric == "recall"]+.$.estimate[.$.metric == "specificity"]-1)) %>% 
  do.call(rbind, .) %>% 
  group_by(.metric) %>% 
  summarise(mean = mean(.estimate),
            sd = sd(.estimate))

# Check predictions
grep("PredictionRegular_", names(.GlobalEnv), value=TRUE) %>% 
  lapply(., get) %>% 
  do.call(rbind, .) %>% 
  # filter(vir_bin == 0) %>% 
  group_by(Species1, Species2, vir_bin) %>% 
  summarise(mean_prob = mean(preds),
            sd_prob = sd(preds),
            pred_bin = ifelse(mean_prob >= 0.5, 1, 0)) -> FinalPredictions

combined_df <- mapply(function(x, i) {
  x$id <- i
  x
}, lapply(grep("PredictionRegular_", names(.GlobalEnv), value = TRUE), get), 
seq_along(grep("PredictionRegular_", names(.GlobalEnv), value = TRUE)), 
SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

saveRDS(FinalPredictions, "Data/FinalPreds.rds")

# Check predictions (citations held at median)
grep("PredictionMedianCit_", names(.GlobalEnv), value=TRUE) %>% 
  lapply(., get) %>% 
  do.call(rbind, .) %>% 
  # filter(vir_bin == 0) %>% 
  group_by(Species1, Species2, vir_bin) %>% 
  summarise(mean_prob = mean(preds),
            sd_prob = sd(preds),
            pred_bin = ifelse(mean_prob >= 0.5, 1, 0)) -> FinalPredictionsMedianCit

combined_df2 <- mapply(function(x, i) {
  x$id <- i
  x
}, lapply(grep("PredictionMedianCit_", names(.GlobalEnv), value = TRUE), get), 
seq_along(grep("PredictionMedianCit_", names(.GlobalEnv), value = TRUE)), 
SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

