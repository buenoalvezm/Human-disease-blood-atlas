
# Functions for visualization
library(limma)
library(tidymodels)
library(themis)

# Run differential expression using limma
de_limma_disease <-
  function(data_wide, 
           metadata,
           disease,
           correct = T,
           cutoff = 0) {
    
    # Select current disease
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>% 
      rename(Group = Disease) %>% 
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
    
    
    if(correct == T) {
     dat <- 
       dat |> 
       filter(!is.na(Sex),
              !is.na(Age))
    } else {
     
      dat <- dat
    }
    

    # Design a model - add Group, and Sex, Age, BMI
    if(correct == T ) {
      design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age ) 
      colnames(design) <- c("control", "case",  "Sex", "Age") 
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    # Fit linear model to each protein assay
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -BMI, -Group)  %>% 
      column_to_rownames("DAid") %>% 
      t()
    
    fit <- lmFit(dat_fit, design = design,  method = "robust", maxit = 10000)
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit)
    
    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n = nrow(ebays_fit$p.value), 
               adjust.method = "fdr", 
               confint = TRUE)
    
    DE_res <- 
      DE_results %>% 
      as_tibble(rownames = "Assay") %>% 
      mutate(Disease = disease,
             sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
                             adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
                             T ~ "not significant"))
    
    return(DE_res)
  }

## Binary classification against healthy
disease_against_healthy <-  
  function(disease, 
           split_train, 
           split_test) {
    
    # Prepare data - make custom split for current disease (+ Healthy samples) based on master split (including all diseases)
    training_dat <- 
      split_train |> 
      filter(Disease %in% c(disease, "Healthy")) |> 
      mutate(Disease = case_when(Disease == disease ~ paste0("0_", disease),
                                 Disease == "Healthy" ~ "1_Healthy")) |> 
      mutate(Disease = factor(Disease))
    
    testing_dat <- 
      split_test |> 
      filter(Disease %in% c(disease, "Healthy")) |> 
      mutate(Disease = case_when(Disease == disease ~ paste0("0_", disease),
                                 Disease == "Healthy" ~ "1_Healthy"))|> 
      mutate(Disease = factor(Disease))
    
    
    # For male and female diseases, filter only male and female samples respectively 
    if(disease %in% female_diseases) {
      
      training_dat <- 
        training_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "F") |> 
        select(-Sex)
      
      testing_dat <- 
        testing_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "F") |> 
        select(-Sex)
    } else if (disease %in% male_diseases) {
      
      training_dat <- 
        training_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "M") |> 
        select(-Sex)
      
      testing_dat <- 
        testing_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "M") |> 
        select(-Sex)
    } else {
      training_dat <- training_dat
      testing_dat <- testing_dat
      
    }
    
    cancer_split_custom <- make_splits(training_dat, testing_dat)
    
    # Recipe with ML steps
    disease_recipe <- 
      recipe(Disease ~ ., data = training_dat) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) |> 
      step_downsample(Disease)
    
    # LASSO model specifications
    glmnet_specs <- 
      logistic_reg() |> 
      set_mode("classification") |> 
      set_engine("glmnet") |> 
      set_args(penalty = tune(), 
               mixture = 1) 
    
    # ML workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(disease_recipe) |> 
      add_model(glmnet_specs) 
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 10)
    
    # Define the resamples (CV)
    set.seed(213)
    cancer_rs <- vfold_cv(training_dat, v = 5, strata = Disease)
    
    # Define the evaluation metrics (add brier)
    eval_metrics <- metric_set(roc_auc)
    
    # Define control_grid
    set.seed(213)
    ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything") # look at extract = identity
    
    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = cancer_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics
      )
    #autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    # Select best hyperparameter
    best_glmnet <- 
      select_best(glmnet_res, metric = "roc_auc") |> 
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <- 
      glmnet_wflow |>  
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow, cancer_split_custom, metrics = eval_metrics) 
    
    # Extract model performance
    performance <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      select(-.config, -.estimator)
    
    glmnet_auc <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      filter(.metric == "roc_auc") |> 
      pull(.estimate) |> 
      round(2)
    
    # Extract protein importance
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0) |> 
      select(-penalty)
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Confusion matrix
    cm <-
      predictions |>
      mutate(pred = ifelse(.pred_1_Healthy > 0.5, "Healthy", disease),
             Disease = ifelse(Disease == "1_Healthy", "Healthy", disease)) |>
      mutate(Disease = factor(Disease, levels = c(disease, "Healthy")),
             pred = factor(pred, levels = c(disease, "Healthy"))) |> 
      conf_mat(Disease, pred)
    
    # ROC curve
    roc <- 
      predictions |>
      roc_curve(truth = Disease, paste0(".pred_0_", disease)) 
    
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "confusion_matrix" = cm,
                "roc_curve" = roc, 
                "important_proteins" = important_proteins))
  }

## Multiclassification against class
disease_against_class <- 
  function(class, 
           split_train, 
           split_test) {
    
    
    disease_train <- 
      split_train |> 
      left_join(resource_meta |> 
                  select(DAid, Class), by = "DAid") |>
      filter(Class == class,
             Disease %in% include_diseases) |> 
      select(-Class)
    
    disease_test <-   
      split_test |> 
      left_join(resource_meta |> 
                  select(DAid, Class), by = "DAid") |>
      filter(Class == class,
             Disease %in% include_diseases) |> 
      select(-Class)
    
    disease_split <- make_splits(disease_train, disease_test)
    
    
    ### Define general recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = disease_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) 
    
    # Generate resamples
    set.seed(213)
    disease_rs <- vfold_cv(disease_train, v = 5, strata = Disease)
    
    # Define evaluation metrics for all workflows
    eval_metrics <- metric_set(roc_auc)
    
    # Define control grid
    set.seed(213)
    ctrl <- control_grid(verbose = TRUE, 
                         allow_par = TRUE,
                         save_pred = TRUE, 
                         parallel_over = "everything") 
    
    # Tidymodels lasso multiclassification recipe
    glmnet_lasso_specs <-
      multinom_reg() |>
      set_mode("classification") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # Set up lasso workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs) 
    
    # Define hyperparameter tuning grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 15)
    
    # Hyperparameter tuning
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = disease_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics)
    
    # Explore results and select the best performing hyperparameter combination
    autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    best_glmnet <- 
      glmnet_res |> 
      select_best("roc_auc")
    
    # Final fit
    final_glmnet <- 
      glmnet_wflow |> 
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- 
      last_fit(final_glmnet, disease_split)
    
    # Explore performance
    final_predictions <- 
      final_glmnet_fit |> 
      collect_predictions() |> 
      mutate(DAid = disease_test$DAid) |> 
      relocate(DAid) |> 
      select(-id)
    
    final_metrics <- 
      final_glmnet_fit |> 
      collect_metrics()
    
    # ROC
    dat <-
      disease_test %>% 
      select(DAid, Disease) %>% 
      mutate(value = 1) %>% 
      spread(Disease,value, fill= 0) 
    
    true_dat <- 
      dat %>% 
      set_names(paste(names(dat), "_true", sep = "")) %>%
      rename(DAid = `DAid_true`)
    
    dat_prob <- 
      final_predictions %>% 
      rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
      select(-".row", -"class", -"Disease", -".config") 
    
    prob_data <- 
      dat_prob %>% 
      set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
      rename(DAid = DAid_pred_glmnet)
    
    final_df <- 
      true_dat %>% 
      left_join(prob_data, by = "DAid") %>% 
      select(-DAid) |> 
      as.data.frame()
    
    roc_res <- multi_roc(final_df, force_diag=T)
    plot_roc_df <- plot_roc_data(roc_res)
    
    roc_dat <- 
      plot_roc_df %>%
      filter(!Group %in% c("Macro","Micro")) %>% 
      mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
      arrange(-AUC)
    
    roc_plot <- 
      roc_dat %>% 
      arrange(Performance) %>% 
      ggplot(aes(x = 1-Specificity, y=Sensitivity)) +
      geom_path(size=1, show.legend = F, color = pal_class[class]) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey",
                   linetype = 'dotdash',
                   show.legend = F) +
      geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
      theme_hpa() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 1)) +
      facet_wrap(~Group, nrow = 5) 
    
    # Confusion matrix
    cm <- 
      final_glmnet_fit %>% 
      collect_predictions() %>%
      conf_mat(truth = Disease, estimate = .pred_class) 
    
    cm_plot <- 
      cm |>
      autoplot(type = "heatmap") +
      theme_hpa(angled = T)
    
    # Performance scores  
    cm_2 <- confusionMatrix(table(final_predictions$.pred_class, disease_test$Disease))
    
    performance <- 
      as_tibble(cm_2$byClass, rownames = "Class") %>% 
      mutate(Class = gsub("Class: ", "", Class))
    
    # Probabilities
    
    # Important proteins
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0)
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_fit" = final_glmnet_fit,
                "predictions" = final_predictions,
                "final_metrics" = final_metrics,
                "performance" = performance,
                "confusion_matrix" = cm,
                "confusion_matrix_plot" = cm_plot,
                "roc_curve" = roc_dat,
                "roc_plot" = roc_plot,
                "important_proteins" = important_proteins))
    
  }

## Multiclassification against all other diseases
disease_against_all<-  
  function(split_train, 
           split_test) {                    
    
    disease_train <- 
      split_train |> 
      filter(Disease %in% include_diseases)
    
    disease_test <-   
      split_test |>
      filter(Disease %in% include_diseases) 
    
    disease_split <- make_splits(disease_train, disease_test)
    
    
    ### Define general recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = disease_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) 
    
    # Generate resamples
    set.seed(213)
    disease_rs <- vfold_cv(disease_train, v = 5, strata = Disease)
    
    # Define evaluation metrics for all workflows
    eval_metrics <- metric_set(roc_auc)
    
    # Define control grid
    set.seed(213)
    ctrl <- control_grid(verbose = TRUE, 
                         allow_par = TRUE,
                         save_pred = TRUE, 
                         parallel_over = "everything") 
    
    # Tidymodels lasso multiclassification recipe
    glmnet_lasso_specs <-
      multinom_reg() |>
      set_mode("classification") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # Set up lasso workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs) 
    
    # Define hyperparameter tuning grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 15)
    
    # Hyperparameter tuning
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = disease_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics)
    
    # Explore results and select the best performing hyperparameter combination
    autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    best_glmnet <- 
      glmnet_res |> 
      select_best("roc_auc")
    
    # Final fit
    final_glmnet <- 
      glmnet_wflow |> 
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- 
      last_fit(final_glmnet, disease_split)
    
    # Explore performance
    final_predictions <- 
      final_glmnet_fit |> 
      collect_predictions() |> 
      mutate(DAid = disease_test$DAid) |> 
      relocate(DAid) |> 
      select(-id)
    
    final_metrics <- 
      final_glmnet_fit |> 
      collect_metrics()
    
    # ROC
    dat <-
      disease_test %>% 
      select(DAid, Disease) %>% 
      mutate(value = 1) %>% 
      spread(Disease,value, fill= 0) 
    
    true_dat <- 
      dat %>% 
      set_names(paste(names(dat), "_true", sep = "")) %>%
      rename(DAid = `DAid_true`)
    
    dat_prob <- 
      final_predictions %>% 
      rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
      select(-".row", -"class", -"Disease", -".config") 
    
    prob_data <- 
      dat_prob %>% 
      set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
      rename(DAid = DAid_pred_glmnet)
    
    final_df <- 
      true_dat %>% 
      left_join(prob_data, by = "DAid") %>% 
      select(-DAid) |> 
      as.data.frame()
    
    roc_res <- multi_roc(final_df, force_diag=T)
    plot_roc_df <- plot_roc_data(roc_res)
    
    roc_dat <- 
      plot_roc_df %>%
      filter(!Group %in% c("Macro","Micro")) %>% 
      mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
      arrange(-AUC)
    
    roc_plot <- 
      roc_dat %>% 
      left_join(resource_meta |> 
                  distinct(Disease, Class), by = c("Group" = "Disease")) |> 
     # mutate(Group = factor(Group, levels = disease_class_order_roc)) |> 
      ggplot(aes(x = 1-Specificity, y=Sensitivity, color= Class)) +
      geom_path(size=1, show.legend = F) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey",
                   linetype = 'dotdash',
                   show.legend = F) +
      geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
      theme_hpa() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 1)) +
      scale_color_manual(values = pal_class) +
      facet_wrap(~Group, nrow = 5) +
      coord_fixed()
    
    # Confusion matrix
    cm <- 
      final_glmnet_fit %>% 
      collect_predictions() %>%
      conf_mat(truth = Disease, estimate = .pred_class) 
    
    cm_plot <- 
      cm |>
      autoplot(type = "heatmap") +
      theme_hpa(angled = T)
    
    # Performance scores  
    cm_2 <- confusionMatrix(table(final_predictions$.pred_class, disease_test$Disease))
    
    performance <- 
      as_tibble(cm_2$byClass, rownames = "Class") %>% 
      mutate(Class = gsub("Class: ", "", Class))
    
    # Important proteins
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0)
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = final_glmnet,
                "final_fit" = final_glmnet_fit,
                "predictions" = final_predictions,
                "final_metrics" = final_metrics,
                "performance" = performance,
                "confusion_matrix" = cm,
                "confusion_matrix_plot" = cm_plot,
                "roc_curve" = roc_dat,
                "roc_plot" = roc_plot,
                "important_proteins" = important_proteins))
    
  }

