
savepath <- 
  function(savename) { 
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_folder <- 
  function(folder, savename) { 
    result_folder <- paste0("results/", Sys.Date(), "/",folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_data <- 
  function(folder, savename) { 
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_results <- 
  function(folder, savename) { 
    
    dir.create("results/", showWarnings = FALSE)
    result_folder <- paste0("results/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    return(savename)
    
  }

do_umap <- 
  function(wide_data, 
           seed = 42, 
           n_neighbors = 15,
           n_components = 2, 
           ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(uwot))
    
    set.seed(seed)
    
    umap_res <- 
      umap(wide_data, 
           n_neighbors = n_neighbors,
           n_components = n_components,
           ...) %>% 
      as.data.frame()
    
    rownames(umap_res) <- rownames(wide_data)
    colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
    
    umap_res
  }
impute_values <- 
  function(data, ID, wide_data = F) {
    
    if(wide_data == F) {
      data_wide <- 
        data %>% 
        select(ID, Assay, NPX) %>% 
        spread(Assay,NPX) 
      
    } else {
      data_wide <- 
        data
    }
    
    data_imputed <- 
      data_wide %>% 
      column_to_rownames(ID) %>% 
      as.matrix() %>% 
      t() %>% 
      impute.knn() 
    
    final_data <- 
      data_imputed$data %>% 
      t() %>% 
      as_tibble(rownames = ID)
    
    return(final_data)
    
  }
