
# Functions to perform Olink QC 
generate_qc_params <-
  function(data,
           manifest) {
    
    # Define control assays and samples
    control_sample_types <-  c("SAMPLE_CONTROL", "PLATE_CONTROL", "NEGATIVE_CONTROL")
    control_assay_types <- c("ext_ctrl", "inc_ctrl", "amp_ctrl")
    
    # Extract control samples from manifest    
    control_samples <-
      manifest |>
      filter(Cohort == "CTRL") |>
      pull(DAid)
    
    # Generate Olink metadata file
    olink_meta <-
      data |>
      filter(!AssayType %in% control_assay_types) |>
      distinct(Assay, OlinkID, UniProt, AssayType, Block)
    
    # Map SampleIDs to DA identifiers 
    id_mapping <- 
      data_ht |> 
      filter(SampleType == "SAMPLE") |> 
      distinct(SampleID) |> 
      separate(SampleID, into = c("DAid", "Batch", "Plate"), sep = "-", remove = F) 
    
    # Map data
    mapped_data <- 
      data |>
      left_join(id_mapping, by = "SampleID") |>
      relocate(DAid)
    
    # Calculate LOD
    lod_data <- olink_lod(data = mapped_data, lod_method = "NCLOD")
    
    # Calculate number of control measurements per protein
    n_measurements_control <- 
      lod_data |> 
      filter(DAid %in% control_samples)  |> 
      count(Assay) |> 
      distinct(n) |> 
      pull()
    
    n_under_lod_controls <- 
      lod_data |> 
      filter(DAid %in% control_samples) |> 
      mutate(under_LOD = ifelse(NPX < LOD, "Yes", "No")) |>
      group_by(Assay) |> 
      count(under_LOD) |> 
      arrange(-n)
    
    # Exclude protens below LOD in 50% of control measurements
    exclude_proteins <- 
      n_under_lod_controls |> 
      filter(under_LOD == "Yes") |> 
      ungroup() |> 
      filter(n >= n_measurements_control/2)
    
    # Calculate number of proteins, plates and samples in the dataset
    n_samples_run <-
      data |>
      filter(!SampleType %in% control_sample_types) |>
      distinct(SampleID) |>
      nrow()
    
    n_samples <- 
      mapped_data |> 
      filter(!SampleType %in% control_sample_types) |> 
      distinct(DAid) |> 
      nrow()
    
    n_proteins <-
      olink_meta$Assay |>
      unique() |>
      length()
    
    n_plates <- 
      id_mapping |> 
      distinct(Plate) |> 
      nrow()
    
    # Filtered manifest
    current_ids <- 
      mapped_data |> 
      filter(!SampleType %in% control_sample_types) |> 
      distinct(DAid) 
    
    manifest_filtered <- 
      manifest |> 
      filter(DAid %in% current_ids$DAid) |> 
      mutate(Cohort = paste(Cohort, `PI - sample owner`, sep = "_")) 
    
    # Manifest cohorts
    cohort_number <- 
      manifest_filtered|> 
      distinct(Cohort) |> 
      nrow()
    
    # Disease order
    disease_order <- 
      manifest_filtered |> 
      mutate(Disease = ifelse(Disease == "Healthy", paste(Disease, Cohort, sep = " - "), Disease)) |> 
      distinct(Disease, Cohort, Class) |> 
      arrange(desc(Cohort))
    
    return(
      list(
        control_sample_types = control_sample_types,
        control_assay_types = control_assay_types,
        control_samples = control_samples,
        olink_meta = olink_meta,
        id_mapping = id_mapping,
        mapped_data = mapped_data,
        lod_data = lod_data,
        exclude_proteins = exclude_proteins,
        n_samples_run = n_samples_run,
        n_samples = n_samples,
        n_proteins = n_proteins,
        n_plates = n_plates,
        manifest_filtered = manifest_filtered,
        cohort_number = cohort_number,
        disease_order = disease_order
      ))
  }


calculate_n_samples_assays <-
  function(data) {
    
    n_samples <-
      data |>
      distinct(DAid) |>
      nrow()
    
    n_proteins <-
      data |>
      distinct(Assay) |>
      nrow()
    
    return(list(n_samples = n_samples,
                n_proteins = n_proteins))
  }

calculate_sample_warnings <- 
  function(data) {
    data |> 
      filter(SampleQC != "PASS") |> 
      count(DAid) |> 
      mutate(perc_warning = n / qc$n_proteins) |> 
      filter(perc_warning > 0.5) 
  }

calculate_assay_warnings <- 
  function(data) {
    data |> 
      filter(!DAid %in% qc$control_samples) |> 
      group_by(Assay) |> 
      mutate(n_total = n_distinct(DAid)) |> 
      filter(AssayQC == "WARN") |> 
      mutate(n = n_distinct(DAid)) |> 
      ungroup() |> 
      distinct(Assay, n, n_total) |> 
      mutate(perc_warning = n / n_total) |> arrange(-n) |> 
      filter(perc_warning > 0.5)
  }

calculate_sample_detectability <- 
  function(data) {
    data |> 
      filter(above_LOD == "Yes") |> 
      group_by(DAid) |> 
      count()
  }

plot_cohort_n <- 
  function(plot_data,
           save_plot = NULL) {
    
    plot <- 
      plot_data |> 
      count(Cohort) |> 
      ggplot(aes(Cohort, n, fill = Cohort)) +
      geom_col() +
      scale_fill_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3")) +  
      geom_text(aes(label = n), vjust = -0.5) +
      theme_hpa(angled = T) +
      ylab("Number of samples")
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")), plot, height = 6, width = 6)
    }
    
    return(plot)
    
  }

plot_disease_n <- 
  function(plot_data,
           save_plot = NULL) {
    
    plot <- 
      plot_data |> 
      count(Cohort, Disease) |> 
      ggplot(aes(Disease, n, fill = Cohort)) +
      geom_col() +
      scale_fill_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3")) +  
      geom_text(aes(label = n), vjust = -0.5) +
      theme_hpa(angled = T) +
      ylab("Number of samples")
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")), plot, height = 6, width = 6)
    }
    
    return(plot)
    
  }

plot_plate_distribution <- 
  function(plot_data,
           save_plot = NULL)  {
    
    plot <- 
      plot_data |> 
      distinct(DAid, SampleID, Plate, SampleType) |> 
      mutate(Plate = gsub("P", "", Plate),
             Plate = as.numeric(Plate),
             Plate = as.factor(Plate)) |>
      filter(!is.na(Plate)) |> 
      left_join(qc$manifest_filtered, by = "DAid") |> 
      mutate(Cohort = ifelse(SampleType %in% qc$control_sample_types, "CONTROLS", Cohort))  |> 
      group_by(Plate, Cohort) |> 
      count() |> 
      ggplot(aes(Plate, n, fill = Cohort)) +
      geom_col() +
      scale_fill_manual(values = brewer.pal(n = 11, name = "Set3")) +  
      theme_hpa(angled = T) +
      ylab("Number of samples")
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")), plot, height = 6, width = 6)
    }
    
    return(plot)
  }

plot_NPX_missingness <- 
  function(plot_data, 
           plot_title,
           samples = NULL,
           save_plot = NULL) {
    
    
    if (!is.null(samples)) {
      plot_data_filt <- 
        plot_data |> 
        separate(SampleID, into = c("DAid", "Batch", "Plate"), sep = "-", remove = F) |> 
        filter(DAid %in% samples)
      
      plot <- 
        plot_data_filt |> 
        mutate(NPX_NA = ifelse(is.na(NPX), "Yes", "No"),
               Block = paste("Block", Block)) |> 
        count(DAid, PlateID, NPX_NA, Block) |>
        ggplot(aes(PlateID, n, fill = NPX_NA)) +
        geom_col() +
        facet_grid(Block~DAid, scales = "free_y") +
        scale_fill_manual(values = c("Yes" = "darkred", "No" = "grey")) +
        theme_hpa(angled = T) +
        ggtitle(plot_title)
      
    } else {
      plot <- 
        plot_data |> 
        mutate(NPX_NA = ifelse(is.na(NPX), "Yes", "No"),
               Block = paste("Block", Block)) |> 
        count(PlateID, NPX_NA, Block) |>
        ggplot(aes(PlateID, n, fill = NPX_NA)) +
        geom_col() +
        facet_wrap(~Block, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = c("Yes" = "darkred", "No" = "grey")) +
        theme_hpa(angled = T) +
        ggtitle(plot_title)
      
    }
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")),
             plot,
             height = 8,
             width = 8)
    }
    
    return(plot)
    
  }



plot_qc_warnings <-
  function(plot_data,
           plot_type,
           save_plot = NULL) {
    if (grepl("barplot", plot_type)) {
      qc_data_fail <-
        plot_data |>
        group_by(DAid) |>
        count(SampleQC) |>
        filter(!SampleQC %in% c("NA", "PASS"), !is.na(DAid)) |>
        arrange(-n) |>
        left_join(qc$manifest_filtered, by = "DAid")
      
      base_plot <-
        qc_data_fail |>
        ggplot(aes(fct_reorder(DAid, -n), n)) +
        geom_col() +
        geom_hline(
          yintercept = qc$n_proteins / 2,
          linetype = "dashed",
          color = "grey"
        ) +
        theme_hpa(angled = T) +
        xlab("Sample") +
        theme(axis.text.x = element_text(size = 5),
              axis.ticks.x = element_blank()) +
        ggtitle("Samples with QC warnings")
      
      if (plot_type == "barplot_cohort") {
        plot <-
          base_plot +
          geom_col(aes(fill = Cohort)) +
          geom_hline(
            yintercept = qc$n_proteins / 2,
            linetype = "dashed",
            color = "grey"
          ) +
          scale_fill_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3"))
        
        
      } else if (plot_type == "barplot_type") {
        plot <-
          base_plot +
          geom_col(aes(fill = SampleQC)) +
          geom_hline(
            yintercept = qc$n_proteins / 2,
            linetype = "dashed",
            color = "grey"
          )
        
      }
      
      if (!is.null(save_plot)) {
        ggsave(savepath(paste0(save_plot, ".png")),
               plot,
               height = 4,
               width = 10)
      }
      
    }
    
    else if (plot_type == "beeswarm") {
      qc_data_pass <-
        plot_data |>
        group_by(DAid) |>
        filter(SampleQC == "PASS") |>
        count(SampleQC) |>
        left_join(qc$manifest_filtered, by = "DAid") |>
        filter(Cohort != "CTRL")
      
      outliers <- 
        qc_data_pass |> 
        filter(n <  qc$n_proteins / 2)
      
      plot <- 
        qc_data_pass |>
        ggplot(aes(Cohort, n, color = Cohort)) +
        geom_quasirandom() +
        geom_hline(
          yintercept = qc$n_proteins / 2,
          color = "grey40",
          linetype = "dashed"
        ) +
        geom_text_repel(aes(label = DAid), data = outliers, show.legend = F) +
        scale_color_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3")) +
        theme_hpa(angled = T) +
        ggtitle("Number of proteins that pass QC per sample")
      
      if (!is.null(save_plot)) {
        ggsave(savepath(paste0(save_plot, ".png")),
               plot,
               height = 8,
               width = 8)
      }
      
    } else {
      stop("Invalid plot type specified.")
    }
    
    return(plot)
  }


plot_assay_warnings <-
  function(plot_data,
           save_plot = NULL) {
    plot <- 
      plot_data |>
      group_by(Assay) |>
      count(AssayQC) |>
      filter(!AssayQC %in% c("NA", "PASS")) |>
      arrange(-n) |>
      ggplot(aes(fct_reorder(Assay, -n), n)) +
      geom_col() +
      geom_text(aes(label = n), vjust = -0.5) +
      theme_hpa(angle = T) +
      xlab("Assay")
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")),
             plot,
             height = 5,
             width = 10)
    }
    
    return(plot)
  }

plot_sample_detectability <-
  function(plot_data,
           save_plot = NULL) {
    
    plot <- 
      plot_data |>
      left_join(qc$manifest_filtered, by = "DAid") |>
      mutate(
        Disease = ifelse(Disease == "Healthy", paste(Disease, Cohort, sep = " - "), Disease),
        Disease = factor(Disease, levels = qc$disease_order$Disease)
      ) |>
      ggplot(aes(Disease, n, fill = Cohort, color = Cohort)) +
      geom_quasirandom() +
      geom_boxplot(outlier.color = NA,
                   alpha = 0.5,
                   color = "black") +
      scale_fill_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3")) +
      scale_color_manual(values = brewer.pal(n = qc$cohort_number, name = "Set3")) +
      theme_hpa(angled = T)
    
    if (!is.null(save_plot)) {
      ggsave(savepath(paste0(save_plot, ".png")),
             plot,
             height = 6,
             width = 8)
    }
    
    return(plot)
  }

plot_patient_heatmap <-
  function(plot_data,
           select_id) {
    data <-
      plot_data |>
      left_join(qc$id_mapping, by = "SampleID") |>
      filter(DAid %in% select_id)
    
    if (length(select_id) == 1) {
      plot <-
        data |>
        select(Plate, OlinkID, NPX) |>
        pivot_wider(names_from = OlinkID, values_from = NPX) |>
        column_to_rownames("Plate") |>
        t() |>
        cor(method = "spearman", use = "complete.obs") |>
        pheatmap()
      
    } else {
      ann <-
        qc$id_mapping |>
        distinct(DAid, Plate) |>
        mutate(ID = paste(DAid, Plate, sep = "_")) |>
        column_to_rownames("ID")
      
      data |>
        mutate(ID = paste(DAid, Plate, sep = "_")) |>
        select(ID, OlinkID, NPX) |>
        pivot_wider(names_from = OlinkID, values_from = NPX) |>
        column_to_rownames("ID") |>
        t() |>
        cor(method = "spearman", use = "complete.obs") |>
        pheatmap(
          annotation_row  = ann,
          show_colnames = F,
          show_rownames = F
        )
    }
    
    
  }

plot_correlogram <-
  function(data,
           title,
           n_plates) {
    data |>
      select(Plate, OlinkID, NPX) |>
      pivot_wider(names_from = Plate, values_from = NPX) |>
      ggpairs(c(2:n_plates + 1),
              lower = list(continuous = "points", size = 0.05)) +
      ggtitle(title)
  }

plot_dim_reduction <- function(plot_data,
                               plot_type,
                               plot_title,
                               save_plot = NULL,
                               method = c("umap", "pca")) {
  
  method <- match.arg(method)
  
  # Determine which axis names to use
  axis1 <- sym(ifelse(method == "umap", "UMAP1", "PC1"))
  axis2 <- sym(ifelse(method == "umap", "UMAP2", "PC2"))
  
  plot_data <- plot_data |>
    rename(DAid = Sample) |>
    left_join(qc$manifest_filtered, by = "DAid")
  
  if (plot_type == "class") {
    plot <- plot_data |>
      ggplot(aes(x = !!axis1, y = !!axis2, color = Class)) +
      geom_point(alpha = 0.7, size = 0.8) +
      scale_color_manual(values = pal_class_2) +
      theme_hpa() +
      ggtitle(plot_title)
    
  } else if (plot_type == "disease") {
    plot <- plot_data |>
      ggplot(aes(x = !!axis1, y = !!axis2, color = Disease)) +
      geom_point(alpha = 0.7, size = 0.8) +
      theme_hpa() +
      ggtitle(plot_title)
    
  } else if (plot_type == "plate") {
    plot <- plot_data |>
      left_join(qc$id_mapping, by = "DAid") |>
      ggplot(aes(x = !!axis1, y = !!axis2, color = Plate)) +
      geom_point(alpha = 0.7, size = 0.8) +
      theme_hpa() +
      ggtitle(plot_title)
    
  } else if (plot_type == "detectability") {
    plot <- plot_data |>
      left_join(sample_detected_proteins, by = "DAid") |>
      ggplot(aes(x = !!axis1, y = !!axis2, color = n)) +
      geom_point(alpha = 0.7, size = 0.8) +
      scale_color_viridis_c() +
      theme_hpa() +
      ggtitle(plot_title)
    
  } else {
    stop("Invalid plot type specified. Choose from 'class', 'disease', 'plate', or 'detectability'.")
  }
  
  if (!is.null(save_plot)) {
    ggsave(savepath(paste0(save_plot, ".png")),
           plot,
           height = 8,
           width = 8)
  }
  
  return(plot)
}
