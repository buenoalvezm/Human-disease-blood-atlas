# Data Preprocessing – Human Disease Blood Atlas

This internal repository is used to track the **preprocessing and data preparation steps for the Olink and Somalogic proteomics data** included in the **Human Disease Blood Atlas**.

The repository includes:
- 🧾 Preprocessing scripts for Olink/Somalogic data
- 🧪 Sample and assay-level quality control
- 🔍 Exploratory analyses and investigations of assay characteristics
- 🗂️ Metadata preparation and harmonization
- 📂 Output files and intermediate results from processing steps

This internal repository documents the preprocessing and preparation of Olink proteomics data included in the Human Disease Blood Atlas.

## 📂 Folder Structure (scripts/)

The main working folders include:

```bash
functions/              # Custom R functions used across multiple scripts
assay_exploration/      # Analysis of assay characteristics (e.g., detection rates, dynamic range, plate effects)
sample_overview/        # Summaries and visualizations of sample metadata and inclusion criteria
data_preprocessing_phase2_olink/ # Scripts for normalization, filtering, and final Olink dataset preparation (phase 2)
data_preprocessing_phase2_somalogic/ # Scripts for normalization, filtering, and final Somalogic dataset preparation (phase 2)
data_exploration/       # Global data exploration including PCA/UMAP, distributions, and missingness
disease_explorations/   # Exploratory analyses for some disease groups
hpa_data_generation/    # Code related to the creation of HPA-ready datasets
internal_qc_scripts/    # Internal QC routines for sample, plate, and batch evaluation
```

## 📝 Notes
- This repository is for internal use only
- Data is not public
- The repository serves as both workflow tracking and documentation hub for the project
