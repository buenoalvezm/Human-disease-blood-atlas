# Olink Data Preprocessing – Human Disease Blood Atlas

This internal repository is used to track the **preprocessing and data preparation steps for the Olink proteomics data** included in the **Human Disease Blood Atlas**.

The repository includes:
- Preprocessing scripts for NPX data
- Sample and assay-level quality control
- Exploratory analyses and investigations of assay characteristics
- Metadata preparation and harmonization
- Output files and intermediate results from processing steps

HTML reports generated from R Markdown will also be included to support reproducibility and documentation of key steps.

## Folder Structure (scripts/)

The main working folders include:

- `functions/` – Custom R functions used across multiple scripts
- `assay_exploration/` – Analysis of assay characteristics (e.g., detection rates, dynamic range, plate effects)
- `sample_overview/` – Summaries and visualizations of sample metadata and inclusion criteria
- `data_preprocessing_phase2/` – Scripts for normalization, filtering, and final dataset preparation (phase 2)
- `data_exploration/` – Global data exploration including PCA/UMAP, distributions, and missingness
- `disease_explorations/` – Exploratory analyses for some disease groups
- `hpa_data_generation/` – Code related to the creation of HPA-ready datasets
- `internal_qc_scripts/` – Internal QC routines for sample, plate, and batch evaluation

---

This repository is intended for internal use only.
