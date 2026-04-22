# NeuroDE-LMM  
### Differential Gene Expression Analysis with Linear Mixed Models for Neuropsychiatric Transcriptomic Data

**NeuroDE-LMM** is a reproducible bioinformatics pipeline for **differential gene expression (DGE)** analysis in human brain transcriptomic datasets.  
The project focuses on identifying genes associated with neuropsychiatric conditions by combining **data preprocessing**, **normalization**, and **per-gene Linear Mixed Models (LMMs)** with careful control of biological and technical confounders.

This repository was built as a personal project to reproduce and extend transcriptomic analyses used in neuropsychiatric research, with particular emphasis on:
- **confounder-aware statistical modeling**
- **reproducible data analysis**
- **clear result visualization**
- **downstream integration with co-expression networks and pathway enrichment**

## Why this project?
Gene expression datasets are high-dimensional, noisy, and strongly affected by confounding factors such as **batch effects**, **brain bank differences**, **age**, **sex**, and other sample-specific variables.  
This project addresses that problem by using **Linear Mixed Models** to estimate diagnosis-related effects at the gene level while accounting for nuisance variability.

The result is a pipeline that produces:
- ranked lists of differentially expressed genes
- multiple-testing corrected significance values
- effect sizes and confidence intervals
- publication-style volcano plots
- outputs ready for downstream analyses such as **WGCNA**, **GO enrichment**, and **imaging-transcriptomics**

## Main components
- **Preprocessing and normalization (R)**  
  Cleaning metadata, harmonizing samples, and preparing normalized expression matrices.

- **Differential expression with LMM (Python)**  
  Per-gene mixed-effects modeling using `statsmodels` to estimate diagnosis effects while controlling for relevant covariates.

- **Visualization and export**  
  Automated generation of volcano plots and export of structured results for downstream analyses.

## Repository structure
- `code/` → analysis scripts  
- `Exports/` → processed data tables  
- `LMMResults/` → differential expression outputs  
- `plots/` → generated visualizations  

## Technical stack
**Python:** `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`  
**R:** preprocessing and normalization scripts for transcriptomic data  

## Project goal
The long-term goal of this repository is to provide a transparent and extensible framework for **transcriptomic differential expression analysis in neuropsychiatric research**, and to connect those results with **network biology** and **brain imaging approaches**.
