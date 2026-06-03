# Exploratory Analysis of L1000 Gene Expression Signatures for Perturbation Mapping

**Medical Bioinformatics Internship Project**  
Department of Medical Bioinformatics, University Medical Center Göttingen  
**Author:** Shayan Shokrani  
**Supervisor:** Fazel Amirvahedi  
**PI:** Prof. Dr. Michael Altenbuchinger  

---

## Overview

Large-scale perturbation datasets like LINCS L1000 profile how small molecules reshape gene expression across diverse cell lines, offering a rich resource for drug discovery and repurposing. However, extracting true drug-induced signals is challenging: cell-line identity dominates expression variance, batch effects between studies introduce technical noise, and compound annotations are heterogeneous across databases.

This project builds a clean, annotated, batch-corrected L1000 resource from two GEO releases and evaluates machine learning models for perturbation prediction. The central question is whether a **delta (residual) modeling strategy** — explicitly subtracting cell-line baselines before training — can help models focus on the smaller, drug-specific signal rather than learning cell identity.

---

## Repository Structure

```
.
├── notebooks/
│   ├── <01_preprocessing notebook>        # Phase 1: data acquisition and preprocessing
│   ├── <02_moa_annotation notebook>       # Phase 2: ChEMBL MoA annotation
│   ├── <03_eda_batch_correction notebook> # Phase 3: PCA/UMAP, splitting, ComBat
│   └── <04_modeling notebook>             # Phase 4: baseline and delta modeling
│
├── figures/                               # PCA, UMAP, and other exploratory plots
│
├── pearson_Best_RF_comparison.jpg         # Pearson scores – best Random Forest folds
├── pearson_Worst_RF_comparison.jpg        # Pearson scores – worst Random Forest folds
├── r2_Best_RF_comparison.jpg              # R² scores – best Random Forest folds
└── r2_Worst_RF_comparison.jpg            # R² scores – worst Random Forest folds
```

> **Note:** Notebook filenames in angle brackets are placeholders — replace them with the actual filenames from your `notebooks/` folder.

---

## Data Access

### Raw Input Data (L1000 Level 3 from GEO)

This project uses **Level 3 (Q2NORM)** L1000 gene expression data from two GEO series. Level 3 data represents quantile-normalized expression values for the 978 directly measured landmark genes. The LINCS User Guide describing all data files and metadata columns is available here:

> [GEO CMap LINCS User Guide v2.1](https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit?usp=sharing)

| Series | Name | Samples |
|--------|------|---------|
| [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) | LINCS L1000 Phase II | 345,976 |
| [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) | LINCS L1000 Phase I | 1,319,138 |

Download the `level3_Q2NORM` matrix files and the `sig_info`, `pert_info`, and `gene_info` metadata tables from each GEO accession page.

### Preprocessed & Combined Dataset

The final combined, harmonized, and annotated dataset (1,038,163 samples × 978 landmark genes) is available as an AnnData `.h5ad` file — ready for direct use with [Scanpy](https://scanpy.readthedocs.io/):

> [Download combined L1000 dataset (Google Drive)](https://drive.google.com/file/d/16_T0AaxKJsjhZiXManpSVhR6O5sxqhV6/view?usp=drivesdk)

This file contains the expression matrix, all sample metadata (cell line, compound, dose, treatment time, MoA, action type), and gene annotations in a single file. If you are only interested in the modeling or EDA phases, start here rather than re-running Phase 1.

### External Annotation: ChEMBL

Mechanism-of-action (MoA) and action-type annotations were integrated from [ChEMBL](https://www.ebi.ac.uk/chembl/) via InChIKey matching. These annotations are already embedded in the `.h5ad` metadata above, so you do not need to re-download ChEMBL unless you want to reproduce Phase 2 from scratch.

---

## Pipeline

The analysis is organized into four sequential phases, each corresponding to one notebook.

### Phase 1 — Data Acquisition and Preprocessing

The two GEO datasets are merged after restricting features to the 978 landmark genes and filtering to small-molecule treatment profiles (`trt_cp`) and DMSO vehicle controls (`ctl_vehicle`). Dose strings are harmonized to micromolar units, compounds with fewer than five samples are removed, and SMILES validity is checked with RDKit. The result is stored as an AnnData object in H5AD format, which consolidates the expression matrix, sample metadata, and gene annotations into a single file. The combined dataset comprises 1,038,163 samples from 82 cell lines and 17,936 unique compounds.

### Phase 2 — MoA Annotation via ChEMBL

Each compound is annotated with its mechanism of action and action type by joining the L1000 compound table to ChEMBL using InChIKey as a standardized molecular identifier — more reliable than name-based matching due to synonym inconsistencies. Where a compound maps to multiple MoAs in ChEMBL, the dominant MoA/action-type combination by frequency is retained. After annotation, 750,371 samples carry an assigned MoA spanning 35 distinct action types.

### Phase 3 — Exploratory Analysis, Splitting Strategy, and ComBat Correction

PCA and UMAP visualizations reveal that cell-line identity is the dominant source of variance, with additional systematic structure introduced by dataset phase (GSE70138 vs. GSE92742) and treatment time (especially 24 h incubations). To correct for these technical effects, ComBat batch-effect correction is applied — but crucially, **only within training folds after splitting**, to prevent leakage of test-set information into the correction.

The splitting strategy is biologically informed: a stratification key combining MoA, cell line, and discretized dose identifies "rare" contexts (appearing only once), which are confined entirely to training data to ensure the test set reflects only well-represented biological conditions.

### Phase 4 — Baseline and Delta Modeling

Four models are evaluated in two settings each — raw expression and delta (residual) expression. The **delta formulation** subtracts the cell-specific DMSO mean before training, forcing models to predict drug-induced deviations rather than cell identity. A simple cell-specific mean baseline already achieves Pearson ≈ 0.877 (R² ≈ 0.32), showing how much variance is explained by cell identity alone. Ridge regression shows no improvement under delta modeling, while Random Forest benefits substantially (Pearson ≈ 0.885, R² ≈ 0.37), suggesting the residual drug signal has non-linear structure that linear models cannot exploit.

---

## Key Results

| Model | Setting | Pearson | R² |
|-------|---------|---------|-----|
| Cell-specific mean (baseline) | — | 0.877 | 0.320 |
| Ridge Regression | Raw expression | ~0.877 | ~0.320 |
| Ridge Regression | Delta (residual) | ~0.877 | ~0.320 |
| Random Forest | Raw expression | < 0.877 | < 0.320 |
| **Random Forest** | **Delta (residual)** | **0.885** | **0.370** |

Result comparison plots per fold are available in the root of the repository (`pearson_*` and `r2_*` figures).

---

## Requirements

This project was developed with Python 3.10+. Install all dependencies with:

```bash
pip install scanpy anndata pandas numpy scikit-learn matplotlib seaborn rdkit
```

---

## References

Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell.* 2017;171(6):1437-1452. doi:10.1016/j.cell.2017.10.049

Davies M, et al. ChEMBL: an open source database for small molecule bioactivity data. *Nucleic Acids Research.* 2025;53(D1):D218-D225. doi:10.1093/nar/gkae1085

Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics.* 2007;8(1):118-127. doi:10.1093/biostatistics/kxj037

Goodman JM, et al. InChI version 1.06: now more than 99.99% reliable. *J Cheminform.* 2021;13:40. doi:10.1186/s13321-021-00517-z
