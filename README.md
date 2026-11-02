# Cardiogenic Shock: Multimodal Analysis of SCAI Stage D vs C

Analysis code for the cardiogenic shock (CS) manuscript (Bin Chen Lab, Michigan State University).
The study compares SCAI shock stage **D** with stage **C** using three data modalities:

- **Clinical and laboratory data:** admission vs. shock-onset labs, demographics (Table 1)
- **Bulk RNA-seq:** differential expression, pathway enrichment, ssGSEA, WGCNA, CIBERSORTx immune deconvolution
- **Cytokine panel:** univariate tests, ROC models, sex and age interaction analyses

These are then integrated into cross-modal prediction and association models (elastic net, Firth logistic regression, LOOCV).

> **Data are not included in this repository.** The clinical, cytokine, and sample metadata contain patient-level information and are not publicly available. See [Data availability](#data-availability).

---

## Repository structure

```
scripts/
├── analysis/
│   ├── 1_lab_data.R             # Clinical/lab cleaning, admission vs. shock, Table 1
│   ├── 2_bulk_data.R            # DESeq2, ORA, ssGSEA, WGCNA, GEO export
│   ├── 3_cytokine_data.R        # Cytokine C vs D tests, ROC, confounder models
│   ├── 4_prediction.R           # Integrated clinical + cytokine + RNA prediction
│   ├── 5_association_model.R    # Elastic-net association model, feature correlations
│   └── 6_cross_modal_data.R     # Multimodal PCA, correlations, ORs, LOOCV by modality
└── figures/
    ├── figure_2.R               # Clinical / lab figures (+ supplementary)
    ├── figure_3.R               # Transcriptomic figures (+ supplementary)
    ├── figure_4.R               # Cytokine figures
    └── figure_5.R               # Multimodal integration figures
```

Figure 1 (study overview) is not generated in R.

---

## Setup

All scripts set the working directory to the folder that contains `cs_project/`:

```r
setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
```

To run them elsewhere, set the environment variable before starting R:

```bash
export CS_ROOT=/path/to/folder/containing/cs_project
```

---

## Input data

All paths are relative to `CS_ROOT`.

| File | Description | Patient-level | Used by |
|---|---|---|---|
| `cs_project/clinical_data/CS_full_clinical_data_set.csv` | Raw clinical and laboratory data | Yes | 1 |
| `cs_project/bulk_data/count_data.RData` | Bulk RNA-seq matrices: `log2.read.count.matrix`, `log2.tpm.matrix`, `log2.fpkm.matrix` | No (expression only) | 2, 4, 5, 6, Fig 3 |
| `cs_project/bulk_data/metadata_bulk.csv` | Sample metadata linking specimens to SCAI stage (`SR.specimen.ID`, `SCAI.Stage`) | Yes | 2, 4, 5, 6, Fig 3 |
| `cs_project/cytokine_data/CScytokineData.xlsx` | Cytokine panel measurements per specimen | Yes | 3, 4, 5, 6, Fig 4 |
| `cs_project/gene_info_biomart.csv` | Ensembl gene ID to gene name mapping (`gene_id`, `gene_name`) from BioMart | No | 2, 4, 5, 6, Fig 3 |
| `cs_project/bulk_data/pathway sets/*.gmt` | Custom gene sets for ssGSEA | No | 2 |
| `cs_project/bulk_data/CIBERSORTx_Job3_Results.csv` | Immune cell fractions from the CIBERSORTx web tool (run outside R) | No | Fig 3 |

Samples across modalities are linked by the normalized SR specimen ID.

---

## Run order

```
1_lab_data.R ─────────────┐
2_bulk_data.R ────────────┤
3_cytokine_data.R ────────┼──> 4_prediction.R
                          ├──> 5_association_model.R
                          └──> 6_cross_modal_data.R

1 ──> figure_2.R    2 + CIBERSORTx ──> figure_3.R
3 ──> figure_4.R    6 ──> figure_5.R
```

Outputs are written to `results/` subfolders inside `cs_project/` (e.g. `bulk_data/results/`, `multi_modality/results/`) and to `cs_project/figures/FigureN_outputs/`.

---

## R packages

**CRAN:** dplyr, tidyr, tibble, purrr, stringr, forcats, readxl, openxlsx, writexl, tidyverse, ggplot2, ggrepel, ggtext, ggpubr, ggh4x, ggnewscale, ggplotify, patchwork, pheatmap, RColorBrewer, tidytext, broom, caret, glmnet, logistf, pROC, FactoMineR, factoextra, WGCNA, msigdbr

**Bioconductor:** DESeq2, apeglm, clusterProfiler, enrichplot, ReactomePA, reactome.db, org.Hs.eg.db, AnnotationDbi, GSVA, KEGGREST

Package versions: see `sessionInfo.txt` (to be added).

---

## Data availability

- **Bulk RNA-seq:** deposited in GEO ([GSE349137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349137)). `2_bulk_data.R` contains the GEO export.
- **Clinical, cytokine, and sample metadata:** contain patient-level information and are not publicly available. De-identified data are available from the corresponding author on reasonable request, subject to IRB approval and a data use agreement.

---

## Contact

Shubham Koirala ([ORCID 0000-0002-4399-2041](https://orcid.org/0000-0002-4399-2041))
Bin Chen Lab
