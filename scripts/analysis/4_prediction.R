setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Cross-Modality Prediction of SCAI Stage D vs C
# Clinical + Cytokine + Transcriptomic Gene-Set Scores
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(glmnet)
library(pROC)
library(caret)
library(readxl)
library(ggplot2)
library(pheatmap)

set.seed(123)

output_dir <- "cs_project/integrated_prediction/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD SAMPLE METADATA
# ============================================================

metadata <- read.csv("cs_project/bulk_data/metadata_bulk.csv") %>%
  mutate(
    specimen_id = str_trim(SR.specimen.ID),
    scai_stage = toupper(str_trim(SCAI.Stage))
  ) %>%
  filter(scai_stage %in% c("C", "D")) %>%
  mutate(
    outcome = ifelse(scai_stage == "D", 1, 0)
  ) %>%
  dplyr::select(specimen_id, scai_stage, outcome)

# ============================================================
# 2. LOAD AND PREPARE CLINICAL FEATURES
# ============================================================

clinical <- read.csv("cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv") %>%
  mutate(
    specimen_id = str_replace(specimen_id, "^SH-", "SR-"),
    scai_stage = toupper(str_trim(scai_stage))
  ) %>%
  filter(scai_stage %in% c("C", "D"))

clinical_features <- clinical %>%
  dplyr::select(
    specimen_id,
    
    # admission/metabolic/end-organ dysfunction
    matches("admission_labs_lactate"),
    matches("admission_labs_bicarbonate"),
    matches("admission_labs_creatinine"),
    matches("admission_labs_mdrd_egfr"),
    matches("admission_labs_bilirubin"),
    matches("admission_labs_albumin"),
    matches("admission_labs_glucose"),
    matches("admission_labs_nt_pro_bnp"),
    matches("admission_labs_troponin"),
    matches("admission_labs_c_reactive_protein"),
    
    # shock onset labs
    matches("shock_onset_labs_lactate"),
    matches("shock_onset_labs_bicarbonate"),
    matches("shock_onset_labs_creatinine"),
    matches("shock_onset_labs_mdrd_egfr"),
    matches("shock_onset_labs_bilirubin"),
    matches("shock_onset_labs_albumin"),
    matches("shock_onset_labs_glucose"),
    matches("shock_onset_labs_nt_pro_bnp"),
    matches("shock_onset_labs_troponin"),
    matches("shock_onset_labs_c_reactive_protein"),
    
    # hemodynamics
    matches("shock_onset_rhc_fick_ci"),
    matches("shock_onset_rhc_fick_co"),
    matches("shock_onset_rhc_svr"),
    matches("shock_onset_rhc_pvr"),
    matches("shock_onset_rhc_pcwp"),
    matches("shock_onset_pac_fick_cardiac_index"),
    matches("shock_onset_pac_fick_cardiac_output"),
    matches("shock_onset_pac_systemic_vascular_resistance"),
    matches("shock_onset_pac_pulmonary_artery_pulsatility_index"),
    
    # vasoactive support
    matches("norepinephrine_max"),
    matches("epinephrine_max"),
    matches("dobutamine_max"),
    matches("milrinone_max"),
    matches("vasopressin_max")
  )

clinical_features <- clinical_features %>%
  mutate(across(-specimen_id, ~ suppressWarnings(as.numeric(.))))

# Remove variables with too much missingness
clinical_features <- clinical_features %>%
  dplyr::select(
    specimen_id,
    where(~ mean(is.na(.)) < 0.50)
  )
clinical_features <- clinical_features %>%
  mutate(
    specimen_id = specimen_id %>%
      str_replace("^SR-24-0", "SR24") %>%
      str_replace("^SR-25-0", "SR25") %>%
      str_replace("-", "")
  )

# ============================================================
# 3. LOAD AND PREPARE CYTOKINE FEATURES
# ============================================================
cytokine_raw <- read_excel("cs_project/cytokine_data/CScytokineData.xlsx") %>%
  rename_with(~ str_replace_all(., "-", "")) %>%
  rename_with(~ str_remove_all(., " \\(pg/ml\\)")) %>%
  rename_with(~ str_replace_all(., " ", "")) %>%
  mutate(
    SCAIStage = toupper(SCAIStage),
    specimen_id = str_trim(SRspecimenID)
  ) %>%
  filter(SCAIStage %in% c("C", "D"))

meta_cols <- c("StudyID", "Age", "Sex", "Race", "SCAIStage", "SRspecimenID", "specimen_id")

cytokine_features <- cytokine_raw %>%
  dplyr::select(
    specimen_id,
    where(is.numeric)
  )

cytokine_features <- cytokine_features %>%
  mutate(across(-specimen_id, ~ log1p(.)))

cytokine_features <- cytokine_features %>%
  dplyr::select(
    specimen_id,
    where(~ mean(is.na(.)) < 0.50),
    where(~ length(unique(na.omit(.))) > 3)
  )

# Optional manually prioritized cytokines
priority_cytokines <- c(
  "ST2", "IL6", "IL8", "S100A8", "VEGFA",
  "TroponinI", "TNF", "IL1B", "MCP1"
)

cytokine_features_priority <- cytokine_features %>%
  dplyr::select(
    specimen_id,
    any_of(priority_cytokines)
  )
cytokine_features_priority <- cytokine_features_priority %>%
  mutate(
    specimen_id = specimen_id %>%
      str_replace("^SR-24-0", "SR24") %>%
      str_replace("^SR-25-0", "SR25") %>%
      str_replace("-", "")
  )
# ============================================================
# 4. LOAD RNA EXPRESSION AND BUILD GENE-SET SCORES
# ============================================================

load("cs_project/bulk_data/count_data.RData")

expr <- log2.tpm.matrix

mapping <- read.csv("cs_project/gene_info_biomart.csv") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)

expr <- expr[, intersect(colnames(expr), metadata$specimen_id)]

expr_df <- as.data.frame(expr) %>%
  rownames_to_column("gene_id") %>%
  left_join(mapping, by = "gene_id") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

expr_mat <- expr_df %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# test <- (expr_df)
# dim(test)
# write.xlsx(test, file.path(output_dir, "rna_data_nabas.xlsx"))

# Z-score genes across samples
expr_z <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0

gene_sets <- list(
  metabolic_stress_hyperglycemia = c(
    "SLC2A1", "SLC2A3", "HK2", "PFKFB3", "LDHA",
    "PDK1", "HIF1A", "DDIT3", "XBP1", "ATF4",
    "FOS", "JUN", "NR4A1", "TXNIP", "G6PD",
    "ENO1", "ALDOA", "PKM"
  ),
  
  inflammation_capillary_leak = c(
    "IL6", "IL1B", "TNF", "NFKBIA", "NFKB1",
    "RELA", "CXCL8", "CCL2", "CCL3", "CCL4",
    "SOCS3", "S100A8", "S100A9", "PTGS2",
    "MMP9", "LCN2", "OSM", "IL1RN", "STAT3"
  ),
  
  endothelial_activation_permeability = c(
    "VWF", "ICAM1", "VCAM1", "SELE", "SELP",
    "PECAM1", "CDH5", "CLDN5", "KDR", "FLT1",
    "TEK", "ANGPT1", "ANGPT2", "VEGFA",
    "SERPINE1", "THBD", "ENG", "ESAM"
  ),
  
  vasoplegia_vascular_tone = c(
    "NOS2", "NOS3", "ARG1", "ARG2", "GCH1",
    "EDN1", "ADM", "CALCRL", "RAMP2", "RAMP3",
    "PTGS2", "ADRA1A", "ADRA1B", "ADRA1D",
    "ADRB1", "ADRB2", "AVPR1A", "AVPR1B",
    "AGTR1", "AGTR2"
  ),
  
  neutrophil_activation = c(
    "S100A8", "S100A9", "MPO", "ELANE", "MMP8",
    "MMP9", "CEACAM8", "FCGR3B", "CXCR1",
    "CXCR2", "OLFM4", "LCN2", "DEFA1",
    "DEFA3", "AZU1", "CTSG"
  ),
  
  monocyte_macrophage_activation = c(
    "CD14", "FCN1", "LST1", "VCAN", "CCR2",
    "CCL2", "CD68", "MSR1", "ITGAM", "TYROBP",
    "FCGR1A", "FCGR2A", "LYZ", "S100A12", "TREM1"
  ),
  
  hypoxia_mitochondrial_stress = c(
    "HIF1A", "EGLN1", "EGLN3", "VEGFA", "BNIP3",
    "NDUFA1", "NDUFB8", "COX4I1", "ATP5F1A",
    "SOD2", "TXNIP", "UCP2", "PGK1", "ENO1",
    "CA9", "ADM"
  ),
  
  coagulation_complement = c(
    "C3", "C5", "C5AR1", "CFD", "CFB", "F3",
    "SERPINE1", "THBD", "PLAU", "PLAUR",
    "ITGA2B", "PF4", "PPBP", "GP1BA",
    "GP9", "CLU", "C1QA", "C1QB", "C1QC"
  )
)

score_gene_set <- function(expr_z, genes) {
  genes_present <- intersect(genes, rownames(expr_z))
  
  if (length(genes_present) < 2) {
    return(rep(NA_real_, ncol(expr_z)))
  }
  
  colMeans(expr_z[genes_present, , drop = FALSE], na.rm = TRUE)
}

rna_scores <- map_dfc(
  gene_sets,
  ~ score_gene_set(expr_z, .x)
)

rna_scores <- as.data.frame(rna_scores) %>%
  mutate(specimen_id = colnames(expr_z)) %>%
  relocate(specimen_id)

# ============================================================
# 5. MERGE MODALITIES
# ============================================================
# clinical_features$specimen_id
# cytokine_features$specimen_id
# rna_scores$specimen_id
# head(integrated_data)

integrated_data <- metadata %>%
  left_join(clinical_features, by = "specimen_id") %>%
  left_join(cytokine_features_priority, by = "specimen_id") %>%
  left_join(rna_scores, by = "specimen_id")

write.csv(
  integrated_data,
  file.path(output_dir, "integrated_feature_matrix.csv"),
  row.names = FALSE
)

# ============================================================
# 6. MODELING FUNCTIONS
# ============================================================

impute_train_apply <- function(train_x, test_x) {
  medians <- apply(train_x, 2, function(z) median(z, na.rm = TRUE))
  
  for (j in seq_len(ncol(train_x))) {
    train_x[is.na(train_x[, j]), j] <- medians[j]
    test_x[is.na(test_x[, j]), j] <- medians[j]
  }
  
  list(train = train_x, test = test_x)
}

scale_train_apply <- function(train_x, test_x) {
  means <- apply(train_x, 2, mean)
  sds <- apply(train_x, 2, sd)
  sds[sds == 0 | is.na(sds)] <- 1
  
  train_scaled <- sweep(train_x, 2, means, "-")
  train_scaled <- sweep(train_scaled, 2, sds, "/")
  
  test_scaled <- sweep(test_x, 2, means, "-")
  test_scaled <- sweep(test_scaled, 2, sds, "/")
  
  list(train = train_scaled, test = test_scaled)
}

fit_loocv_elastic_net <- function(df, feature_cols, alpha = 0.5) {
  df <- df %>%
    dplyr::select(specimen_id, outcome, all_of(feature_cols))
  
  df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
  
  feature_cols <- setdiff(colnames(df), c("specimen_id", "outcome"))
  
  # remove near-zero variance
  nzv <- nearZeroVar(df[, feature_cols, drop = FALSE])
  
  if (length(nzv) > 0) {
    feature_cols <- feature_cols[-nzv]
  }
  
  if (length(feature_cols) < 2) {
    stop("Not enough usable features.")
  }
  
  predictions <- data.frame()
  
  selected_features <- list()
  
  for (i in seq_len(nrow(df))) {
    train <- df[-i, ]
    test <- df[i, ]
    
    x_train <- as.matrix(train[, feature_cols, drop = FALSE])
    x_test <- as.matrix(test[, feature_cols, drop = FALSE])
    y_train <- train$outcome
    
    imp <- impute_train_apply(x_train, x_test)
    x_train <- imp$train
    x_test <- imp$test
    
    sc <- scale_train_apply(x_train, x_test)
    x_train <- sc$train
    x_test <- sc$test
    
    # class weighting
    weights <- ifelse(
      y_train == 1,
      sum(y_train == 0) / sum(y_train == 1),
      1
    )
    
    cvfit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = alpha,
      weights = weights,
      type.measure = "auc",
      nfolds = min(5, length(unique(y_train)) * 3)
    )
    
    pred <- as.numeric(
      predict(
        cvfit,
        newx = x_test,
        s = "lambda.min",
        type = "response"
      )
    )
    
    coef_df <- as.matrix(coef(cvfit, s = "lambda.min"))
    selected <- rownames(coef_df)[coef_df[, 1] != 0]
    selected <- setdiff(selected, "(Intercept)")
    
    selected_features[[i]] <- selected
    
    predictions <- bind_rows(
      predictions,
      tibble(
        specimen_id = test$specimen_id,
        observed = test$outcome,
        predicted_probability = pred
      )
    )
  }
  
  roc_obj <- roc(
    predictions$observed,
    predictions$predicted_probability,
    quiet = TRUE
  )
  
  threshold <- 0.5
  
  predictions <- predictions %>%
    mutate(
      predicted_class = ifelse(predicted_probability >= threshold, 1, 0)
    )
  
  sensitivity <- sum(predictions$observed == 1 & predictions$predicted_class == 1) /
    sum(predictions$observed == 1)
  
  specificity <- sum(predictions$observed == 0 & predictions$predicted_class == 0) /
    sum(predictions$observed == 0)
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  feature_frequency <- table(unlist(selected_features)) %>%
    as.data.frame() %>%
    arrange(desc(Freq))
  
  colnames(feature_frequency) <- c("feature", "selection_frequency")
  
  list(
    predictions = predictions,
    auc = as.numeric(auc(roc_obj)),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_accuracy,
    feature_frequency = feature_frequency
  )
}

# ============================================================
# 7. DEFINE MODALITY-SPECIFIC FEATURE SETS
# ============================================================

clinical_cols <- setdiff(
  colnames(clinical_features),
  "specimen_id"
)

cytokine_cols <- setdiff(
  colnames(cytokine_features_priority),
  "specimen_id"
)

rna_cols <- setdiff(
  colnames(rna_scores),
  "specimen_id"
)

model_feature_sets <- list(
  Clinical = clinical_cols,
  Cytokine = cytokine_cols,
  RNA_GeneSet = rna_cols,
  Clinical_Cytokine = c(clinical_cols, cytokine_cols),
  Clinical_RNA = c(clinical_cols, rna_cols),
  Cytokine_RNA = c(cytokine_cols, rna_cols),
  Integrated_All = c(clinical_cols, cytokine_cols, rna_cols)
)

# ============================================================
# 8. RUN MODELS
# ============================================================

model_results <- list()

for (model_name in names(model_feature_sets)) {
  message("Running model: ", model_name)
  
  features <- intersect(
    model_feature_sets[[model_name]],
    colnames(integrated_data)
  )
  
  fit <- tryCatch(
    fit_loocv_elastic_net(
      df = integrated_data,
      feature_cols = features,
      alpha = 0.5
    ),
    error = function(e) {
      message("Model failed: ", model_name, " | ", e$message)
      NULL
    }
  )
  
  if (!is.null(fit)) {
    model_results[[model_name]] <- fit
  }
}

# ============================================================
# 9. SUMMARIZE MODEL PERFORMANCE
# ============================================================

performance_table <- bind_rows(
  lapply(names(model_results), function(model_name) {
    fit <- model_results[[model_name]]
    
    tibble(
      Model = model_name,
      AUC = fit$auc,
      Sensitivity = fit$sensitivity,
      Specificity = fit$specificity,
      Balanced_Accuracy = fit$balanced_accuracy
    )
  })
) %>%
  arrange(desc(AUC))

write.csv(
  performance_table,
  file.path(output_dir, "model_performance_loocv.csv"),
  row.names = FALSE
)

# Save predictions and selected feature frequencies
for (model_name in names(model_results)) {
  write.csv(
    model_results[[model_name]]$predictions,
    file.path(output_dir, paste0("predictions_", model_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_results[[model_name]]$feature_frequency,
    file.path(output_dir, paste0("selected_features_", model_name, ".csv")),
    row.names = FALSE
  )
}

# ============================================================
# 10. PLOT MODEL PERFORMANCE
# ============================================================

p_auc <- performance_table %>%
  mutate(Model = factor(Model, levels = rev(Model))) %>%
  ggplot(aes(x = Model, y = AUC)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  ylim(0, 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Cross-validated discrimination of SCAI D vs C",
    x = NULL,
    y = "Leave-one-out CV AUC"
  )

ggsave(
  file.path(output_dir, "model_auc_comparison.png"),
  p_auc,
  width = 7,
  height = 5,
  dpi = 600
)

# ============================================================
# 11. INTEGRATED SIGNATURE HEATMAP
# ============================================================

top_features <- model_results$Integrated_All$feature_frequency %>%
  slice_head(n = 20) %>%
  pull(feature)

heatmap_df <- integrated_data %>%
  dplyr::select(specimen_id, scai_stage, all_of(top_features))

heatmap_mat <- heatmap_df %>%
  dplyr::select(-specimen_id, -scai_stage) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_df$specimen_id

heatmap_mat <- as.matrix(heatmap_mat)
heatmap_mat <- apply(heatmap_mat, 2, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
})
heatmap_mat <- scale(heatmap_mat)

annotation_row <- data.frame(
  SCAI = heatmap_df$scai_stage
)
rownames(annotation_row) <- heatmap_df$specimen_id

pheatmap(
  t(heatmap_mat),
  annotation_col = annotation_row,
  filename = file.path(output_dir, "integrated_signature_heatmap.png"),
  width = 12,
  height = 6,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "predictive signature: SCAI D vs C"
)
cat("\nCross-modality prediction complete.\n")
cat("Results saved in:", output_dir, "\n")


