setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Cross-Modality Analysis SCAI Stage D vs C
# Clinical + Transcriptomic + Cytokine
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
library(patchwork)
library(logistf)

output_dir <- "cs_project/multi_modality/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stage_colors <- c("Stage C" = "steelblue", "Stage D" = "darkred")
direction_colors <- c("Higher in Stage C" = "steelblue", "Higher in Stage D" = "darkred", "Not significant" = "gray75")
theme_pub <- function(base_size = 12) {theme_classic(base_size = base_size) + 
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9, hjust = 0),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10), axis.title = element_text(size = 11), 
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
          legend.position = "none", legend.justification = "left", legend.title = element_blank())
}


# 1. Cinical
clinical_features <- read.csv("cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv") %>%
  mutate(specimen_id = str_replace(specimen_id, "^SH-", "SR-"), scai_stage = toupper(str_trim(scai_stage))) %>% filter(scai_stage %in% c("C", "D")) %>%
  dplyr::select(specimen_id, matches("^(admission|shock_onset)_labs_(lactate|glucose|albumin|calcium|nt_pro_bnp|troponin)")) %>%
  mutate(specimen_id = specimen_id %>% str_replace("^SR-24-0", "SR24") %>% str_replace("^SR-25-0", "SR25") %>% str_replace("-", ""), across(-specimen_id, ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::select(specimen_id, where(~ mean(is.na(.)) < 0.90)) %>% mutate(across(where(is.numeric), ~ replace_na(.x, median(.x, na.rm = TRUE))))

# 2. Cytokines
cytokine_features <- read_excel("cs_project/cytokine_data/CScytokineData.xlsx") %>%
  rename_with(~ str_replace_all(., "-", "")) %>% rename_with(~ str_remove_all(., " \\(pg/ml\\)")) %>% rename_with(~ str_replace_all(., " ", "")) %>%
  mutate(SCAIStage = toupper(SCAIStage), specimen_id = str_trim(SRspecimenID)) %>% filter(SCAIStage %in% c("C", "D")) %>%
  dplyr::select(specimen_id, any_of(c("IL6", "ST2", "VEGFA", "S100A8", "IL8", "TroponinI"))) %>%
  mutate(specimen_id = specimen_id %>% str_replace("^SR-24-0", "SR24") %>% str_replace("^SR-25-0", "SR25") %>% str_replace("-", ""), across(-specimen_id, ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::select(specimen_id, where(~ mean(is.na(.) | . == 0) < 0.90)) %>%
  mutate(across(where(is.numeric), ~ replace(.x, is.na(.x) | .x == 0, median(.x[!is.na(.x) & .x > 0], na.rm = TRUE))))

# 3. Transcriptomics
metadata <- read.csv("cs_project/bulk_data/metadata_bulk.csv") %>% mutate(specimen_id = str_trim(SR.specimen.ID), scai_stage = toupper(str_trim(SCAI.Stage))) %>%
  filter(scai_stage %in% c("C", "D")) %>% mutate(outcome = ifelse(scai_stage == "D", 1, 0)) %>% dplyr::select(specimen_id, scai_stage, outcome)
mapping <- read.csv("cs_project/gene_info_biomart.csv") %>% dplyr::select(gene_id, gene_name) %>% distinct(gene_id, .keep_all = TRUE)
genes <- c("PEX11G", "ARG1", "FN1", "OSBPL1A", "SEC14L5", "RSPH14",
           "KLKB1", "FGFR2", "LRRN3", "TRGV10", "ACOT11")
gene_ids <- mapping %>% filter(gene_name %in% genes) %>% distinct(gene_id, gene_name)

load("cs_project/bulk_data/count_data.RData")
expr <- 2^log2.tpm.matrix - 1
expr <- expr[rownames(expr) %in% gene_ids$gene_id, intersect(colnames(expr), metadata$specimen_id), drop = FALSE]
rna_features <- expr %>% as.data.frame() %>% rownames_to_column("gene_id") %>% left_join(gene_ids, by = "gene_id") %>% dplyr::select(-gene_id) %>% group_by(gene_name) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% column_to_rownames("gene_name") %>% as.matrix() %>% t() %>% as.data.frame() %>%
  rownames_to_column("specimen_id")

# 4. Merge modalities
# summary(clinical_features)
# summary(cytokine_features)
# summary(rna_features)

id_cols <- c("specimen_id", "scai_stage", "outcome")
integrated_data <- metadata %>% left_join(clinical_features, by = "specimen_id") %>% left_join(cytokine_features, by = "specimen_id") %>% left_join(rna_features, by = "specimen_id") %>%
  mutate(across(-all_of(id_cols), ~ suppressWarnings(as.numeric(.x)))) %>% mutate(across(-all_of(id_cols), ~ replace_na(.x, median(.x, na.rm = TRUE)))) %>% mutate(across(-all_of(id_cols), ~ as.numeric(scale(log1p(.x)))))
# summary(integrated_data)
# write.csv(integrated_data, file.path(output_dir, "multimodal_integrated_scaled_dataset.csv"), row.names = FALSE)

# 5. PCA based separation of C and D
clinical_vars <- clinical_features %>% dplyr::select(-specimen_id) %>% names()
cytokine_vars <- cytokine_features %>% dplyr::select(-specimen_id) %>% names()
rna_vars <- rna_features %>% dplyr::select(-specimen_id) %>% names()

modalities <- list(
  "Clinical" = clinical_vars,
  "Cytokine" = cytokine_vars,
  "RNA" = rna_vars,
  "Clinical + Cytokine" = c(clinical_vars, cytokine_vars),
  "Clinical + RNA" = c(clinical_vars, rna_vars),
  "Cytokine + RNA" = c(cytokine_vars, rna_vars),
  "Integrated" = c(clinical_vars, cytokine_vars, rna_vars)
)

# get_pca_roc <- function(vars, label) {
#   pca <- prcomp(integrated_data %>% dplyr::select(all_of(vars)) %>% as.matrix(), center = FALSE, scale. = FALSE)
#   var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
#   pca_df <- as.data.frame(pca$x[, 1:2]) %>% mutate(specimen_id = integrated_data$specimen_id, Stage = integrated_data$scai_stage, modality = label, PC1_var = var_exp[1], PC2_var = var_exp[2])
#   fit <- glm(integrated_data$outcome ~ PC1 + PC2, data = pca_df, family = binomial())
#   roc_obj <- roc(integrated_data$outcome, predict(fit, type = "response"), quiet = TRUE)
#   roc_df <- tibble(specificity = roc_obj$specificities, sensitivity = roc_obj$sensitivities, modality = label, AUC = round(as.numeric(auc(roc_obj)), 3))
#   list(pca = pca_df, roc = roc_df)
# }
# res <- imap(modalities, get_pca_roc)
# pca_scores <- bind_rows(map(res, "pca"))
# roc_scores <- bind_rows(map(res, "roc"))

get_pca_roc <- function(vars, label) {
  pca <- prcomp(integrated_data %>% dplyr::select(all_of(vars)) %>% as.matrix(), center = FALSE, scale. = FALSE)
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  pca_df <- as.data.frame(pca$x[, 1:2]) %>% mutate(specimen_id = integrated_data$specimen_id, Stage = integrated_data$scai_stage, modality = label, PC1_var = var_exp[1], PC2_var = var_exp[2])
  fit <- glm(integrated_data$outcome ~ PC1 + PC2, data = pca_df, family = binomial())
  roc_obj <- roc(integrated_data$outcome, predict(fit, type = "response"), quiet = TRUE)
  ci_obj <- ci.auc(roc_obj)
  roc_df <- tibble(specificity = roc_obj$specificities, sensitivity = roc_obj$sensitivities, modality = label, 
                   AUC = round(as.numeric(auc(roc_obj)), 3),
                   CI_low = round(as.numeric(ci_obj[1]), 3),
                   CI_high = round(as.numeric(ci_obj[3]), 3))
  list(pca = pca_df, roc = roc_df)
}
res <- imap(modalities, get_pca_roc)
pca_scores <- bind_rows(map(res, "pca"))
roc_scores <- bind_rows(map(res, "roc"))

# write.csv(pca_scores, file.path(output_dir, "multimodal_pca_scores.csv"), row.names = FALSE)
# write.csv(roc_scores, file.path(output_dir, "multimodal_roc_scores.csv"), row.names = FALSE)

# 7. Correlation analysis
integrated_data <- read.csv("cs_project/multi_modality/results/multimodal_integrated_scaled_dataset.csv")
clinical_vars <- c("admission_labs_glucose_level_mg_dl", "admission_labs_calcium_level_total_mg_dl", "admission_labs_albumin_level_g_dl",
                  "admission_labs_nt_pro_bnp_ng_l", "admission_labs_lactate_mmol_l", "shock_onset_labs_glucose_level_mg_dl",
                  "shock_onset_labs_calcium_level_total_mg_dl", "shock_onset_labs_albumin_level_g_dl", "shock_onset_labs_nt_pro_bnp_ng_l", "shock_onset_labs_lactate_mmol_l")
cytokine_vars <- c("IL6", "ST2", "VEGFA", "S100A8", "IL8", "TroponinI")
rna_vars <- c("ACOT11", "ARG1", "FGFR2", "FN1", "KLKB1", "LRRN3", "OSBPL1A", "PEX11G", "RSPH14", "SEC14L5", "TRGV10" )

features <- c(clinical_vars, cytokine_vars, rna_vars)
feature_modality <- tibble(feature = features, modality = case_when(feature %in% clinical_vars ~ "Clinical", feature %in% cytokine_vars ~ "Cytokine", feature %in% rna_vars ~ "RNA"))

corr_results <- expand.grid(feature_1 = features, feature_2 = features, stringsAsFactors = FALSE) %>% rowwise() %>%
  mutate(rho = cor(integrated_data[[feature_1]], integrated_data[[feature_2]], method = "spearman"), 
         p_value = ifelse(feature_1 == feature_2, 0, cor.test(integrated_data[[feature_1]], integrated_data[[feature_2]], method = "spearman", exact = FALSE)$p.value)) %>% ungroup() %>% 
  left_join(feature_modality, by = c("feature_1" = "feature")) %>% rename(modality_1 = modality) %>%
  left_join(feature_modality, by = c("feature_2" = "feature")) %>% rename(modality_2 = modality) %>%
  mutate(modality_pair = paste(modality_1, modality_2, sep = " - "), FDR = p.adjust(p_value, method = "BH"), abs_rho = abs(rho))
# write.csv(corr_results, file.path(output_dir, "integrated_multimodal_correlation_results.csv"), row.names = FALSE)

# 8. Association studies
library(glmnet)
x <- integrated_data %>% dplyr::select(-specimen_id, -scai_stage, -outcome) %>% as.matrix()
y <- integrated_data$outcome
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5, nfolds = nrow(x))
fit <- glmnet(x, y, family = "binomial", alpha = 0.5, lambda = cv_fit$lambda.1se)

coef_df <- as.matrix(coef(fit)) %>% as.data.frame() %>% rownames_to_column("feature") %>% rename(coefficient = s0) %>% filter(feature != "(Intercept)", coefficient != 0) %>%
  mutate(modality = case_when(feature %in% clinical_vars ~ "Clinical", feature %in% cytokine_vars ~ "Cytokine", feature %in% rna_vars ~ "RNA", TRUE ~ "Other")) %>%
  arrange(desc(abs(coefficient)))
coef_df
# write.csv(coef_df, file.path(output_dir, "integrated_elastic_net_coeff_results.csv"), row.names = FALSE)

# 9. Standard Logistic Regression
selected_features <- coef_df$feature
or_results <- purrr::map_dfr(selected_features, function(feat) {
  fit <- glm(as.formula(paste("outcome ~", feat)), data = integrated_data, family = binomial())
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == feat) %>% transmute(feature = feat, OR = estimate, CI_low = conf.low, CI_high = conf.high, p_value = p.value)
}) %>% mutate(FDR = p.adjust(p_value, method = "BH"), OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)) %>% left_join(coef_df %>% dplyr::select(feature, modality), by = "feature")
or_results
# write.csv(or_results, file.path(output_dir, "integrated_multimodal_selected_feature_ORs.csv"), row.names = FALSE)

# # 10. Leave-one-out cross-validated ROC/AUC for selected features
# selected_features <- or_results$feature
# model_formula <- as.formula(paste("outcome ~", paste(selected_features, collapse = " + ")))
# n <- nrow(integrated_data)
# loo_pred <- numeric(n)
# for (i in seq_len(n)) {
#   train_dat <- integrated_data[-i, ]
#   test_dat  <- integrated_data[i, ]
#   fit_loo <- logistf(model_formula, data = train_dat)
#   loo_pred[i] <- predict(fit_loo, newdata = test_dat, type = "response")
# }
# 
# loo_results_df <- integrated_data %>% dplyr::select(specimen_id, outcome) %>% mutate(pred_prob_loo = loo_pred)
# # write.csv(loo_results_df, file.path(output_dir, "integrated_selected_features_LOOCV_predictions.csv"), row.names = FALSE)

modality_sets <- list(Clinical = selected_features[selected_features %in% clinical_vars], Cytokine = selected_features[selected_features %in% cytokine_vars],
                      RNA = selected_features[selected_features %in% rna_vars], Combined = selected_features)
n <- nrow(integrated_data)
run_loocv <- function(feats) {
  form <- as.formula(paste("outcome ~", paste(feats, collapse = " + "))); pred <- numeric(n)
  for (i in seq_len(n)) {train_dat <- integrated_data[-i, ]; test_dat  <- integrated_data[i, ]; fit_loo <- logistf(form, data = train_dat); pred[i] <- predict(fit_loo, newdata = test_dat, type = "response")}
  pred
}
roc_all <- map_dfr(names(modality_sets), function(nm) {
  feats <- modality_sets[[nm]]; pred <- run_loocv(feats); roc_obj <- roc(integrated_data$outcome, pred, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj)); ci_val <- as.numeric(ci.auc(roc_obj))  # c(lower, auc, upper)
  tibble(model = nm, auc = auc_val, ci_low = ci_val[1], ci_high = ci_val[3], specificity = roc_obj$specificities, sensitivity = roc_obj$sensitivities)
}) %>% mutate(model_label = fct_reorder(sprintf("%s (AUC=%.3f, 95%% CI %.3f - %.3f)", model, auc, ci_low, ci_high), auc, .desc = TRUE))
write.csv(roc_all, file.path(output_dir, "integrated_selected_features_LOOCV_predictions.csv"), row.names = FALSE)