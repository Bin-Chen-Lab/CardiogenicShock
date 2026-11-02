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

output_dir <- "cs_project/association_integrated/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv("cs_project/bulk_data/metadata_bulk.csv") %>% mutate(specimen_id = str_trim(SR.specimen.ID), scai_stage = toupper(str_trim(SCAI.Stage))) %>%
  filter(scai_stage %in% c("C", "D")) %>% mutate(outcome = ifelse(scai_stage == "D", 1, 0)) %>% dplyr::select(specimen_id, scai_stage, outcome)

# 1. CLINICAL FEATURES
clinical <- read.csv("cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv") %>% 
  mutate(specimen_id = str_replace(specimen_id, "^SH-", "SR-"), scai_stage = toupper(str_trim(scai_stage))) %>% filter(scai_stage %in% c("C", "D"))

clinical_features <- clinical %>% 
  dplyr::select(specimen_id,
                # admission/metabolic/end-organ dysfunction
                matches("admission_labs_lactate"), matches("admission_labs_bicarbonate"), matches("admission_labs_creatinine"),
                matches("admission_labs_mdrd_egfr"), matches("admission_labs_albumin"), matches("admission_labs_glucose"),
                matches("admission_labs_nt_pro_bnp"), matches("admission_labs_troponin"), matches("admission_labs_c_reactive_protein"),
                # shock onset labs
                matches("shock_onset_labs_lactate"), matches("shock_onset_labs_bicarbonate"), matches("shock_onset_labs_creatinine"),
                matches("shock_onset_labs_mdrd_egfr"), matches("shock_onset_labs_albumin"), matches("shock_onset_labs_glucose"),
                matches("shock_onset_labs_nt_pro_bnp"), matches("shock_onset_labs_troponin"), matches("shock_onset_labs_c_reactive_protein"),
                # hemodynamics
                matches("shock_onset_rhc_fick_ci"), matches("shock_onset_rhc_fick_co"), matches("shock_onset_rhc_svr"), matches("shock_onset_rhc_pvr"),
                matches("shock_onset_rhc_pcwp"), matches("shock_onset_pac_fick_cardiac_index"), matches("shock_onset_pac_fick_cardiac_output"),
                matches("shock_onset_pac_systemic_vascular_resistance"), matches("shock_onset_pac_pulmonary_artery_pulsatility_index"))

clinical_features <- clinical_features %>% mutate(across(-specimen_id, ~ suppressWarnings(as.numeric(.))))
# Remove variables with too much missingness
clinical_features <- clinical_features %>% dplyr::select(specimen_id, where(~ mean(is.na(.)) < 0.50)) %>%
  mutate(specimen_id = specimen_id %>% str_replace("^SR-24-0", "SR24") %>% str_replace("^SR-25-0", "SR25") %>% str_replace("-", ""))


# 2. CYTOKINE FEATURES
cytokine_raw <- read_excel("cs_project/cytokine_data/CScytokineData.xlsx") %>% rename_with(~ str_replace_all(., "-", "")) %>%
  rename_with(~ str_remove_all(., " \\(pg/ml\\)")) %>% rename_with(~ str_replace_all(., " ", "")) %>%
  mutate(SCAIStage = toupper(SCAIStage), specimen_id = str_trim(SRspecimenID)) %>% filter(SCAIStage %in% c("C", "D"))
meta_cols <- c("StudyID", "Age", "Sex", "Race", "SCAIStage", "SRspecimenID", "specimen_id")
cytokine_features <- cytokine_raw %>% dplyr::select(specimen_id, where(is.numeric))
cytokine_features <- cytokine_features %>% mutate(across(-specimen_id, ~ log1p(.)))
cytokine_features <- cytokine_features %>% dplyr::select(specimen_id, where(~ mean(is.na(.)) < 0.50), where(~ length(unique(na.omit(.))) > 3))
priority_cytokines <- c("ST2", "IL6", "IL8", "S100A8", "VEGFA")
cytokine_features <- cytokine_features %>% dplyr::select(specimen_id, any_of(priority_cytokines))
cytokine_features <- cytokine_features %>% mutate(specimen_id = specimen_id %>% str_replace("^SR-24-0", "SR24") %>% str_replace("^SR-25-0", "SR25") %>% str_replace("-", ""))

# 3.RNA EXPRESSION
load("cs_project/bulk_data/count_data.RData")
expr <- log2.tpm.matrix
mapping <- read.csv("cs_project/gene_info_biomart.csv") %>% dplyr::select(gene_id, gene_name) %>% distinct(gene_id, .keep_all = TRUE)
expr <- expr[, intersect(colnames(expr), metadata$specimen_id)]
expr_df <- as.data.frame(expr) %>% rownames_to_column("gene_id") %>% left_join(mapping, by = "gene_id") %>% filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
expr_mat <- expr_df %>% column_to_rownames("gene_name") %>% as.matrix()

gene_sets <- c("PEX11G", "ARG1", "FN1", "OSBPL1A", "SEC14L5",
               "KLKB1", "FGFR2", "LRRN3", "TRGV10", "ACOT11")
rna_features <- as.data.frame(t(expr_mat[gene_sets, ]))
rna_features$specimen_id <- rownames(rna_features)


# 4. MERGE MODALITIES
integrated_data <- metadata %>% left_join(clinical_features, by = "specimen_id") %>% left_join(cytokine_features, by = "specimen_id") %>% left_join(rna_features, by = "specimen_id")
# colnames(integrated_data)
# head(integrated_data, 2)
write.csv2(integrated_data, file.path(output_dir, "association_integrated_feature_matrix.csv"), row.names = FALSE)

# 5. Prepare data and Median imputation
library(tidyverse)
library(glmnet)
library(pROC)

model_data <- integrated_data %>% dplyr::select(-specimen_id, -scai_stage)
y <- as.numeric(as.character(model_data$outcome))
x_data <- model_data %>% dplyr::select(-outcome)

x_data <- x_data %>% mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) # median
x_data <- x_data %>% mutate(across(everything(), ~ as.numeric(scale(.x)))) # scaling
x <- as.matrix(x_data)

# 6. Elastic Net + LOOCV
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5, nfolds = nrow(x), standardize = FALSE, type.measure = "deviance")
plot(cvfit)

lambda_min <- cvfit$lambda.min; lambda_min
lambda_1se <- cvfit$lambda.1se; lambda_1se

# at lambda.min
coef_min <- coef(cvfit, s = "lambda.min")
coef_min_df <- as.matrix(coef_min) %>% as.data.frame() %>% rownames_to_column("marker") %>% rename(coefficient = s1) %>%  filter(coefficient != 0) %>% arrange(desc(abs(coefficient)))
coef_min_df

# at lambda.1se
coef_1se <- coef(cvfit, s = "lambda.1se")
coef_1se_df <- as.matrix(coef_1se) %>% as.data.frame() %>% rownames_to_column("marker") %>% rename(coefficient = s1) %>% filter(coefficient != 0) %>% arrange(desc(abs(coefficient)))
coef_1se_df


# 6. Plot data
clinical <- c("admission_labs_lactate_mmol_l", "admission_labs_bicarbonate_level_mmol_l", "admission_labs_creatinine_mg_dl", "admission_labs_mdrd_egfr",
              "admission_labs_albumin_level_g_dl", "admission_labs_glucose_level_mg_dl", "admission_labs_nt_pro_bnp_ng_l",
              "shock_onset_labs_bicarbonate_level_mmol_l", "shock_onset_labs_creatinine_mg_dl", "shock_onset_labs_mdrd_egfr", "shock_onset_labs_albumin_level_g_dl",
              "shock_onset_labs_glucose_level_mg_dl", "shock_onset_rhc_fick_ci_l_min_bsa", "shock_onset_rhc_fick_co_l_min", "shock_onset_rhc_pcwp_mmhg",
              "shock_onset_pac_fick_cardiac_index_ci_l_min", "shock_onset_pac_fick_cardiac_output_co_l_min", "shock_onset_pac_systemic_vascular_resistance_svr_dynes_sec_cm_5",
              "shock_onset_pac_pulmonary_artery_pulsatility_index_papi")

cytokine <- c("ST2", "IL6", "IL8", "S100A8", "VEGFA")
rna <- c("PEX11G", "ARG1", "FN1", "OSBPL1A", "SEC14L5", "KLKB1", "FGFR2", "LRRN3", "TRGV10", "ACOT11")

coef_min_plot <- coef_min_df %>% mutate(model = "lambda.min")
coef_1se_plot <- coef_1se_df %>% mutate(model = "lambda.1se")
plot_data <- bind_rows(coef_min_df %>% mutate(model = "lambda.min"), coef_1se_df %>% mutate(model = "lambda.1se")) %>% filter(marker != "(Intercept)") %>%
  mutate(feature_group = case_when(marker %in% clinical ~ "Clinical / Hemodynamic", marker %in% cytokine ~ "Cytokine", marker %in% rna ~ "RNA", TRUE ~ "Other"),
         feature_group = factor(feature_group, levels = c("Clinical / Hemodynamic", "Cytokine", "RNA", "Other")), direction = if_else(coefficient > 0, "Positive", "Negative")) %>%
  arrange(model, feature_group, coefficient) %>% mutate(marker_order = factor(marker, levels = unique(marker)))

p_association <- ggplot(plot_data, aes(marker_order, coefficient, fill = feature_group)) + geom_hline(yintercept = 0, linewidth = 0.4) + coord_flip() +
  geom_col(width = 0.75) + geom_text(aes(label = round(coefficient, 2)), hjust = ifelse(plot_data$coefficient > 0, -0.15, 1.15), size = 3.5) +
  facet_wrap(~ model, scales = "free_y") + scale_fill_manual(values = c("Clinical / Hemodynamic" = "steelblue4", "Cytokine" = "darkorange3", "RNA" = "darkgreen", "Other" = "grey50")) +
  theme_minimal(base_size = 14) + labs(title = "Association of markers with severity (Elastic Net)", x = NULL, y = "Standardized coefficient", fill = "Feature group") +
  expand_limits(y = range(plot_data$coefficient) * 1.25)
write.csv(plot_data, 'cs_project/association_integrated/results/plot_df.csv')


ggsave(filename = "cs_project/association_integrated/results/elastic_net_integrated_biomarker_signatures.png", plot = p_association, width = 18, height = 10, dpi = 300)
ggsave(filename = "cs_project/association_integrated/results/elastic_net_integrated_biomarker_signatures.pdf", plot = p_association, width = 18, height = 10)


# 7. Why some features are omitted?
corr_data <- integrated_data %>% dplyr::select(-specimen_id, -outcome, -scai_stage)
corr_data <- corr_data %>% mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))
corr_matrix <- cor(corr_data, method = "spearman", use = "pairwise.complete.obs")

corr_df <- as.data.frame(as.table(corr_matrix))
colnames(corr_df) <- c("Feature1", "Feature2", "Correlation")
write.csv(corr_df, 'cs_project/association_integrated/results/corr_df.csv')

p_corr <- ggplot(corr_df, aes(x = Feature1, y = Feature2, fill = Correlation)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-1, 1)) + coord_fixed() + theme_minimal(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 9), panel.grid = element_blank()) +
  labs(title = "Integrated Feature Correlation Matrix", x = NULL, y = NULL, fill = "Spearman\nCorrelation")

ggsave(filename = "cs_project/association_integrated/results/integrated_feature_correlation_heatmap.png", plot = p_corr, width = 14, height = 12, dpi = 300)
ggsave(filename = "cs_project/association_integrated/results/integrated_feature_correlation_heatmap.pdf", plot = p_corr, width = 14, height = 12)