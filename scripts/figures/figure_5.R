setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Figure 5: Cross-modality integration
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pROC)
library(patchwork)
library(dplyr)
library(tibble)

# ---- Files ----
pca_scores_file <- ("cs_project/multi_modality/results/multimodal_pca_scores.csv")
roc_scores_file <- ("cs_project/multi_modality/results/multimodal_roc_scores.csv")
cross_modality_corr_file <- ("cs_project/multi_modality/results/integrated_multimodal_correlation_results.csv")
elastic_net_file <- ("cs_project/multi_modality/results/integrated_elastic_net_coeff_results.csv")
or_results_file <- ("cs_project/multi_modality/results/integrated_multimodal_selected_feature_ORs.csv")
loo_results_file <- ("cs_project/multi_modality/results/integrated_selected_features_LOOCV_predictions.csv")

# integrated_file <- "cs_project/integrated_prediction/results/integrated_feature_matrix.csv"
# plot_df_file <- "cs_project/association_integrated/results/plot_df.csv"
# corr_file <- "cs_project/association_integrated/results/corr_df.csv"
# performance_file <- "cs_project/integrated_prediction/results/model_performance_loocv.csv"
# 
# prediction_files <- list(
#   "Clinical" = "cs_project/integrated_prediction/results/predictions_Clinical.csv",
#   "Cytokine" = "cs_project/integrated_prediction/results/predictions_Cytokine.csv",
#   "RNA Gene set" = "cs_project/integrated_prediction/results/predictions_RNA_GeneSet.csv",
#   "Clinical + Cytokine" = "cs_project/integrated_prediction/results/predictions_Clinical_Cytokine.csv",
#   "Clinical + RNA" = "cs_project/integrated_prediction/results/predictions_Clinical_RNA.csv",
#   "Cytokine + RNA" = "cs_project/integrated_prediction/results/predictions_Cytokine_RNA.csv",
#   "Integrated all" = "cs_project/integrated_prediction/results/predictions_Integrated_All.csv"
# )

output_dir <- "cs_project/figures/Figure5_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dpi_out <- 600
panel_width <- 6.5

stage_colors <- c("Stage C" = "steelblue", "Stage D" = "darkred")
# modality_colors <- c("Clinical / Hemodynamic" = "#0072B2", "Cytokine" = "#D55E00", "RNA" = "#009E73", "Other" = "grey50")
modality_colors <- c("Clinical" = "#0072B2", "Cytokine" = "#D55E00", "RNA" = "#009E73")
pretty_feature <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_replace_all("admission labs ", "Admission: ") %>% str_replace_all("shock onset labs ", "Shock Onset: ") %>% 
    str_replace_all("level ", "") %>% str_replace_all("total ", "") %>% str_replace_all("mg dl", "mg/dL") %>% str_replace_all("mmol l", "mmol/L") %>%
    str_replace_all("ng l", "ng/L") %>% str_replace_all("g dl", "g/dL") %>% str_squish()
}

theme_pub <- function(base_size = 12) {theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10), axis.title = element_text(size = 11),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11), 
          # legend.position = "right", 
          legend.justification = "left", legend.title = element_blank())
}

# ============================================================
# Load data
# ============================================================
pca_scores <- read.csv(pca_scores_file)
roc_scores <- read.csv(roc_scores_file)

corr_results <- read.csv(cross_modality_corr_file)
coeff_results <- read.csv(elastic_net_file)
or_results <- read.csv(or_results_file)
# ============================================================
# Figure 5A: PCA based multi-modal
# ============================================================
plot_pca <- function(label) {
  df <- pca_scores %>% filter(modality == label) %>% mutate(Stage = factor(Stage, levels = c("C", "D"), labels = c("Stage C", "Stage D")))
  ggplot(df, aes(PC1, PC2, color = Stage)) + geom_point(size = 2.5, alpha = 0.9) + stat_ellipse(linewidth = 0.7) + scale_color_manual(values = stage_colors) + theme_pub() + 
    theme(legend.position = "none") + labs(title = label, x = paste0("PC1 (", unique(df$PC1_var), "%)"), y = paste0("PC2 (", unique(df$PC2_var), "%)"))
}

modality_order <- c("Clinical", "Cytokine", "RNA", "Clinical + Cytokine", "Clinical + RNA", "Cytokine + RNA", "Integrated")
roc_plot_df <- roc_scores %>% mutate(modality_auc = paste0(modality, " (AUC=", AUC, ")"), modality = factor(modality, levels = modality_order))
label_levels <- roc_plot_df %>% distinct(modality, modality_auc) %>% arrange(modality) %>% pull(modality_auc)
roc_plot_df <- roc_plot_df %>% mutate(modality_auc = factor(modality_auc, levels = label_levels))
p_roc <- ggplot(roc_plot_df, aes(1 - specificity, sensitivity, color = modality_auc)) + geom_path(linewidth = 0.9) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") + coord_equal() + theme_pub() + theme(legend.position = "right") + 
  labs(title = "Principle Component based ROC curves", x = "1 - Specificity", y = "Sensitivity", color = NULL)

p_fig5a <- (plot_pca("Clinical") + plot_pca("Cytokine") + plot_pca("RNA")) /
  (plot_pca("Clinical + Cytokine") + plot_pca("Clinical + RNA") + plot_pca("Cytokine + RNA")) /
  (plot_pca("Integrated") + p_roc + plot_spacer()) + 
  plot_annotation(title = "A. Principal component analysis", theme = theme(plot.title = element_text(face = "bold", size = 13)))
p_fig5a
# ggsave(filename = file.path(output_dir, "Figure5A_multimodal_pca_roc.png"), plot = p_fig5a, width = panel_width*2, height = 10, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 5B: Cross modality corr
# ============================================================
feature_order <- corr_results %>% dplyr::select(feature_1, modality_1) %>% rename(feature = feature_1, modality = modality_1) %>% 
  bind_rows(corr_results %>% dplyr::select(feature_2, modality_2) %>% rename(feature = feature_2, modality = modality_2)) %>% distinct() %>% 
  arrange(factor(modality, levels = c("Clinical", "Cytokine", "RNA")), feature) %>% pull(feature)
feature_labels <- corr_results %>% dplyr::select(feature_1, modality_1) %>% rename(feature = feature_1, modality = modality_1) %>%
  bind_rows(corr_results %>% dplyr::select(feature_2, modality_2) %>% rename(feature = feature_2, modality = modality_2)) %>% distinct() %>% 
  filter(feature %in% feature_order) %>% arrange(match(feature, feature_order))

x_colors <- modality_colors[feature_labels$modality]
y_colors <- rev(x_colors)
heatmap_df <- corr_results %>% mutate(feature_1 = factor(feature_1, levels = rev(feature_order)), feature_2 = factor(feature_2, levels = feature_order))
p_fig5b <- ggplot(heatmap_df, aes(feature_2, feature_1, fill = rho)) + geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-1, 1)) + scale_x_discrete(labels = pretty_feature) + scale_y_discrete(labels = pretty_feature) +
  theme_pub() + theme(legend.position = "right", axis.text.x = element_text(angle = 50, hjust = 1, size = 7.5, color = x_colors), axis.text.y = element_text(size = 9, color = y_colors)) +
  labs(title = "B. Multimodal feature correlation", x = "", y = "", fill = "Spearman ρ")
p_fig5b
# ggsave(file.path(output_dir, "Figure5B_integrated_multimodal_correlation.png"), p_fig5b, width = panel_width + 0.5, height = 6.5, dpi = dpi_out, bg = "white")
# ggsave(file.path(output_dir, "Figure5B_integrated_multimodal_correlation.pdf"), p_fig5b, width = panel_width + 0.5, height = 6.5, device = "pdf", bg = "white")

# ============================================================
# Figure 5C: Association with stage D
# ============================================================
coef_df <- coeff_results %>% mutate(feature_label = pretty_feature(feature), modality = factor(modality, levels = c("Clinical", "Cytokine", "RNA"))) %>%
  arrange(modality, desc(abs(coefficient))) %>% mutate(feature_label = factor(feature_label, levels = rev(feature_label)))
p_fig5c <- ggplot(coef_df, aes(coefficient, feature_label, fill = modality)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = modality_colors) + theme_pub() + theme(legend.position = "right") + labs(title = "C. Elastic-net multimodal association", x = "Coefficient", y = "", fill = NULL)
p_fig5c
# ggsave(file.path(output_dir, "Figure5C_elastic_net_coefficients.png"), p_fig5c, width = panel_width, height = 4, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 5D: Stage D OR
# ============================================================
or_df <- or_results %>% mutate(feature_label = pretty_feature(feature), modality = factor(modality, levels = c("Clinical", "Cytokine", "RNA"))) %>% arrange(modality, desc(OR)) %>%
  mutate(feature_label = factor(feature_label, levels = rev(feature_label)))
p_fig5d <- ggplot(or_df, aes(OR, feature_label, fill = modality)) + geom_col(width = 0.7) + geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") + scale_fill_manual(values = modality_colors) + theme_pub() +
  theme(legend.position = "none") + labs(title = "D. Elastic-net selected features", x = "Odds ratio (95% CI)", y = "", fill = NULL)
p_fig5d
# ggsave(file.path(output_dir, "Figure5D_multimodal_odds_ratio.png"), p_fig5d, width = panel_width, height = 4, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 5E: LOOCV ROC curve plot
# ============================================================
# loo_results_df <- read_csv(loo_results_file, show_col_types = FALSE)
# roc_obj_loo <- roc(loo_results_df$outcome, loo_results_df$pred_prob_loo, quiet = TRUE)
# auc_val_loo <- round(as.numeric(auc(roc_obj_loo)), 3)
# roc_5e_df <- tibble(specificity = roc_obj_loo$specificities, sensitivity = roc_obj_loo$sensitivities) %>% arrange(specificity)
# p_fig5e <- ggplot(roc_5e_df, aes(1 - specificity, sensitivity)) + geom_path(linewidth = 0.9, color = "#1f78b4") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") + coord_equal() + theme_pub() +
#   labs(title = paste0("E. LOOCV ROC for selected features \n (AUC=", sprintf("%.3f", auc_val_loo), ")"), x = "1 - Specificity", y = "Sensitivity")
# p_fig5e
# ggsave(file.path(output_dir, "Figure5E_LOOCV_ROC.png"), p_fig5e, width = panel_width, height = 2.5, dpi = dpi_out, bg = "white")

loo_results_df <- read_csv(loo_results_file, show_col_types = FALSE) %>% mutate(model_label = fct_reorder(model_label, auc, .desc = TRUE))
p_fig5e_compare <- ggplot(loo_results_df, aes(1 - specificity, sensitivity, color = model_label)) + geom_path(linewidth = 0.9) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  coord_equal(expand = FALSE, xlim = c(0, 1), ylim = c(0, 1)) + theme_pub() + labs(title = "E. LOOCV ROC for selected features", x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme(legend.position = "right", legend.text = element_text(size = 8))
p_fig5e_compare
ggsave(file.path(output_dir, "Figure5E_LOOCV_ROC_by_modality.png"), p_fig5e_compare, width = panel_width, height = 3, dpi = dpi_out, bg = "white")

