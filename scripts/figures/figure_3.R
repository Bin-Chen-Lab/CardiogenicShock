setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Figure 3: Bulk RNA-seq analysis of SCAI Stage D vs C shock
# ============================================================
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggnewscale)
library(dplyr)
library(enrichplot)
library(patchwork)
library(ggplotify)

# ---- Files ----
expr_file <- "cs_project/bulk_data/count_data.RData"
metadata_file <- "cs_project/bulk_data/metadata_bulk.csv"
gene_info_file <- "cs_project/gene_info_biomart.csv"
dge_file <- "cs_project/bulk_data/results/dge_results_annotated.csv"
dge_sig_file <- "cs_project/bulk_data/results/dge_results_significant_1.csv"
targeted_file <- "cs_project/bulk_data/results/targeted_DE_genes_by_clinical_axis.csv"
pathway_file <- "cs_project/bulk_data/results/pathway_Enrichment_ORA.xlsx"
custom_pathawy_file <- "cs_project/bulk_data/results/gsea_res_custom_metabolic.rds"
custom_pathawy_gsea_tpm_res <- "cs_project/bulk_data/results/ssgsea_tpm_results.rds"
ciber_file <- "cs_project/bulk_data/CIBERSORTx_Job3_Results.csv"
co_exp_file <- "cs_project/bulk_data/results/wgcna_results.rds"

output_dir <- "cs_project/figures/Figure3_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Settings ----
dpi_out <- 600
panel_width <- 6.5
stage_colors <- c("Stage C" = "steelblue", "Stage D" = "darkred")
direction_colors <- c("Higher in Stage C" = "steelblue", "Higher in Stage D" = "darkred", "Not significant" = "gray75")
pathway_colors <- c("Up_sig" = "darkred", "Up_ns" = "#d9a3a3", "Down_sig" = "steelblue", "Down_ns" = "#b8cfe0")
theme_pub <- function(base_size = 12) {theme_classic(base_size = base_size) + 
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9, hjust = 0),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10), axis.title = element_text(size = 11),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11), 
          legend.position = "none", legend.justification = "left", legend.title = element_blank())
}
theme_pub_small <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = 7), plot.subtitle = element_text(size = 6), plot.caption = element_text(size = 5, hjust = 0),
          axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 6), axis.title = element_text(size = 6),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6), legend.position = "none")
}

sig_star <- function(p) {dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}

# ============================================================
# Load data
# ============================================================
load(expr_file)
metadata <- read.csv(metadata_file)
gene_info <- read.csv(gene_info_file)
metadata_cd <- metadata %>% dplyr::mutate(sample_id = SR.specimen.ID, scai_stage = SCAI.Stage, Stage = factor(scai_stage, levels = c("C", "D"), labels = c("Stage C", "Stage D"))) %>%
  dplyr::filter(scai_stage %in% c("C", "D"))
common_samples <- intersect(colnames(log2.tpm.matrix), metadata_cd$sample_id)
expr_cd <- log2.tpm.matrix[, common_samples, drop = FALSE]
metadata_cd <- metadata_cd %>% dplyr::filter(sample_id %in% common_samples) %>% dplyr::arrange(match(sample_id, colnames(expr_cd)))
expr_cd <- expr_cd[, metadata_cd$sample_id, drop = FALSE]
gene_annot <- gene_info %>% dplyr::select(gene_id, gene_name) %>% dplyr::filter(!is.na(gene_name), gene_name != "") %>% dplyr::distinct(gene_id, .keep_all = TRUE)

# ============================================================
# Figure 3A: Heatmap of top DE genes
# ============================================================
dge_sig <- read.csv(dge_sig_file)
heat_labels <- dge_sig %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 40) %>% dplyr::mutate(gene_label = dplyr::if_else(is.na(gene_name) | gene_name == "", gene, gene_name))
heat_mat <- expr_cd[heat_labels$gene, , drop = FALSE]
rownames(heat_mat) <- heat_labels$gene_label
# # heat_mat_z <- heat_mat
# heat_mat_z <- t(scale(t(heat_mat)))
# heat_mat_z <- heat_mat_z[complete.cases(heat_mat_z), , drop = FALSE]
ann_df <- metadata_cd %>% dplyr::filter(sample_id %in% colnames(heat_mat_z)) %>% dplyr::select(sample_id, Stage)
cluster_ids <- function(ids) {
  if (length(ids) <= 1) return(ids)
  hclust(dist(t(heat_mat_z[, ids, drop = FALSE])))$labels[hclust(dist(t(heat_mat_z[, ids, drop = FALSE])))$order]
}
row_order <- hclust(dist(heat_mat_z))$labels[hclust(dist(heat_mat_z))$order]
sample_order <- ann_df %>% dplyr::mutate(Stage = factor(Stage, levels = c("Stage C", "Stage D"))) %>% dplyr::arrange(Stage) %>% dplyr::group_by(Stage) %>%
  dplyr::summarise(sample_id = list(cluster_ids(sample_id)), .groups = "drop") %>% dplyr::pull(sample_id) %>% unlist()
heat_df <- heat_mat_z[row_order, sample_order, drop = FALSE] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% 
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "Expression") %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(row_order)), sample_id = factor(sample_id, levels = sample_order))
stage_df <- ann_df %>% dplyr::mutate(sample_id = factor(sample_id, levels = sample_order), y = "SCAI")

p_fig3a <- ggplot() + geom_tile(data = stage_df, aes(sample_id, y, fill = Stage), height = 0.8) + scale_fill_manual(values = stage_colors, name = "SCAI") +
  ggnewscale::new_scale_fill() + geom_tile(data = heat_df, aes(sample_id, gene, fill = Expression)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Expression") +
  labs(title = "A. Transcriptomic heterogeneity", x = NULL, y = NULL) + theme_pub() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(), 
        legend.position = "right", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9))
p_fig3a
ggsave(filename = file.path(output_dir, "Figure3A_top_DEG_heatmap_grouped_by_stage_testII.png"), plot = p_fig3a, width = panel_width, height = 6.5, dpi = dpi_out, bg = "white")


# ============================================================
# Figure 3B: Differentially expressed genes volcano plot
# ============================================================
dge <- read.csv(dge_file)
volcano_dat <- dge %>% dplyr::mutate(gene_label = dplyr::case_when(is.na(gene_name) | gene_name == "" ~ gene, TRUE ~ gene_name), neg_log10_fdr = -log10(padj),
                                     direction = dplyr::case_when(!is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Higher in Stage D", 
                                                                  !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Higher in Stage C", TRUE ~ "Not significant"))
top_labels <- volcano_dat %>% dplyr::filter(!is.na(padj), padj < 0.05) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 20)

p_fig3b <- ggplot(volcano_dat, aes(x = log2FoldChange, y = neg_log10_fdr, color = direction)) + geom_point(size = 1.5, alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray45") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "gray45") +
  ggrepel::geom_text_repel(data = top_labels, aes(label = gene_label), size = 3.5, box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = direction_colors) + labs(title = "B. Differentially expressed genes", x = "log2 fold-change: Stage D vs Stage C", y = "-log10(FDR)") + theme_pub()
p_fig3b
ggsave(filename = file.path(output_dir, "Figure3B_differential_expression_volcano.png"), plot = p_fig3b, width = panel_width, height = 5.2, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 3C: Relevant cardiogenic shock genes forest-style plot
# Use ENSG ID if gene_name is NA
# ============================================================
dge_all <- read.csv(dge_file)

gene_comments <- tibble::tribble(
  ~gene_label, ~comment,
  "GDF15", "myocardial stress",
  "IL1RL1", "cardiac remodeling",
  "ADM", "vascular stress",
  "ANGPT2", "endothelial injury",
  "VWF", "endothelial activation",
  "ICAM1", "vascular inflammation",
  "VCAM1", "vascular inflammation",
  "SELE", "activated endothelium",
  "HIF1A", "hypoxia response",
  "VEGFA", "ischemia signaling",
  "EDN1", "vasoconstriction",
  "NLRP3", "inflammasome",
  "IL1B", "sterile inflammation",
  "TNF", "systemic inflammation",
  "CXCL8", "neutrophil recruitment",
  "S100A8", "neutrophil activation",
  "S100A9", "innate immunity",
  "MMP9", "matrix remodeling",
  "LCN2", "tissue injury",
  "PPBP", "platelet activation",
  "PF4", "platelet degranulation",
  "THBS1", "platelet-endothelial interaction",
  "TGFB1", "fibrosis",
  "TIMP1", "matrix regulation",
  "FN1", "fibrosis/remodeling"
)

relevant_de <- dge_all %>% dplyr::mutate(gene_label = dplyr::if_else(is.na(gene_name) | gene_name == "", gene, gene_name), direction = dplyr::if_else(log2FoldChange > 0, "Higher in Stage D", "Higher in Stage C")) %>%
  dplyr::inner_join(gene_comments, by = "gene_label") %>% dplyr::mutate(gene_plot = forcats::fct_reorder(paste0("**", gene_label, "**", " (", comment, ")"), log2FoldChange))
setdiff(gene_comments$gene_label, relevant_de$gene_label)
p_fig3c <- ggplot(relevant_de, aes(log2FoldChange, gene_plot, color = direction)) + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray45") +
  geom_segment(aes(x = 0, xend = log2FoldChange, yend = gene_plot), linewidth = 0.8) + geom_point(size = 3) + scale_color_manual(values = direction_colors) +
  labs(title = "C. Relevant genes", x = "log2 fold-change: Stage D vs Stage C", y = NULL) + 
  theme_pub() + theme(axis.text.y = ggtext::element_markdown(size = 10), legend.position = "none")
p_fig3c
ggsave(file.path(output_dir, "Figure3C_relevant_CS_genes_forest.png"), p_fig3c, width = panel_width, height = 6, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 3C: Heatmap of top DE genes
# ============================================================
dge_sig <- read.csv(dge_sig_file)
heat_labels <- dge_sig %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 40) %>% dplyr::mutate(gene_label = dplyr::if_else(is.na(gene_name) | gene_name == "", gene, gene_name))
heat_mat <- expr_cd[heat_labels$gene, , drop = FALSE]
rownames(heat_mat) <- heat_labels$gene_label
heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z <- heat_mat_z[complete.cases(heat_mat_z), , drop = FALSE]

ann_df <- metadata_cd %>% dplyr::filter(sample_id %in% colnames(heat_mat_z)) %>% dplyr::select(sample_id, Stage)
cluster_ids <- function(ids) {
  if (length(ids) <= 1) return(ids)
  hclust(dist(t(heat_mat_z[, ids, drop = FALSE])))$labels[hclust(dist(t(heat_mat_z[, ids, drop = FALSE])))$order]
}
row_order <- hclust(dist(heat_mat_z))$labels[hclust(dist(heat_mat_z))$order]
sample_order <- ann_df %>% dplyr::mutate(Stage = factor(Stage, levels = c("Stage C", "Stage D"))) %>% dplyr::arrange(Stage) %>% dplyr::group_by(Stage) %>%
  dplyr::summarise(sample_id = list(cluster_ids(sample_id)), .groups = "drop") %>% dplyr::pull(sample_id) %>% unlist()
heat_df <- heat_mat_z[row_order, sample_order, drop = FALSE] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% 
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "Expression") %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(row_order)), sample_id = factor(sample_id, levels = sample_order))
stage_df <- ann_df %>% dplyr::mutate(sample_id = factor(sample_id, levels = sample_order), y = "SCAI")

p_fig3c <- ggplot() + geom_tile(data = stage_df, aes(sample_id, y, fill = Stage), height = 0.8) + scale_fill_manual(values = stage_colors, name = "SCAI") +
  ggnewscale::new_scale_fill() + geom_tile(data = heat_df, aes(sample_id, gene, fill = Expression)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Expression") +
  labs(title = "C. Top differentially expressed genes", x = NULL, y = NULL) + theme_pub() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(), 
        legend.position = "right", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9))
p_fig3c
ggsave(filename = file.path(output_dir, "Figure3C_top_DEG_heatmap_grouped_by_stage.png"), plot = p_fig3c, width = panel_width, height = 6.5, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 3D: Pathway enrichment
# ============================================================
cap_first <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}
pathway <- readxl::read_excel(pathway_file, sheet = "pathway_enrichment_all")
head(pathway)
pathway_plot_all <- pathway %>% dplyr::mutate(Direction = factor(Direction, levels = c("Up", "Down")), sig_status = dplyr::case_when(
  Sig_FDR_0.05 == TRUE & Direction == "Up" ~ "Up_sig", Sig_FDR_0.05 != TRUE & Direction == "Up" ~ "Up_ns",
  Sig_FDR_0.05 == TRUE & Direction == "Down" ~ "Down_sig", Sig_FDR_0.05 != TRUE & Direction == "Down" ~ "Down_ns"),
  pathway_label = stringr::str_wrap(cap_first(Description), 36)) %>% dplyr::group_by(Database, Direction) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = 5) %>% 
  dplyr::ungroup() %>% dplyr::mutate(pathway_label = forcats::fct_reorder(pathway_label, logFDR))

# pathway_plot_all <- pathway %>% dplyr::filter(Sig_FDR_0.05 == TRUE) %>% dplyr::mutate(Direction = factor(Direction, levels = c("Up", "Down")), sig_status = dplyr::case_when(
#   Sig_FDR_0.05 == TRUE & Direction == "Up" ~ "Up_sig", Sig_FDR_0.05 != TRUE & Direction == "Up" ~ "Up_ns",
#   Sig_FDR_0.05 == TRUE & Direction == "Down" ~ "Down_sig", Sig_FDR_0.05 != TRUE & Direction == "Down" ~ "Down_ns"),
#   pathway_label = stringr::str_wrap(cap_first(Description), 36)) %>% dplyr::group_by(Database, Direction) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% 
#   dplyr::slice_head(n = 5) %>% dplyr::ungroup() %>% dplyr::mutate(pathway_label = forcats::fct_reorder(pathway_label, logFDR))

p_fig3d_combined <- ggplot(pathway_plot_all, aes(x = logFDR, y = pathway_label, fill = sig_status)) + geom_col(width = 0.7) +
  facet_grid(Database ~ Direction, scales = "free_y", space = "free_y") + scale_fill_manual(values = pathway_colors) +
  labs(title = "D. Pathway enrichment (ORA)", x = "-log10(FDR)", y = NULL, fill = NULL) + theme_pub() +
  theme(axis.text.y = element_text(size = 8), strip.text = element_text(face = "bold", size = 10), legend.position = "top", legend.justification = "left",
        panel.spacing = unit(0.15, "lines"), panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4))
p_fig3d_combined
ggsave(filename = file.path(output_dir, "Figure3D_pathway_top5_up_down_by_database.png"), plot = p_fig3d_combined, width = panel_width, height = 8.5, dpi = dpi_out, bg = "white")
# ggsave(filename = file.path(output_dir, "Figure3D_pathway_top5_up_down_by_database_testI.png"), plot = p_fig3d_combined, width = panel_width, height = 5, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 3E: Univariate analysis
# ============================================================
library(broom)
dge_sig <- read.csv(dge_sig_file) %>% mutate(gene_label = ifelse(is.na(gene_name) | gene_name == "", gene, gene_name), rank_score = abs(log2FoldChange) * -log10(padj)) %>%
  arrange(desc(rank_score))
de_gene_ids <- dge_sig$gene
expr_long <- expr_cd[rownames(expr_cd) %in% de_gene_ids, , drop = FALSE] %>% as.data.frame() %>% rownames_to_column("gene") %>%
  left_join(dge_sig %>% dplyr::select(gene, gene_label), by = "gene") %>% pivot_longer(-c(gene, gene_label), names_to = "sample_id", values_to = "expression") %>%
  left_join(metadata_cd %>% dplyr::select(sample_id, Stage), by = "sample_id") %>% group_by(gene_label) %>%
  mutate(expression_z = as.numeric(scale(expression)), Stage_binary = ifelse(Stage == "Stage D", 1, 0)) %>% ungroup()

uni_results <- expr_long %>% group_by(gene_label) %>% group_modify(~{
  fit <- glm(Stage_binary ~ expression_z, data = .x, family = binomial())
  tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  }) %>% ungroup() %>% filter(term == "expression_z") %>% transmute(gene_label, OR = estimate, CI_lower = conf.low, CI_upper = conf.high, p_value = p.value, FDR = p.adjust(p_value, method = "BH"),
                                                                    direction = ifelse(OR > 1, "Higher in Stage D", "Higher in Stage C"), OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, CI_lower, CI_upper)) %>%
  arrange(FDR)
uni_results
# write.csv(uni_results, file.path("cs_project/bulk_data/results/univariate_DE_gene_results.csv"), row.names = FALSE)

# fig3e_df <- uni_results %>% filter(OR > 1) %>% arrange(desc(OR)) %>% slice_head(n = 6) %>% mutate(gene_label = factor(gene_label, levels = rev(gene_label)), label_text = sprintf("OR %.2f (p=%.3f)", OR, p_value))
# p_fig3e <- ggplot(fig3e_df, aes(x = gene_label, y = OR, ymin = CI_lower, ymax = CI_upper)) + geom_pointrange(color = stage_colors["Stage D"], size = 0.7) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + geom_text(aes(label = label_text, y = CI_upper), hjust = -0.05, size = 3) + coord_flip(clip = "off") +
#   theme_pub() + theme(legend.position = "none", plot.margin = margin(5.5, 80, 5.5, 5.5)) + labs(title = "E. Stage D associated DE genes", x = "", y = "Odds Ratio (95% CI)")
# p_fig3e

fig3e_df <- uni_results %>% filter(OR > 1) %>% arrange(desc(OR)) %>% slice_head(n = 6) %>% 
  mutate(gene_label = factor(gene_label, levels = rev(gene_label)), label_text = sprintf("OR %.2f %s", OR, sig_star(p_value)), is_sig = p_value < 0.05, label_face = ifelse(is_sig, "bold", "plain"))
p_fig3e <- ggplot(fig3e_df, aes(x = gene_label, y = OR, ymin = CI_lower, ymax = CI_upper)) + geom_pointrange(color = stage_colors["Stage D"], size = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +  geom_text(aes(label = label_text, y = CI_upper, fontface = label_face), hjust = -0.05, size = 3) +  coord_flip(clip = "off") +
  theme_pub() +  theme(legend.position = "none", plot.margin = margin(5.5, 80, 5.5, 5.5), axis.text.y = element_text(face = "bold")) + labs(title = "E. Stage D associated DE genes", x = "", y = "Odds Ratio (95% CI)")
p_fig3e
ggsave(filename = file.path(output_dir, "Figure3E_univariate_analysis.png"), plot = p_fig3e, width = panel_width, height = 3, dpi = dpi_out, bg = "white")


# ============================================================
# Figure 3F: CIBERSORTx inferred immune composition
# ============================================================
ciber <- read.csv(ciber_file, check.names = FALSE)
immune_long <- ciber %>% dplyr::rename(sample_id = Mixture) %>% dplyr::inner_join(metadata_cd %>% dplyr::select(sample_id, Stage), by = "sample_id") %>%
  dplyr::select(sample_id, Stage, `B cells naive`:`Neutrophils`) %>% tidyr::pivot_longer(cols = -c(sample_id, Stage), names_to = "Cell_type", values_to = "Fraction")
immune_tests <- immune_long %>% dplyr::group_by(Cell_type) %>% dplyr::summarise(p_value = tryCatch(wilcox.test(Fraction ~ Stage)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  dplyr::mutate( sig = sig_star(p_value))
top_cells <- immune_long %>% dplyr::group_by(Cell_type) %>% dplyr::summarise(mean_fraction = mean(Fraction, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_fraction)) %>% dplyr::slice_head(n = 8) %>% dplyr::pull(Cell_type)
immune_plot <- immune_long %>% dplyr::filter(Cell_type %in% top_cells) %>% dplyr::left_join(immune_tests, by = "Cell_type") %>% 
  dplyr::mutate(Cell_type_plot = stringr::str_wrap(Cell_type, 22), Cell_type_plot = factor(Cell_type_plot, levels = stringr::str_wrap(top_cells, 22)))
immune_sig_pos <- immune_plot %>% dplyr::group_by(Cell_type_plot) %>% dplyr::summarise(y_pos = max(Fraction, na.rm = TRUE) * 1.10, sig = unique(sig)[1], .groups = "drop") %>% dplyr::filter(sig != "")

p_fig3f <- ggplot(immune_plot, aes(x = Stage, y = Fraction, fill = Stage)) + geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.25) + geom_jitter(width = 0.07, size = 1.4, alpha = 0.75) +
  geom_text(data = immune_sig_pos, aes(x = 1.5, y = y_pos, label = sig), inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  facet_wrap(~ Cell_type_plot, scales = "free_y", ncol = 4) + scale_fill_manual(values = stage_colors) + 
  labs(title = "F. Inferred immune-cell composition", x = NULL, y = "Estimated fraction", fill = NULL, caption = "*p < 0.05, **p < 0.01, ***p < 0.001") + theme_pub() +
  theme(strip.text = element_text(face = "bold", size = 9), axis.text.x = element_text(size = 8))
p_fig3f
ggsave(filename = file.path(output_dir, "Figure3F_CIBERSORTx_cell_fractions.png"), plot = p_fig3f, width = 7, height = 5.5, dpi = dpi_out, bg = "white")



# # ============================================================
# # Supplementary Figure: CS-relevant pathway landscape
# # ============================================================
# pathway <- readxl::read_excel(pathway_file, sheet = "pathway_enrichment_all")
# cap_first <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}
# cs_theme_lookup <- tibble::tribble(
#   ~Description,                                                                              ~Theme,
#   "Hemostasis",                                                                               "Coagulation & platelet",
#   "Platelet activation, signaling and aggregation",                                           "Coagulation & platelet",
#   "Platelet degranulation",                                                                   "Coagulation & platelet",
#   "Formation of Fibrin Clot (Clotting Cascade)",                                               "Coagulation & platelet",
#   "Complement and coagulation cascades",                                                       "Coagulation & platelet",
#   "platelet activation",                                                                       "Coagulation & platelet",
#   "blood coagulation",                                                                         "Coagulation & platelet",
#   "hemostasis",                                                                                "Coagulation & platelet",
#   "complement activation",                                                                     "Complement",
#   "regulation of complement activation",                                                       "Complement",
#   "complement receptor mediated signaling pathway",                                            "Complement",
#   "pattern recognition receptor signaling pathway",                                            "PRR & innate immune",
#   "toll-like receptor signaling pathway",                                                      "PRR & innate immune",
#   "innate immune response-activating signaling pathway",                                       "PRR & innate immune",
#   "regulation of innate immune response",                                                      "PRR & innate immune",
#   "Neutrophil extracellular trap formation",                                                   "Neutrophil & NETs",
#   "neutrophil mediated cytotoxicity",                                                          "Neutrophil & NETs",
#   "neutrophil chemotaxis",                                                                     "Neutrophil & NETs",
#   "humoral immune response",                                                                   "Immunoglobulin & humoral",
#   "immunoglobulin production",                                                                 "Immunoglobulin & humoral",
#   "immunoglobulin mediated immune response",                                                   "Immunoglobulin & humoral",
#   "interferon-mediated signaling pathway",                                                     "Interferon",
#   "response to type II interferon",                                                            "Interferon",
#   "mitochondrion organization",                                                                "Mitochondrial & apoptosis",
#   "regulation of mitochondrial membrane permeability",                                         "Mitochondrial & apoptosis",
#   "apoptotic mitochondrial changes",                                                            "Mitochondrial & apoptosis",
#   "release of cytochrome c from mitochondria",                                                 "Mitochondrial & apoptosis",
#   "Cell surface interactions at the vascular wall",                                            "Vascular, endothelial & ECM",
#   "ECM-receptor interaction",                                                                  "Vascular, endothelial & ECM",
#   "leukocyte adhesion to vascular endothelial cell",                                           "Vascular, endothelial & ECM",
#   "vascular associated smooth muscle cell proliferation",                                      "Vascular, endothelial & ECM",
#   "cellular response to vascular endothelial growth factor stimulus",                          "Vascular, endothelial & ECM",
#   "wound healing",                                                                             "Wound healing",
#   "response to wounding",                                                                      "Wound healing",
#   "regulation of wound healing",                                                               "Wound healing",
#   "RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function",      "Megakaryocyte & TF programs",
#   "megakaryocyte differentiation",                                                             "Megakaryocyte & TF programs"
# )
# pathway_joined <- pathway %>% dplyr::inner_join(cs_theme_lookup, by = "Description") %>% dplyr::mutate(pathway_label = stringr::str_wrap(cap_first(Description), 34)) %>%
#   dplyr::distinct(pathway_label, Direction, .keep_all = TRUE) %>% dplyr::mutate(Direction = factor(Direction, levels = c("Down", "Up")))
# pathway_points <- pathway_joined %>% dplyr::select(pathway_label, Theme, Database, Direction, logFDR, Sig_FDR_0.05)
# pathway_seg <- pathway_points %>% tidyr::pivot_wider(id_cols = c(pathway_label, Theme, Database), names_from = Direction, values_from = logFDR) %>% dplyr::filter(!is.na(Up), !is.na(Down))
# order_key <- pathway_points %>% dplyr::group_by(pathway_label, Theme) %>% dplyr::summarise(order_val = max(logFDR, na.rm = TRUE), .groups = "drop")
# pathway_points <- pathway_points %>% dplyr::left_join(order_key, by = c("pathway_label", "Theme"))
# pathway_seg    <- pathway_seg    %>% dplyr::left_join(order_key, by = c("pathway_label", "Theme"))
# lvl_order <- order_key %>% dplyr::arrange(Theme, dplyr::desc(order_val)) %>% dplyr::pull(pathway_label) %>% unique()
# pathway_points$pathway_label <- factor(pathway_points$pathway_label, levels = lvl_order)
# pathway_seg$pathway_label    <- factor(pathway_seg$pathway_label, levels = lvl_order)
# sig_thresh <- -log10(0.05)
# db_colors <- c("GO_BP" = "#1B9E77", "KEGG" = "#D95F02", "Reactome" = "#7570B3")
# 
# write.csv(pathway_points, file = 'cs_project/bulk_data/results/pathway_supp_figure.csv')
# 
# p_dumbbell <- ggplot() + geom_segment(data = pathway_seg, aes(x = Down, xend = Up, y = pathway_label, yend = pathway_label, color = Database), linewidth = 0.6, alpha = 0.5) +
#   geom_vline(xintercept = sig_thresh, linetype = "dashed", color = "grey40", linewidth = 0.4) + 
#   geom_point(data = pathway_points, aes(x = logFDR, y = pathway_label, color = Database, shape = Direction, alpha = Sig_FDR_0.05), size = 3) +  facet_grid(Theme ~ ., scales = "free_y", space = "free_y") +
#   scale_color_manual(values = db_colors, name = "Database") + scale_shape_manual(values = c("Down" = 16, "Up" = 17), name = "Direction") + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
#   labs(title = "Supp. 3A. Categorized pathway enrichment", x = "-log10(FDR)", y = NULL) +
#   theme_pub_small() + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 7, hjust = 0), axis.text.y = element_text(size = 6), panel.spacing = unit(0.4, "lines"), panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4), legend.position = "top", legend.justification = "left")
# p_dumbbell
# ggsave(filename = file.path(output_dir, "SuppFigure_pathway_landscape.png"), plot = p_dumbbell, width = panel_width, height = 8, dpi = dpi_out, bg = "white")



# ============================================================
# Supplementary Figure: CS-relevant pathway landscape - II ssGSEA using TPM
# ============================================================
known_pathway_gsea_tpm_res <- "cs_project/bulk_data/results/ssgsea_known_pathways_tpm_results.rds"
gsea_res_known   <- readRDS(known_pathway_gsea_tpm_res)
stat_res_known   <- gsea_res_known$stat_res
target_pathways  <- gsea_res_known$target_pathways
cap_first <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}
ssgsea_points <- stat_res_known %>% dplyr::left_join(target_pathways %>% dplyr::select(Description, Theme, Database), by = c("pathway" = "Description")) %>%
  dplyr::mutate(pathway_label_raw = cap_first(pathway), logFDR = -log10(p_adj), diff = mean_D - mean_C, Direction = factor(dplyr::if_else(diff > 0, "Up", "Down"), levels = c("Down", "Up")), sig = p_wilcox < 0.05)
ssgsea_points <- ssgsea_points %>% 
  dplyr::mutate(pathway_label_raw = dplyr::if_else(pathway_label_raw %in% pathway_label_raw[duplicated(pathway_label_raw)], paste0(pathway_label_raw, " (", Database, ")"), pathway_label_raw),
                pathway_label_wrapped = stringr::str_replace_all(stringr::str_wrap(pathway_label_raw, 32), "\n", "<br>"), pathway_label = dplyr::if_else(sig, paste0("<b>", pathway_label_wrapped, "</b>"), pathway_label_wrapped))
order_key <- ssgsea_points %>% dplyr::group_by(pathway_label, Theme) %>% dplyr::summarise(order_val = max(logFDR, na.rm = TRUE), .groups = "drop")
lvl_order <- order_key %>% dplyr::arrange(Theme, dplyr::desc(order_val)) %>% dplyr::pull(pathway_label) %>% unique()
ssgsea_points$pathway_label <- factor(ssgsea_points$pathway_label, levels = lvl_order)
sig_thresh <- -log10(0.05)

# p_ssgsea_summary <- ggplot(ssgsea_points, aes(x = diff, y = pathway_label, color = Direction, alpha = p_adj < 0.05)) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
#   geom_point(size = 3) + geom_text(aes(x = Inf, label = paste0("p = ", signif(p_wilcox, 2)), fontface = ifelse(sig, "bold", "plain")), color = "black", hjust = 1.05, size = 2.2) +
#   facet_grid(Theme ~ ., scales = "free_y", space = "free_y") + scale_color_manual(values = c("Down" = "steelblue", "Up" = "darkred"), name = "Direction") + theme_pub_small() +
#   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") + labs(title = "Supp. 3A. Known pathway enrichment (ssGSEA)", x = "Mean ssGSEA score difference (Stage D - Stage C)", y = NULL) +
#   theme(strip.text.y = element_text(angle = 0, face = "bold", size = 7, hjust = 0), axis.text.y = ggtext::element_markdown(size = 6), panel.spacing = unit(0.4, "lines"), panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4), legend.position = "top", legend.justification = "left")
# p_ssgsea_summary
# ggsave(filename = file.path(output_dir, "SuppFigure_known_pathways_ssGSEA_summary_dumbell.png"), plot = p_ssgsea_summary, width = panel_width, height = 8, dpi = dpi_out, bg = "white")

library(ggh4x)
library(ggtext)
scores_long_known <- gsea_res_known$scores_long %>% dplyr::mutate(stage_label = dplyr::recode(scai_stage_bulk, "C" = "Stage C", "D" = "Stage D")) %>% dplyr::left_join(ssgsea_points %>% dplyr::select(pathway, pathway_label), by = "pathway")
# pval_labels <- ssgsea_points %>% dplyr::select(pathway, pathway_label, Theme, sig) %>% dplyr::left_join(stat_res_known %>% dplyr::select(pathway, p_wilcox), by = "pathway") %>% dplyr::mutate(label = paste0("p = ", signif(p_wilcox, 2)))
pval_labels <- ssgsea_points %>% dplyr::select(pathway, pathway_label, Theme, sig) %>% dplyr::left_join(stat_res_known %>% dplyr::select(pathway, p_wilcox), by = "pathway") %>% dplyr::mutate(star = sig_star(p_wilcox), label = ifelse(star == "", "ns", star))
p_theme_box <- ggplot(scores_long_known, aes(x = es, y = pathway_label, fill = stage_label)) + geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(0.75), width = 0.6) + theme_pub_small() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 0.8, alpha = 0.5) + 
  geom_text(data = pval_labels, aes(x = Inf, y = pathway_label, label = label, fontface = ifelse(sig, "bold", "plain")), inherit.aes = FALSE, hjust = 2, size = 2.2) +
  scale_fill_manual(values = stage_colors) + ggh4x::facet_wrap2(vars(Theme), scales = "free_y", ncol = 2, axes = "all") + 
  labs(y = NULL, x = "ssGSEA enrichment score", fill = NULL, title = "Supp. 3A. Known pathway enrichment (ssGSEA)") + 
  theme(axis.text.y = ggtext::element_markdown(size = 5.5), axis.text.x = element_text(size = 5.5), legend.position = "top", legend.justification = "left", strip.text = element_text(size = 6.5)) +
  theme(panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6), panel.background = element_rect(fill = NA))
p_theme_box
ggsave(file.path(output_dir, "SuppFigure_known_pathways_ssGSEA_summary_boxplots.png"), p_theme_box, width = panel_width, height = 8, dpi = dpi_out, bg = "white")


# # ============================================================
# # Supplementary Figure: Custom pathway enrichment (GSEA)
# # ============================================================
# gsea_res <- readRDS(custom_pathawy_file)
# gsea_df <- as.data.frame(gsea_res) %>% dplyr::mutate(Sig_FDR_0.05 = p.adjust < 0.05, Direction = dplyr::if_else(NES > 0, "Higher in Stage D", "Higher in Stage C"), logFDR = -log10(p.adjust)) %>% dplyr::arrange(p.adjust)
# all_sets <- gsea_df$ID
# 
# gsea_plot_list <- lapply(all_sets, function(set_id) {
#   stats <- gsea_df %>% dplyr::filter(ID == set_id)
#   set_color <- direction_colors[[stats$Direction]]
#   title_wrapped <- stringr::str_replace(set_id, "\\s*(\\([A-Z]{2,3}:\\d+\\))$", "\n\\1")
#   p <- enrichplot::gseaplot2(gsea_res, geneSetID = set_id, title = title_wrapped, pvalue_table = FALSE, subplots = 1:2, color = set_color)
#   p[[1]] <- p[[1]] + labs(y = "Enrichment score") + annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 2.25, color = "black", label = sprintf("p = %.3g\np.adj = %.3g", stats$pvalue, stats$p.adjust)) +
#     theme_pub_small() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
#   p[[2]] <- p[[2]] + theme_pub_small() + theme(plot.title = element_blank(), plot.subtitle = element_blank())
#   ggplotify::as.ggplot(p)
# })
# p_gsea_grid <- patchwork::wrap_plots(gsea_plot_list, ncol = 2) + patchwork::plot_annotation(title = "Supp. 4A. Selected Pathway (GSEA)", theme = theme_pub_small())
# p_gsea_grid
# ggsave(filename = file.path(output_dir, "SuppFigure_custom_pathway_GSEA.png"), plot = p_gsea_grid, width = panel_width, height = 8, dpi = dpi_out, bg = "white")

# ============================================================
# Supplementary Figure: Custom pathway enrichment (GSEA) - II using TPM
# ============================================================
library(ggpubr)
gsea_res_tpm <- readRDS(custom_pathawy_gsea_tpm_res)
stat_res    <- gsea_res_tpm$stat_res
scores_long <- gsea_res_tpm$scores_long %>% dplyr::mutate(stage_label = dplyr::recode(scai_stage_bulk, "C" = "Stage C", "D" = "Stage D"), pathway = factor(pathway, levels = stat_res$pathway))
p_labels <- stat_res %>% dplyr::mutate(pw_name = stringr::str_remove(pathway, "\\s*\\([A-Z]{2,3}:\\d+\\)\\s*$"), pw_id   = stringr::str_extract(pathway, "\\([A-Z]{2,3}:\\d+\\)"),
                                       label   = paste0(pw_name, ifelse(is.na(pw_id), "", paste0("\n", pw_id)))) %>% dplyr::select(pathway, label) %>% tibble::deframe()
p_all <- ggplot(scores_long, aes(x = stage_label, y = es, fill = stage_label)) + geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.55) + geom_jitter(width = 0.12, size = 1.3, alpha = 0.6, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2.5, comparisons = list(c("Stage C", "Stage D")), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))) + 
  scale_fill_manual(values = stage_colors) + facet_wrap(~ pathway, scales = "free_y", ncol = 3, labeller = labeller(pathway = p_labels)) +
  labs(x = NULL, y = "Enrichment score") + theme_pub_small() + labs(title = "Supp. 4A. Selected Pathway Enrichment (ssGSEA)") + 
  theme(plot.title = element_text(face = "bold", size = 7), strip.text = element_text(size = 6.5, face = "bold"), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7))
p_all
ggsave(file.path(output_dir, "SuppFigure_custom_pathway_GSEA_TPM.png"), p_all, width = panel_width + 2, height = 5, dpi = dpi_out, bg = "white")

# ============================================================
# Supplementary Figure: Co_expression (WGCNA)
# ============================================================
wgcna_results <- readRDS(co_exp_file)
MEs <- wgcna_results$net$MEs
n_genes_per_module <- 500

hclust_to_segments <- function(hc) {
  n <- nrow(hc$merge); x_pos <- numeric(2 * n + 1); x_pos[1:(n + 1)] <- match(seq_len(n + 1), hc$order); segs <- vector("list", n)
  for (i in seq_len(n)) {
    left <- hc$merge[i, 1]; right <- hc$merge[i, 2]
    xl <- if (left < 0) x_pos[-left] else x_pos[n + 1 + left]; xr <- if (right < 0) x_pos[-right] else x_pos[n + 1 + right]
    yl <- if (left < 0) 0 else hc$height[left]; yr <- if (right < 0) 0 else hc$height[right]; x_pos[n + 1 + i] <- (xl + xr) / 2
    segs[[i]] <- data.frame(x = c(xl, xl, xr), xend = c(xl, xr, xr), y = c(yl, hc$height[i], hc$height[i]), yend = c(hc$height[i], hc$height[i], yr))
  }
  dplyr::bind_rows(segs)
}

# Dendrogram
kIM <- WGCNA::signedKME(t(wgcna_results$expr_var), MEs)
hub_genes <- purrr::map_dfr(unique(wgcna_results$net$colors), function(m) {
  genes_m <- names(wgcna_results$net$colors)[wgcna_results$net$colors == m]
  kME_col <- paste0("kME", m)
  if (length(genes_m) == 0 || !kME_col %in% colnames(kIM)) return(NULL)
  kIM[genes_m, kME_col, drop = FALSE] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% dplyr::arrange(dplyr::desc(.data[[kME_col]])) %>% 
    dplyr::slice_head(n = min(n_genes_per_module, length(genes_m))) %>% dplyr::mutate(module = m)
})
hub_gene_ids <- hub_genes$gene
hc_hub <- hclust(as.dist(1 - cor(t(wgcna_results$expr_var[hub_gene_ids, ]))), method = "average")
seg_data_hub <- hclust_to_segments(hc_hub)
color_strip_hub <- data.frame(x = seq_along(hc_hub$order), color = labels2colors(wgcna_results$net$colors[hub_gene_ids])[hc_hub$order])

p_dendro_hub <- ggplot() + geom_segment(data = seg_data_hub, aes(x = x, y = y, xend = xend, yend = yend), color = "grey20", linewidth = 0.2) +
  geom_tile(data = color_strip_hub, aes(x = x, y = -0.03, fill = color), height = 0.025) + scale_fill_identity() +
  labs(title = "Supp. 5A. Module hub-gene co-expression dendrogram", subtitle = paste0(n_genes_per_module, "/module, color = module"), x = "Modules", y = "Height") +
  theme_pub_small() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
p_dendro_hub
ggsave(file.path(output_dir, "SuppFigure_WGCNA_hubDendrogram.png"), p_dendro_hub, width = panel_width, height = 3, dpi = dpi_out, bg = "white")

# Module relationship network
module_sizes <- as.data.frame(table(wgcna_results$net$colors))
colnames(module_sizes) <- c("label", "size")
module_sizes$module <- paste0("ME", module_sizes$label)

plot_df <- wgcna_results$stage_results %>% dplyr::left_join(module_sizes %>% dplyr::select(module, size), by = "module") %>% 
  dplyr::mutate(direction = dplyr::case_when(cor > 0 & pval < 0.05 ~ "Up_sig", cor > 0 & pval >= 0.05 ~ "Up_ns", cor < 0 & pval < 0.05 ~ "Down_sig", cor < 0 & pval >= 0.05 ~ "Down_ns"), sig = sig_star(pval), module = forcats::fct_reorder(module, cor))
meCor <- cor(MEs)
edges <- as.data.frame(as.table(meCor)) %>% dplyr::rename(from = Var1, to = Var2, weight = Freq) %>% dplyr::filter(from != to, abs(weight) > 0.3) %>%
  dplyr::mutate(pair = purrr::map2_chr(as.character(from), as.character(to), ~paste(sort(c(.x, .y)), collapse = "_"))) %>% dplyr::distinct(pair, .keep_all = TRUE) %>% 
  dplyr::mutate(abs_weight = abs(weight)) %>% dplyr::select(from, to, weight, abs_weight)
g_modules <- graph_from_data_frame(edges, vertices = plot_df %>% dplyr::select(module, size, direction, cor), directed = FALSE)

set.seed(42)
p_moduleGraph <- ggraph(g_modules, layout = "fr", weights = abs_weight) + geom_edge_link(aes(width = abs(weight), color = weight > 0), alpha = 0.5) +
  scale_edge_color_manual(values = c("TRUE" = stage_colors[["Stage D"]], "FALSE" = stage_colors[["Stage C"]]), name = "Eigengene correlation", labels = c("Negative", "Positive")) +
  scale_edge_width_continuous(range = c(0.3, 1.8), guide = "none") + geom_node_point(aes(size = size, color = direction), alpha = 0.9) + geom_node_text(aes(label = name), size = 2.5, fontface = "bold", color = "white") +
  scale_color_manual(values = pathway_colors, drop = FALSE, name = "Module-trait association", labels = c("Up_sig" = "Nominal (Higher in D)", "Up_ns" = "Not sig. (Higher in D)", "Down_sig" = "Nominal (Higher in C)", "Down_ns" = "Not sig. (Higher in C)")) +
  scale_size_continuous(range = c(10, 16), name = "Module size (genes)") + labs(title = "Supp. 5B. Module relationship network") + theme_pub() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), 
        legend.position = "right", legend.box = "vertical", legend.title = element_text(size = 9, face = "plain"), legend.text = element_text(size = 8), plot.title = element_text(face = "bold", size = 13, hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1), size = guide_legend(order = 2), edge_color = guide_legend(order = 3))
p_moduleGraph
ggsave(file.path(output_dir, "SuppFigure_WGCNA_moduleGraph.png"), p_moduleGraph, width = panel_width, height = 6, dpi = dpi_out, bg = "white")
