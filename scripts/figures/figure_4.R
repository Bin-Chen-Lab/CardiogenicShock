setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Figure 4: Cytokine profiling in SCAI Stage D vs C shock
# ============================================================
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggnewscale)


# ---- Files ----
cytokine_file <- "cs_project/cytokine_data/CScytokineData.xlsx"
univ_file <- "cs_project/cytokine_data/results/cytokine_univariate_results_C_vs_D.csv"

output_dir <- "cs_project/figures/Figure4_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Settings ----
dpi_out <- 600
panel_width <- 6.5
p_sig <- 0.05
stage_colors <- c("Stage C" = "steelblue", "Stage D" = "darkred")
direction_colors <- c("Higher in Stage C" = "steelblue", "Higher in Stage D" = "darkred", "Not significant" = "gray75")
theme_pub <- function(base_size = 12) {theme_classic(base_size = base_size) + 
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9, hjust = 0),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10), axis.title = element_text(size = 11), 
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
          legend.position = "none", legend.justification = "left", legend.title = element_blank())
}

sig_star <- function(p) {dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}
clean_cytokine <- function(x) {x %>%
    stringr::str_replace_all("\\s*\\(pg/ml\\)", "") %>% stringr::str_replace_all("-", "") %>% stringr::str_replace_all("_", "") %>% 
    stringr::str_replace_all("\\s+", "") #%>% stringr::str_to_upper()
}
pretty_cytokine <- function(x) {x %>% 
    stringr::str_replace_all("IL6", "IL-6") %>% stringr::str_replace_all("IL8", "IL-8") %>% stringr::str_replace_all("IL10", "IL-10") %>%
    stringr::str_replace_all("IL18", "IL-18") %>% stringr::str_replace_all("IL1B", "IL-1β") %>% stringr::str_replace_all("IL1A", "IL-1α") %>%
    stringr::str_replace_all("TNFA", "TNFα") %>% stringr::str_replace_all("IFNY", "IFNγ")
}

# ============================================================
# Load data
# ============================================================
cyto_raw <- readxl::read_excel(cytokine_file, sheet = 1)
univ <- read.csv(univ_file)
cytokine_cols <- grep("\\(pg/ml\\)", names(cyto_raw), value = TRUE)
cyto_long <- cyto_raw %>% dplyr::filter(`SCAI Stage` %in% c("C", "D")) %>% dplyr::mutate(Stage = factor(`SCAI Stage`, levels = c("C", "D"), labels = c("Stage C", "Stage D"))) %>%
  dplyr::select(`SR specimen ID`, Stage, dplyr::all_of(cytokine_cols)) %>% tidyr::pivot_longer(cols = dplyr::all_of(cytokine_cols), names_to = "cytokine_col", values_to = "value_raw") %>%
  dplyr::mutate(cytokine_key = clean_cytokine(cytokine_col), value_log = log2(value_raw + 1))

univ_plot <- univ %>% dplyr::mutate(cytokine_key = clean_cytokine(Cytokine), 
                                    cytokine_label = pretty_cytokine(Cytokine), log2FC = log2(fold_change_D_vs_C), neg_log10_p = -log10(p_value), p_star = sig_star(p_value), 
                                    direction_label = dplyr::case_when(p_value < p_sig & log2FC > 0 ~ "Higher in Stage D", p_value < p_sig & log2FC < 0 ~ "Higher in Stage C", TRUE ~ "Not significant"))


# ============================================================
# Figure 4A: Cytokine heatmap
# ============================================================
heat_cytokines <- univ_plot %>% dplyr::arrange(p_value) %>% dplyr::pull(cytokine_key)
heat_dat <- cyto_long %>% dplyr::filter(cytokine_key %in% heat_cytokines) %>% dplyr::left_join(univ_plot %>% dplyr::select(cytokine_key, cytokine_label, p_value), by = "cytokine_key") %>%
  dplyr::select(sample_id = `SR specimen ID`, Stage, cytokine_label, value_log)
heat_mat <- heat_dat %>% dplyr::select(sample_id, cytokine_label, value_log) %>% tidyr::pivot_wider(names_from = sample_id, values_from = value_log) %>% as.data.frame()
rownames(heat_mat) <- heat_mat$cytokine_label
heat_mat$cytokine_label <- NULL
# # heat_mat_z <- heat_mat
# heat_mat_z <- t(scale(t(as.matrix(heat_mat))))
heat_mat_z <- heat_mat_z[complete.cases(heat_mat_z), , drop = FALSE]
ann_df <- heat_dat %>% dplyr::distinct(sample_id, Stage) %>% dplyr::filter(sample_id %in% colnames(heat_mat_z))
cluster_ids <- function(ids) {
  if (length(ids) <= 1) return(ids)
  hc <- hclust(dist(t(heat_mat_z[, ids, drop = FALSE])))
  hc$labels[hc$order]
}
row_hc <- hclust(dist(heat_mat_z))
row_order <- row_hc$labels[row_hc$order]
sample_order <- ann_df %>% dplyr::mutate(Stage = factor(Stage, levels = c("Stage C", "Stage D"))) %>% dplyr::arrange(Stage) %>% dplyr::group_by(Stage) %>%
  dplyr::summarise(sample_id = list(cluster_ids(sample_id)), .groups = "drop") %>% dplyr::pull(sample_id) %>% unlist()
heat_df <- heat_mat_z[row_order, sample_order, drop = FALSE] %>% as.data.frame() %>% tibble::rownames_to_column("cytokine") %>%
  tidyr::pivot_longer(-cytokine, names_to = "sample_id", values_to = "Level") %>%
  dplyr::mutate(cytokine = factor(cytokine, levels = rev(row_order)), sample_id = factor(sample_id, levels = sample_order))
stage_df <- ann_df %>% dplyr::mutate(sample_id = factor(sample_id, levels = sample_order), y = "SCAI")

p_fig4a <- ggplot() + geom_tile(data = stage_df, aes(sample_id, y, fill = Stage), height = 0.8) + scale_fill_manual(values = stage_colors, name = "SCAI") +
  ggnewscale::new_scale_fill() + geom_tile(data = heat_df, aes(sample_id, cytokine, fill = Level)) + 
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Level") +
  labs( title = "A. Cytokine inflammatory phenotype", x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(), legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9))
p_fig4a
ggsave(filename = file.path(output_dir, "Figure4A_cytokine_heatmap_grouped_by_stage_testII.png"), plot = p_fig4a, width = panel_width, height = 6.3, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 4B: Cytokine volcano plot
# ============================================================
top_labels <- univ_plot %>% dplyr::arrange(p_value)
p_fig4b <- ggplot(univ_plot, aes(x = log2FC, y = neg_log10_p, color = direction_label)) + geom_point(size = 2.8, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray45") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "gray45") +
  ggrepel::geom_text_repel(data = top_labels, aes(label = cytokine_label), size = 4.5, box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = direction_colors) +labs(title = "B. Differential cytokine level", x = "log2 fold-change: Stage D vs Stage C", y = "-log10(p value)") + theme_pub()
p_fig4b
ggsave(filename = file.path(output_dir, "Figure4B_cytokine_volcano.png"), plot = p_fig4b, width = panel_width, height = 5.2, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 4C: Violin + boxplot + jitter for top cytokines
# ============================================================
top_cytokines <- univ_plot %>% dplyr::arrange(p_value) %>% dplyr::slice_head(n = 12) %>% dplyr::pull(cytokine_key)
violin_dat <- cyto_long %>% dplyr::filter(cytokine_key %in% top_cytokines) %>% dplyr::left_join(univ_plot %>% dplyr::select(cytokine_key, cytokine_label, p_value, p_star), by = "cytokine_key") %>%
  dplyr::mutate(cytokine_label = factor(cytokine_label, levels = univ_plot %>% dplyr::filter(cytokine_key %in% top_cytokines) %>% dplyr::arrange(p_value) %>% dplyr::pull(cytokine_label)))
sig_pos <- violin_dat %>% dplyr::group_by(cytokine_label) %>% dplyr::summarise(y_pos = max(value_log, na.rm = TRUE) * 1.08, p_star = unique(p_star)[1], .groups = "drop") %>% dplyr::filter(p_star != "")
p_fig4c <- ggplot(violin_dat, aes(x = Stage, y = value_log, fill = Stage)) + geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) + geom_jitter(width = 0.08, size = 1.6, alpha = 0.8) +
  geom_text(data = sig_pos, aes(x = 1.5, y = y_pos, label = p_star), inherit.aes = FALSE, size = 4.5, fontface = "bold") + 
  facet_wrap(~ cytokine_label, scales = "free_y", ncol = 6) + scale_fill_manual(values = stage_colors) +
  labs(title = "C. Distribution of cytokines", x = NULL, y = "log2(pg/mL + 1)", fill = NULL, caption = "*p < 0.05, **p < 0.01, ***p < 0.001") + theme_pub() +
  theme(strip.text = element_text(face = "bold", size = 11), axis.text.x = element_text(size = 9))
p_fig4c
ggsave(filename = file.path(output_dir, "Figure4C_cytokines_violin.png"), plot = p_fig4c, width = panel_width*2, height = 6, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 4D: Top cytokines by sex and SCAI stage
# ============================================================
cytokines_to_plot <- c("ST2", "IL6", "VEGFA")
fig4d_dat <- cyto_long %>% dplyr::filter(cytokine_key %in% cytokines_to_plot) %>% dplyr::left_join(cyto_raw %>% dplyr::select(`SR specimen ID`, Sex), by = "SR specimen ID") %>%
  dplyr::mutate(Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")), 
                Stage_Sex = factor(interaction(Stage, Sex, sep = " - "), levels = c("Stage C - Female", "Stage C - Male", "Stage D - Female", "Stage D - Male")), cytokine_label = pretty_cytokine(cytokine_key)) %>%
  dplyr::filter(!is.na(Sex), !is.na(value_log))
stage_sex_colors <- c("Stage C - Female" = "#6BAED6", "Stage C - Male"   = "steelblue", "Stage D - Female" = "#CB6A6A", "Stage D - Male"   = "darkred")
comparisons_stage_sex <- list(#c("Stage C - Female", "Stage C - Male"), c("Stage D - Female", "Stage D - Male"), 
                              c("Stage C - Female", "Stage D - Female"), c("Stage C - Male", "Stage D - Male"))
p_fig4d <- ggplot(fig4d_dat, aes(Stage_Sex, value_log, fill = Stage_Sex)) + geom_violin(trim = FALSE, alpha = 0.65, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) + geom_jitter(aes(color = Stage_Sex), width = 0.08, size = 1.6, alpha = 0.8, show.legend = FALSE) +
  ggpubr::stat_compare_means(comparisons = comparisons_stage_sex, method = "wilcox.test", label = "p.signif", hide.ns = F, size = 4, tip.length = 0.01) +
  facet_wrap(~ cytokine_label, scales = "free_y", nrow = 1) + scale_fill_manual(values = stage_sex_colors) + scale_color_manual(values = stage_sex_colors) +
  labs(title = "D. Selected cytokines by sex and SCAI stage", x = NULL, y = "log2(pg/mL + 1)", caption = "*p < 0.05, **p < 0.01, ***p < 0.001.") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
p_fig4d
ggsave(filename = file.path(output_dir, "Figure4D_selected_cytokines_by_sex_stage.png"), plot = p_fig4d, width = panel_width, height = 4, dpi = dpi_out, bg = "white")

# ============================================================
# Figure 3E: Univariate analysis
# ============================================================
library(broom)
cyto_keep <- cyto_long %>% group_by(cytokine_key, Stage) %>% summarise(n_nonzero = sum(!is.na(value_raw) & value_raw > 0), .groups = "drop") %>%
  pivot_wider(names_from = Stage, values_from = n_nonzero, values_fill = 0) %>% filter(`Stage C` >= 3, `Stage D` >= 3) %>% pull(cytokine_key)

cyto_logit_results <- cyto_long %>% filter(cytokine_key %in% cyto_keep) %>% mutate(Stage_binary = ifelse(Stage == "Stage D", 1, 0), value_log = as.numeric(value_log)) %>%
  filter(!is.na(Stage_binary), is.finite(value_log)) %>% group_by(cytokine_key) %>% filter(n_distinct(Stage_binary) == 2, n_distinct(value_log) > 1) %>%
  mutate(value_z = as.numeric(scale(value_log))) %>% filter(is.finite(value_z)) %>%
  group_modify(~{
    fit <- glm(Stage_binary ~ value_z, data = .x, family = binomial())
    tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  }) %>% ungroup() %>% filter(term == "value_z") %>% 
  transmute(cytokine_key, cytokine_label = pretty_cytokine(cytokine_key), OR = estimate, CI_lower = conf.low, CI_upper = conf.high, 
            p_value = p.value, FDR = p.adjust(p_value, method = "BH"), 
            direction = ifelse(OR > 1, "Stage D", "Stage C"), OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, CI_lower, CI_upper)) %>% arrange(p_value)
# write.csv(cyto_logit_results, file.path("cs_project/cytokine_data/results/cytokine_univariate_logistic_results.csv"), row.names = FALSE)

# fig4e_df <- cyto_logit_results %>% filter(OR > 1) %>% arrange(desc(OR)) %>% slice_head(n = 6) %>%
#   mutate(cytokine_label = factor(cytokine_label, levels = rev(cytokine_label)), label_text = sprintf("OR %.2f (p=%.3f)", OR, p_value))
# p_fig4e <- ggplot(fig4e_df, aes(x = cytokine_label, y = OR, ymin = CI_lower, ymax = CI_upper)) + geom_pointrange(color = stage_colors["Stage D"], size = 0.7) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + geom_text(aes(label = label_text, y = CI_upper), hjust = -0.05, size = 3) + coord_flip(clip = "off") +
#   theme_pub() + theme(legend.position = "none", plot.margin = margin(5.5, 90, 5.5, 5.5)) + labs(title = "E. Stage D associated cytokines", x = "", y = "Odds Ratio (95% CI)")
# p_fig4e
# ggsave(file.path(output_dir, "Figure4E_cytokine_univariate_forestplot.png"), p_fig4e, width = 5.5, height = 3, dpi = dpi_out)

fig4e_df <- cyto_logit_results %>% filter(OR > 1) %>% arrange(desc(OR)) %>% slice_head(n = 6) %>% 
  mutate(cytokine_label = factor(cytokine_label, levels = rev(cytokine_label)), star = sig_star(p_value), star_display = ifelse(star == "", "ns", star), label_text = sprintf("OR %.2f %s", OR, star_display), is_sig = p_value < 0.05, label_face = ifelse(is_sig, "bold", "plain"))
p_fig4e <- ggplot(fig4e_df, aes(x = cytokine_label, y = OR, ymin = CI_lower, ymax = CI_upper)) + geom_pointrange(color = stage_colors["Stage D"], size = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +  geom_text(aes(label = label_text, y = CI_upper, fontface = label_face), hjust = -0.05, size = 3) + coord_flip(clip = "off") +
  theme_pub() + theme(legend.position = "none", plot.margin = margin(5.5, 90, 5.5, 5.5), axis.text.y = element_text(face = "bold")) + labs(title = "E. Stage D associated cytokines", x = "", y = "Odds Ratio (95% CI)")
p_fig4e
ggsave(file.path(output_dir, "Figure4E_cytokine_univariate_forestplot.png"), p_fig4e, width = 5.5, height = 3, dpi = dpi_out)