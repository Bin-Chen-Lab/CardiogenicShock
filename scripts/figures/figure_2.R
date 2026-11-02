setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Figure 2: Laboratory phenotype of SCAI Stage D vs C shock
# ============================================================

library(tidyverse)
library(ggplot2)
library(forcats)
library(stringr)
library(ggtext)
library(tidytext)

# ---- Files ----
input_file <- "cs_project/clinical_data/clinical_results/CS_all_variables_publication_table.csv"
clinical_file <- "cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv"
admshock_file <- "cs_project/clinical_data/clinical_results/CS_Admission_vs_Shock_paired_results.csv"

output_dir <- "cs_project/figures/Figure2_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Settings ----
dpi_out <- 600
p_sig <- 0.05
panel_width <- 6.5

stage_colors <- c("Stage C" = "steelblue", "Stage D" = "darkred")
direction_colors <- c("Higher in Stage C" = "steelblue", "Higher in Stage D" = "darkred")

clean_p <- function(x) {x <- as.character(x); x <- str_replace_all(x, "<", ""); x <- str_replace_all(x, ">", ""); suppressWarnings(as.numeric(x))}
extract_median <- function(x) {x <- as.character(x); suppressWarnings(as.numeric(str_extract(x, "^-?[0-9.]+")))}
extract_sd <- function(x) {x <- as.character(x); suppressWarnings(as.numeric(str_match(x, "±\\s*([-0-9.]+)")[, 2]))}
clean_label <- function(x) {x %>% str_replace_all("_", " ") %>% str_replace_all("\\s+", " ") %>% str_trim()}
format_p <- function(p) {case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}
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
lab_dat <- read_csv(input_file, show_col_types = FALSE) %>% rename_with(~ str_replace_all(.x, "\\s+", "_")) %>% 
  mutate(p_num = clean_p(P_Value), Phase_clean = clean_label(Phase), Category_clean = clean_label(Category), Variable_clean = clean_label(Variable),
         C_median = extract_median(`C_Median_(IQR)`), D_median = extract_median(`D_Median_(IQR)`), C_sd = extract_sd(`C_Mean_±_SD`), D_sd = extract_sd(`D_Mean_±_SD`),
         raw_diff = D_median - C_median, abs_raw_diff = abs(raw_diff), pooled_sd = sqrt((C_sd^2 + D_sd^2) / 2), std_diff = raw_diff / pooled_sd,
         sig_label = if_else(!is.na(p_num) & p_num < p_sig, "p < 0.05", "NS"), p_star = format_p(p_num),
         direction = case_when(raw_diff >= 0 ~ "Higher in Stage D", raw_diff < 0 ~ "Higher in Stage C")) %>%
  dplyr::filter(Variable_Type == "numeric", Phase_clean %in% c("Admission", "Shock onset", "Recovery", "Outcome"), !is.na(std_diff), is.finite(std_diff)) # Category_clean == "Laboratory", 

# common variables measured in 3 phases
# common_var <- lab_dat %>% dplyr::distinct(Variable_clean, Phase_clean) %>% dplyr::count(Variable_clean) %>% dplyr::filter(n == 2) %>% dplyr::pull(Variable_clean)
# lab_dat_common <- lab_dat %>% filter(Variable_clean %in% common_var, Phase_clean %in% c("Admission", "Shock onset", "Recovery")) %>%
#   mutate(Phase_clean = factor(Phase_clean, levels = c("Admission", "Shock onset", "Recovery")), 
#          Category_clean = factor(Category_clean, levels = c("Laboratory", "Pulmonary artery catheter")),
#          Variable_plot = case_when(p_num < p_sig ~ paste0("**", str_wrap(Variable_clean, 34), "**"), TRUE ~ str_wrap(Variable_clean, 34)),
#          Variable_plot = forcats::fct_reorder(Variable_plot, abs_raw_diff), sig_label = factor(sig_label, levels = c("p < 0.05", "NS")),
#          direction = factor(direction, levels = c("Higher in Stage C", "Higher in Stage D")))

# choose only 3 phases: revised one
lab_dat_common <- lab_dat %>% filter(Phase_clean %in% c("Admission", "Shock onset", "Recovery")) %>%
  mutate(Phase_clean = factor(Phase_clean, levels = c("Admission", "Shock onset", "Recovery")), 
         Category_clean = factor(Category_clean, levels = c("Laboratory", "Pulmonary artery catheter")),
         Variable_plot = case_when(p_num < p_sig ~ paste0("**", str_wrap(Variable_clean, 34), "**"), TRUE ~ str_wrap(Variable_clean, 34)),
         Variable_plot = forcats::fct_reorder(Variable_plot, abs_raw_diff), sig_label = factor(sig_label, levels = c("p < 0.05", "NS")),
         direction = factor(direction, levels = c("Higher in Stage C", "Higher in Stage D")))

# ============================================================
# Figure 2A: Variables measured
# ============================================================

# p_fig2a <- ggplot(lab_dat_common, aes(x = std_diff, y = Variable_plot, color = direction, alpha = sig_label)) + 
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "gray45") + geom_segment(aes(x = 0, xend = std_diff, yend = Variable_plot), linewidth = 0.8) + 
#   geom_point(size = 2) + geom_text(aes(label = p_star), nudge_x = 0.08, size = 4, show.legend = FALSE) + 
#   facet_grid(Category_clean ~ Phase_clean, scales = "free_y", space = "free_y") + scale_y_reordered() + 
#   scale_color_manual(values = direction_colors) + scale_alpha_manual(values = c("p < 0.05" = 1, "NS" = 0.35)) + 
#   labs(title = "A. Variables measured", x = "Median difference (Stage D - Stage C)", y = NULL, color = NULL, alpha = NULL) + theme_pub_small() + 
#   theme(axis.text.y = ggtext::element_markdown(size = 6), axis.text.x = element_text(size = 5), strip.text.x = element_text(face = "bold", size = 6), strip.text.y = element_blank(), strip.background = element_blank())
# p_fig2a
# # ggsave(filename = file.path(output_dir, "Figure2A_common_labs_admission_vs_shock_onset.png"), plot = p_fig2a, width = panel_width, height = 4, dpi = dpi_out, bg = "white")

# figure for "choose only 3 phases: revised one"
p_fig2a <- ggplot(lab_dat_common, aes(x = std_diff, y = Variable_plot, color = direction, alpha = sig_label)) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "gray45") + geom_segment(aes(x = 0, xend = std_diff, yend = Variable_plot), linewidth = 0.8) + 
  geom_point(size = 2) + geom_text(aes(label = p_star), nudge_x = 0.08, size = 4, show.legend = FALSE) + 
  facet_grid(Category_clean ~ Phase_clean, scales = "free_y", space = "free_y") + scale_y_reordered() + 
  scale_color_manual(values = direction_colors) + scale_alpha_manual(values = c("p < 0.05" = 1, "NS" = 0.35)) + 
  labs(title = "A. Variables measured", x = "Median difference (Stage D - Stage C)", y = NULL, color = NULL, alpha = NULL) + theme_pub_small() + 
  theme(axis.text.y = ggtext::element_markdown(size = 6), axis.text.x = element_text(size = 5), strip.text.x = element_text(face = "bold", size = 6), strip.text.y = element_blank(), strip.background = element_blank())
p_fig2a
ggsave(filename = file.path(output_dir, "Figure2A_common_labs_admission_vs_shock_onset_revised.png"), plot = p_fig2a, width = panel_width, height = 6.5, dpi = dpi_out, bg = "white")


# ============================================================
# Figure 2B: Admission + Shock onset violin plot
# ============================================================
clinical_data <- read_csv(clinical_file, show_col_types = FALSE)
name_key <- function(x) x %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", "_") %>% str_replace_all("^_|_$", "")
clinical_lookup <- tibble(patient_var = names(clinical_data), patient_key = name_key(names(clinical_data)))
lab_pattern <- "glucose|alkaline|nt_pro|bicarbonate|lactate|protein_total|albumin|calcium"
build_map <- function(phase_label, phase_key) {
  lab_dat %>% filter(Phase_clean == phase_label, Category_clean %in% c("Laboratory")) %>% mutate(var_key = name_key(Variable_clean)) %>% filter(str_detect(var_key, lab_pattern)) %>% arrange(p_num) %>%
    mutate(patient_var = map_chr(var_key, ~ clinical_lookup$patient_var[str_detect(clinical_lookup$patient_key, phase_key) & str_detect(clinical_lookup$patient_key, .x)][1]),
           Lab = case_when(str_detect(var_key, "glucose") ~ "Glucose\n(mg/dL)", str_detect(var_key, "alkaline") ~ "Alkaline phosphatase\n(IU/L)", str_detect(var_key, "nt_pro") ~ "NT-proBNP\n(ng/L)",
                           str_detect(var_key, "bicarbonate") ~ "Bicarbonate\n(mmol/L)", str_detect(var_key, "lactate") ~ "Lactate\n(mmol/L)", str_detect(var_key, "protein_total") ~ "Total protein\n(g/dL)",
                           str_detect(var_key, "albumin") ~ "Albumin\n(g/dL)", str_detect(var_key, "calcium") ~ "Calcium\n(mg/dL)", TRUE ~ Variable_clean)) %>%
    filter(!is.na(patient_var)) %>% distinct(Lab, .keep_all = TRUE) %>% dplyr::select(Lab, patient_var, p_star)
}
admission_map <- build_map("Admission", "admission"); shock_map <- build_map("Shock onset", "shock_onset")
build_dat <- function(map_df, phase_label) {
  clinical_data %>% filter(scai_stage %in% c("C", "D")) %>% mutate(Stage = factor(scai_stage, levels = c("C", "D"), labels = c("Stage C", "Stage D"))) %>% dplyr::select(Stage, all_of(map_df$patient_var)) %>% 
    pivot_longer(-Stage, names_to = "patient_var", values_to = "Value") %>% filter(!is.na(Value)) %>% left_join(map_df, by = "patient_var") %>% mutate(Phase = phase_label)
}
combined_dat <- bind_rows(build_dat(admission_map, "Admission"), build_dat(shock_map, "Shock onset")) %>% 
  mutate(Phase = factor(Phase, levels = c("Admission", "Shock onset")), Lab = factor(Lab, levels = unique(c(as.character(admission_map$Lab), as.character(shock_map$Lab)))))
combined_map <- bind_rows(admission_map %>% mutate(Phase = "Admission"), shock_map %>% mutate(Phase = "Shock onset")) %>% mutate(Phase = factor(Phase, levels = c("Admission", "Shock onset")))
combined_sig_pos <- combined_dat %>% group_by(Lab, Phase) %>% summarise(y_pos = max(Value, na.rm = TRUE) * 1.1, .groups = "drop") %>% left_join(combined_map %>% dplyr::select(Lab, Phase, p_star), by = c("Lab", "Phase"))

p_fig2_combined <- ggplot(combined_dat, aes(Phase, Value, fill = Stage)) + geom_violin(trim = FALSE, alpha = 0.6, color = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.25, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.07, dodge.width = 0.8), size = 1.2, alpha = 0.75) +
  geom_text(data = combined_sig_pos, aes(x = Phase, y = y_pos, label = p_star), inherit.aes = FALSE, size = 4, fontface = "bold") + facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = stage_colors) + labs(title = "B. Admission and Shock Onset Variables", x = NULL, y = NULL, fill = NULL) + theme_pub_small() + 
  theme(strip.text = element_text(face = "bold", size = 7), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
p_fig2_combined
ggsave(file.path(output_dir, "Figure2B_combined_violin.png"), p_fig2_combined, width = panel_width, height = 3, dpi = dpi_out, bg = "white")



# # ============================================================
# # Figure 2B: Admission violin plot
# # ============================================================
# clinical_data <- read_csv(clinical_file, show_col_types = FALSE)
# name_key <- function(x) x %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", "_") %>% str_replace_all("^_|_$", "")
# clinical_lookup <- tibble(patient_var = names(clinical_data), patient_key = name_key(names(clinical_data)))
# 
# admission_map <- lab_dat %>% filter(Phase_clean == "Admission", Category_clean %in% c("Laboratory")) %>% arrange(p_num) %>% slice_head(n = 6) %>%
#   mutate(var_key = name_key(Variable_clean), patient_var = map_chr(var_key, ~ clinical_lookup$patient_var[str_detect(clinical_lookup$patient_key, "admission") & str_detect(clinical_lookup$patient_key, .x)][1]), 
#          Lab = case_when(str_detect(var_key, "glucose") ~ "Glucose\n(mg/dL)", str_detect(var_key, "alkaline") ~ "Alkaline phosphatase\n(IU/L)", str_detect(var_key, "nt_pro") ~ "NT-proBNP\n(ng/L)",
#                          str_detect(var_key, "bicarbonate") ~ "Bicarbonate\n(mmol/L)", str_detect(var_key, "lactate") ~ "Lactate\n(mmol/L)", str_detect(var_key, "protein_total") ~ "Total protein\n(g/dL)", TRUE ~ Variable_clean)) %>% 
#   filter(!is.na(patient_var)) %>% dplyr::select(Lab, patient_var, p_star)
# admission_dat <- clinical_data %>% filter(scai_stage %in% c("C", "D")) %>% mutate(Stage = factor(scai_stage, levels = c("C", "D"), labels = c("Stage C", "Stage D"))) %>%
#   dplyr::select(Stage, all_of(admission_map$patient_var)) %>% pivot_longer(-Stage, names_to = "patient_var", values_to = "Value") %>% filter(!is.na(Value)) %>%
#   left_join(admission_map, by = "patient_var") %>% mutate(Lab = factor(Lab, levels = admission_map$Lab))
# admission_sig_pos <- admission_dat %>% group_by(Lab) %>% summarise(y_pos = max(Value, na.rm = TRUE) * 1.1, .groups = "drop") %>% left_join(admission_map %>% dplyr::select(Lab, p_star), by = "Lab")
# 
# p_fig2b <- ggplot(admission_dat, aes(Stage, Value, fill = Stage)) + geom_violin(trim = FALSE, alpha = 0.6, color = NA) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.25) +
#   geom_jitter(width = 0.07, size = 1.2, alpha = 0.75) + geom_text(data = admission_sig_pos, aes(x = 1.5, y = y_pos, label = p_star), inherit.aes = FALSE, size = 4, fontface = "bold") +
#   facet_wrap(~ Lab, scales = "free_y", ncol = 3) + scale_fill_manual(values = stage_colors) + labs(title = "B. Admission", x = NULL, y = NULL, fill = NULL) + theme_pub() +
#   theme(strip.text = element_text(face = "bold", size = 7), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
# p_fig2b
# # ggsave(file.path(output_dir, "Figure2B_admission_violin.png"), p_fig2b, width = panel_width, height = 4, dpi = dpi_out, bg = "white")
# 
# # ============================================================
# # Figure 2C: Shock onset violin plot
# # ============================================================
# shock_map <- lab_dat %>% filter(Phase_clean == "Shock onset", Category_clean %in% c("Laboratory")) %>% arrange(p_num) %>% slice_head(n = 4) %>%
#   mutate(var_key = name_key(Variable_clean), patient_var = map_chr(var_key, ~ clinical_lookup$patient_var[str_detect(clinical_lookup$patient_key, "shock_onset") & str_detect(clinical_lookup$patient_key, .x)][1]), 
#          Lab = case_when(str_detect(var_key, "albumin") ~ "Albumin\n(g/dL)", str_detect(var_key, "calcium") ~ "Calcium\n(mg/dL)", str_detect(var_key, "glucose") ~ "Glucose\n(mg/dL)",
#                          str_detect(var_key, "protein_total") ~ "Total protein\n(g/dL)", TRUE ~ Variable_clean)) %>% 
#   filter(!is.na(patient_var)) %>% dplyr::select(Lab, patient_var, p_star)
# shock_dat <- clinical_data %>% filter(scai_stage %in% c("C", "D")) %>% mutate(Stage = factor(scai_stage, levels = c("C", "D"), labels = c("Stage C", "Stage D"))) %>%
#   dplyr::select(Stage, all_of(shock_map$patient_var)) %>% pivot_longer(-Stage, names_to = "patient_var", values_to = "Value") %>% filter(!is.na(Value)) %>%
#   left_join(shock_map, by = "patient_var") %>% mutate(Lab = factor(Lab, levels = shock_map$Lab))
# shock_sig_pos <- shock_dat %>% group_by(Lab) %>% summarise(y_pos = max(Value, na.rm = TRUE) * 1.1, .groups = "drop") %>% left_join(shock_map %>% dplyr::select(Lab, p_star), by = "Lab")
# 
# p_fig2c <- ggplot(shock_dat, aes(Stage, Value, fill = Stage)) + geom_violin(trim = FALSE, alpha = 0.6, color = NA) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.25) +
#   geom_jitter(width = 0.07, size = 1.2, alpha = 0.75) + geom_text(data = shock_sig_pos, aes(x = 1.5, y = y_pos, label = p_star), inherit.aes = FALSE, size = 4, fontface = "bold") +
#   facet_wrap(~ Lab, scales = "free_y", ncol = 2) + scale_fill_manual(values = stage_colors) + labs(title = "C. Shock onset", x = NULL, y = NULL, fill = NULL) + theme_pub() +
#   theme(strip.text = element_text(face = "bold", size = 7), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
# p_fig2c
# # ggsave(file.path(output_dir, "Figure2C_shock_onset_violin.png"), p_fig2c, width = panel_width -1.75, height = 4, dpi = dpi_out, bg = "white")


################################
traj_vars <- lab_dat %>% filter(Phase_clean %in% c("Admission", "Shock onset", "Recovery")) %>% group_by(Variable_clean) %>% filter(n_distinct(Phase_clean) >= 2) %>%
  ungroup() %>% arrange(p_num) %>% distinct(Variable_clean)  %>% pull(Variable_clean)

traj_dat <- lab_dat %>% filter(Variable_clean %in% traj_vars, Phase_clean %in% c("Admission", "Shock onset", "Recovery")) %>%
  mutate(Phase_clean = factor(Phase_clean, levels = c("Admission", "Shock onset", "Recovery")), Variable_plot = str_wrap(Variable_clean, 25))

p_traj_line <- ggplot(traj_dat, aes(x = Phase_clean, y = std_diff, group = Variable_plot, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + geom_line(linewidth = 0.7, alpha = 0.7) + geom_point(size = 2) +
  facet_wrap(~ Variable_plot, scales = "free_y", ncol = 5) + scale_color_manual(values = direction_colors) + theme_pub(base_size = 8) +
  labs(title = "Trajectory of Stage D vs C differences", x = "", y = "Standardized difference")
p_traj_line

p_traj_heatmap <- ggplot(traj_dat, aes(x = Phase_clean, y = fct_reorder(Variable_plot, abs_raw_diff), fill = std_diff)) +
  geom_tile(color = "white", linewidth = 0.4) + geom_text(aes(label = p_star), size = 3) + 
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = 0) + theme_pub() + theme(legend.position = "right", axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "Pseudo-temporal laboratory trajectory", x = "", y = "", fill = "Std. diff\nD − C")
p_traj_heatmap

p_traj_dot <- ggplot(traj_dat, aes(x = std_diff, y = fct_reorder(Variable_plot, abs_raw_diff), color = Phase_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") + geom_point(size = 2.5, alpha = 0.85) + theme_pub() + theme(legend.position = "top") +
  labs(title = "Phase-specific Stage D vs C separation", x = "Standardized difference", y = "", color = "")
p_traj_dot



# ============================================================
# Supplementary Figure: Admission vs Shock
# ============================================================
clinical_data <- read_csv(clinical_file, show_col_types = FALSE)
admshock_results <- read.csv(admshock_file)
head(clinical_data)
head(admshock_results)

# Volcano stype progression plot
split_label <- function(var) {
  var <- sub("^labs_", "", var); unit_pattern <- "_(mg_dl|mmol_l|mmhg|iu_l|g_dl|ng_l|f)$"; has_unit <- grepl(unit_pattern, var)
  unit_part <- sub(paste0(".*", unit_pattern), "\\1", var); name_part <- sub(unit_pattern, "", var)
  name_clean <- clean_label(ifelse(has_unit, name_part, var)); unit_clean <- clean_label(unit_part)
  ifelse(has_unit, paste0(name_clean, "\n(", unit_clean, ")"), name_clean)
}

plot_df <- admshock_results %>% dplyr::mutate(direction = dplyr::case_when(Adjusted_P_Value < 0.05 & Median_diff > 0 ~ "Increased with shock progression", Adjusted_P_Value < 0.05 & Median_diff < 0 ~ "Decreased with shock progression", TRUE ~ "Not significant"),
                                              label = split_label(Variable), neglog10padj = -log10(Adjusted_P_Value))
progression_colors <- c("Increased with shock progression" = stage_colors[["Stage D"]], "Decreased with shock progression" = stage_colors[["Stage C"]], "Not significant" = "gray75")
p_admshock <- ggplot(plot_df, aes(x = Median_diff, y = neglog10padj, color = direction, label = label)) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) + geom_point(size = 2, alpha = 0.85) + ggrepel::geom_text_repel(size = 2.5, max.overlaps = 10, min.segment.length = 0, lineheight = 0.85) +
  scale_color_manual(values = progression_colors) + theme_pub() + theme(legend.position = "top", legend.justification = "left") +
labs(title = "Supp. 2A. Clinical trajectory: Admission to Shock Onset", x = "Change with shock progression", y = expression(-log[10]("adjusted p-value")), color = NULL)
p_admshock
ggsave(file.path(output_dir, "SuppFigure_Admission_to_Shock_progression.png"), p_admshock, width = panel_width, height = 5.2, dpi = dpi_out, bg = "white")

# Paired violin + trajectory plot, top significant variables
cap_first <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}
top_vars <- admshock_results %>% dplyr::arrange(Adjusted_P_Value) %>% dplyr::pull(Variable)
paired_long <- purrr::map_dfr(top_vars, function(v) {
  data.frame(record_id = clinical_data$record_id, Admission = clinical_data[[paste0("admission_", v)]], Shock = clinical_data[[paste0("shock_onset_", v)]], Variable = cap_first(split_label(v)))
}) %>% dplyr::filter(!is.na(Admission), !is.na(Shock)) %>% tidyr::pivot_longer(c(Admission, Shock), names_to = "Phase", values_to = "Value") %>% dplyr::mutate(Phase = factor(Phase, levels = c("Admission", "Shock")))
p_violin <- ggplot(paired_long, aes(x = Phase, y = Value)) + geom_violin(aes(fill = Phase), alpha = 0.5, color = NA, trim = FALSE) + geom_line(aes(group = record_id), color = "grey50", alpha = 0.35, linewidth = 0.3) +
  geom_point(aes(color = Phase), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = stage_colors %>% setNames(c("Admission", "Shock"))) + scale_color_manual(values = stage_colors %>% setNames(c("Admission", "Shock"))) + facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  labs(title = "Supp. 1A. Paired clinical trajectories: Admission to Shock Onset", subtitle = "Lines connect the same patient", x = NULL, y = NULL) + theme_pub_small() + theme(legend.position = "none", strip.text = element_text(size = 6))
p_violin
ggsave(file.path(output_dir, "SuppFigure_Admission_Shock_violins.png"), p_violin, width = panel_width, height = 8, dpi = dpi_out, bg = "white")

# Radar chart: standardized change by organ-system domain
domain_map <- c(
  labs_creatinine_mg_dl = "Renal", labs_blood_urea_nitrogen_bun_mg_dl = "Renal", labs_mdrd_egfr = "Renal",
  labs_bilirubin_total_mg_dl = "Hepatic", labs_bilirubin_direct_mg_dl = "Hepatic", labs_bilirubin_indirect_mg_dl = "Hepatic",
  labs_alanine_aminotransferase_iu_l = "Hepatic", labs_alkaline_phosphatase_iu_l = "Hepatic",
  labs_albumin_level_g_dl = "Hepatic", labs_protein_total_g_dl = "Hepatic",
  labs_hemoglobin_g_dl = "Hematologic", labs_platelet = "Hematologic",
  labs_white_blood_cell = "Hematologic", labs_neutrophil_absolute_count = "Hematologic",
  labs_po2_arterial_mmhg = "Respiratory", labs_pco2_arterial_mmhg = "Respiratory", labs_ph_arterial = "Respiratory",
  labs_fio2 = "Respiratory", labs_po2_fio2_ratio = "Respiratory", labs_measured_o2_saturation_arterial = "Respiratory",
  labs_measured_o2_saturation_venous = "Respiratory", labs_a_a_o2_ratio = "Respiratory",
  labs_aa_o2_difference_mmhg = "Respiratory", labs_spo2 = "Respiratory",
  labs_bicarbonate_arterial_mmol_l = "Respiratory", labs_base_excess_arterial_mmol_l = "Respiratory",
  labs_systolic_blood_pressure_mmhg = "Hemodynamic", labs_diastolic_blood_pressure_mmhg = "Hemodynamic",
  labs_ekg_rate = "Hemodynamic", labs_respirations = "Hemodynamic", labs_temperature_f = "Hemodynamic",
  labs_glucose_level_mg_dl = "Metabolic", labs_sodium_level_mmol_l = "Metabolic", labs_potassium_level_mmol_l = "Metabolic",
  labs_chloride_mmol_l = "Metabolic", labs_calcium_level_total_mg_dl = "Metabolic", labs_phosphorus_level_mg_dl = "Metabolic",
  labs_anion_gap_mmol_l = "Metabolic", labs_bicarbonate_level_mmol_l = "Metabolic",
  labs_lactate_mmol_l = "Perfusion/Stress", labs_nt_pro_bnp_ng_l = "Perfusion/Stress", labs_c_reactive_protein_mg_l = "Perfusion/Stress"
)

radar_df <- purrr::map_dfr(names(domain_map), function(v) {
  adm <- clinical_data[[paste0("admission_", v)]]; shk <- clinical_data[[paste0("shock_onset_", v)]]
  complete <- !is.na(adm) & !is.na(shk)
  if (sum(complete) < 5 || sd(adm[complete], na.rm = TRUE) == 0) return(NULL)
  z <- (median(shk[complete]) - median(adm[complete])) / sd(adm[complete], na.rm = TRUE)
  data.frame(Variable = v, Domain = domain_map[[v]], abs_z = abs(z))
})
domain_summary <- radar_df %>% dplyr::group_by(Domain) %>% dplyr::summarise(mean_abs_z = mean(abs_z), n_vars = dplyr::n(), .groups = "drop") %>% dplyr::arrange(Domain)
n_domains <- nrow(domain_summary)
domain_summary <- domain_summary %>% dplyr::arrange(Domain)
domain_summary$angle <- seq(pi/2, pi/2 - 2*pi, length.out = n_domains + 1)[1:n_domains]

polygon_pts <- domain_summary %>% dplyr::mutate(x = mean_abs_z * cos(angle), y = mean_abs_z * sin(angle)) %>% dplyr::bind_rows(dplyr::slice(., 1))
r_max <- max(domain_summary$mean_abs_z) * 1.4
grid_rings <- pretty(c(0, r_max), n = 4); grid_rings <- grid_rings[grid_rings > 0]
circle_pts <- purrr::map_dfr(grid_rings, function(r) {th <- seq(0, 2*pi, length.out = 200); data.frame(x = r * cos(th), y = r * sin(th), r = r)})
spoke_pts <- domain_summary %>% dplyr::mutate(x = r_max * cos(angle), y = r_max * sin(angle), x0 = 0, y0 = 0)
label_pts <- domain_summary %>% dplyr::mutate(x = (r_max * 1.12) * cos(angle), y = (r_max * 1.12) * sin(angle), label = paste0(Domain, "\n(", n_vars, " vars)"))
ring_labels <- data.frame(x = 0.05, y = grid_rings, label = grid_rings)
  
p_radar <- ggplot() + geom_path(data = circle_pts, aes(x = x, y = y, group = r), color = "grey85", linewidth = 0.3) + geom_segment(data = spoke_pts, aes(x = x0, y = y0, xend = x, yend = y), color = "grey85", linewidth = 0.3) + 
  geom_polygon(data = polygon_pts, aes(x = x, y = y), fill = stage_colors[["Stage D"]], alpha = 0.25, color = stage_colors[["Stage D"]], linewidth = 0.9) +
  geom_point(data = dplyr::filter(polygon_pts, !duplicated(Domain)), aes(x = x, y = y), color = stage_colors[["Stage D"]], size = 2.2) +
  geom_text(data = label_pts, aes(x = x, y = y, label = label), size = 2.4, fontface = "bold", lineheight = 0.8) + geom_text(data = ring_labels, aes(x = x, y = y, label = label), size = 2, color = "grey50", hjust = 0) +
  coord_equal(clip = "off") + theme_pub() + 
  labs(title = "Supp. 2B. Systemic derangement, Admission to Shock Onset", subtitle = "Mean |standardized change| per organ-system domain (z relative to admission SD)") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), panel.background = element_blank(),
                      plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(15, 45, 15, 45))
p_radar
ggsave(file.path(output_dir, "SuppFigure_Admission_Shock_radar.png"), p_radar, width = panel_width, height = 6, dpi = dpi_out, bg = "white")




admshock_res <- readRDS("cs_project/clinical_data/clinical_results/CS_Admission_vs_Shock_byStage.rds")
admshock_by_stage <- admshock_res$admshock_by_stage
clinical_data <- admshock_res$clinical_data

stage_col <- "consent_scai_stage_at_time_of_informed_consent"
split_label <- function(var) {
  var <- sub("^labs_", "", var); unit_pattern <- "_(mg_dl|mmol_l|mmhg|iu_l|g_dl|ng_l|f)$"; has_unit <- grepl(unit_pattern, var)
  unit_part <- sub(paste0(".*", unit_pattern), "\\1", var); name_part <- sub(unit_pattern, "", var)
  name_clean <- clean_label(ifelse(has_unit, name_part, var)); unit_clean <- clean_label(unit_part)
  ifelse(has_unit, paste0(name_clean, "\n(", unit_clean, ")"), name_clean)
}
cap_first <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x}

# Volcano plot, faceted by Stage
set.seed(42)
min_n_volcano <- 5
plot_df <- admshock_by_stage %>% dplyr::filter(N_paired >= min_n_volcano) %>%
  dplyr::mutate(direction = dplyr::case_when(Adjusted_P_Value < 0.05 & Median_diff > 0 ~ "Increased with shock progression", Adjusted_P_Value < 0.05 & Median_diff < 0 ~ "Decreased with shock progression", TRUE ~ "Not significant"),
                label = split_label(Variable_), neglog10padj = -log10(Adjusted_P_Value), y_nudge = neglog10padj + runif(dplyr::n(), -0.03, 0.03))
progression_colors <- c("Increased with shock progression" = stage_colors[["Stage D"]], "Decreased with shock progression" = stage_colors[["Stage C"]], "Not significant" = "gray75")
p_admshock <- ggplot(plot_df, aes(x = Median_diff, y = neglog10padj, color = direction, label = label)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) + geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) + geom_point(size = 2.5, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(y = y_nudge), size = 3, max.overlaps = Inf, min.segment.length = 0.1, box.padding = 0.6, point.padding = 0.4, force = 6, force_pull = 0.2, max.time = 3, max.iter = 30000, segment.size = 0.3, segment.color = "grey50", segment.curvature = -0.1, segment.ncp = 3, seed = 42, lineheight = 1) +
  scale_color_manual(values = progression_colors) + scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) + facet_wrap(~ Stage, ncol = 2, scales = "free_x") +
  theme_pub() + theme(legend.position = "top", legend.justification = "left") +
  labs(title = "Supp. 2A. Clinical trajectory: Admission to Shock Onset, by SCAI Stage", x = "Change with shock progression", y = expression(-log[10]("adjusted p-value")), color = NULL)
p_admshock
ggsave(file.path(output_dir, "SuppFigure_Admission_to_Shock_progression_byStage.png"), p_admshock, width = panel_width * 2, height = 5.2, dpi = dpi_out, bg = "white")

# Paired violin + trajectory plot, top significant variables, split by Stage (top significant vars, per stage's own ranking)
build_paired_long <- function(stage_label, top_vars) {
  ids <- clinical_data %>% dplyr::filter(.data[[stage_col]] == stage_label) %>% dplyr::pull(record_id)
  purrr::map_dfr(top_vars, function(v) {
    data.frame(record_id = ids, Admission = clinical_data[[paste0("admission_", v)]][clinical_data$record_id %in% ids], Shock = clinical_data[[paste0("shock_onset_", v)]][clinical_data$record_id %in% ids], Variable = cap_first(split_label(v)))}) %>%
    dplyr::filter(!is.na(Admission), !is.na(Shock)) %>% tidyr::pivot_longer(c(Admission, Shock), names_to = "Phase", values_to = "Value") %>% dplyr::mutate(Phase = factor(Phase, levels = c("Admission", "Shock")), Stage = stage_label)
}
top_vars_C <- admshock_res$admshock_C %>% dplyr::arrange(Adjusted_P_Value) %>% dplyr::pull(Variable_)
top_vars_D <- admshock_res$admshock_D %>% dplyr::arrange(Adjusted_P_Value) %>% dplyr::pull(Variable_)
top_vars_both <- intersect(top_vars_C, top_vars_D)
paired_long_all <- dplyr::bind_rows(build_paired_long("C", top_vars_both), build_paired_long("D", top_vars_both)) %>% dplyr::mutate(Stage = dplyr::recode(Stage, "C" = "Stage C", "D" = "Stage D"), record_id_stage = paste0(record_id, "_", Stage))
variables <- unique(paired_long_all$Variable)
pval_lookup <- admshock_by_stage %>% dplyr::filter(Variable_ %in% top_vars_both) %>% 
  dplyr::mutate(Variable = cap_first(split_label(Variable_)), Stage = dplyr::recode(Stage, "C" = "Stage C", "D" = "Stage D"), label = dplyr::if_else(is.na(Adjusted_P_Value) | is.na(P_Value), "ns", sig_star(Adjusted_P_Value))) %>%
  dplyr::mutate(label = dplyr::if_else(label == "", "ns", label)) %>% dplyr::select(Variable, Stage, label)
all_combos <- tidyr::expand_grid(Variable = variables, Stage = c("Stage C", "Stage D"))
pval_lookup <- all_combos %>% dplyr::left_join(pval_lookup, by = c("Variable", "Stage")) %>% dplyr::mutate(label = dplyr::if_else(is.na(label), "ns", label))

make_var_panel <- function(var_name) {
  df <- paired_long_all %>% dplyr::filter(Variable == var_name)
  pvals <- pval_lookup %>% dplyr::filter(Variable == var_name)
  ggplot(df, aes(x = Phase, y = Value)) + geom_violin(aes(fill = Phase), alpha = 0.5, color = NA, trim = FALSE) + geom_line(aes(group = record_id), color = "grey50", alpha = 0.35, linewidth = 0.3) +
    geom_point(aes(color = Phase), size = 1, alpha = 0.7) + geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label, fontface = ifelse(label == "ns", "plain", "bold")), inherit.aes = FALSE, vjust = 1.3, size = 2) +
    scale_fill_manual(values = stage_colors %>% setNames(c("Admission", "Shock"))) + scale_color_manual(values = stage_colors %>% setNames(c("Admission", "Shock"))) +
    facet_wrap(~ Stage, scales = "free_y") + labs(title = var_name, x = NULL, y = NULL) + theme_pub_small() +
    theme(legend.position = "none", strip.text = element_text(size = 6), plot.title = element_text(size = 7, hjust = 0.5), plot.background = element_rect(color = "grey40", fill = NA, linewidth = 0.5), plot.margin = margin(4, 4, 4, 4))
}
var_panels <- purrr::map(variables, make_var_panel)
p_violin_combined <- patchwork::wrap_plots(var_panels, ncol = 3) + patchwork::plot_annotation(title = "Supp. 1A. Paired clinical trajectories: Admission to Shock Onset", subtitle = "Lines connect the same patient", theme = theme_pub_small())
p_violin_combined
ggsave(file.path(output_dir, "SuppFigure_Admission_Shock_violins_byStage.png"), p_violin_combined, width = panel_width, height = 10, dpi = dpi_out, bg = "white")


# Radar chart: standardized change by organ-system domain, Stages overlaid on one plot
domain_map <- c(
  labs_creatinine_mg_dl = "Renal", labs_blood_urea_nitrogen_bun_mg_dl = "Renal", labs_mdrd_egfr = "Renal",
  labs_bilirubin_total_mg_dl = "Hepatic", labs_bilirubin_direct_mg_dl = "Hepatic", labs_bilirubin_indirect_mg_dl = "Hepatic",
  labs_alanine_aminotransferase_iu_l = "Hepatic", labs_alkaline_phosphatase_iu_l = "Hepatic",
  labs_albumin_level_g_dl = "Hepatic", labs_protein_total_g_dl = "Hepatic",
  labs_hemoglobin_g_dl = "Hematologic", labs_platelet = "Hematologic",
  labs_white_blood_cell = "Hematologic", labs_neutrophil_absolute_count = "Hematologic",
  labs_po2_arterial_mmhg = "Respiratory", labs_pco2_arterial_mmhg = "Respiratory", labs_ph_arterial = "Respiratory",
  labs_fio2 = "Respiratory", labs_po2_fio2_ratio = "Respiratory", labs_measured_o2_saturation_arterial = "Respiratory",
  labs_measured_o2_saturation_venous = "Respiratory", labs_a_a_o2_ratio = "Respiratory",
  labs_aa_o2_difference_mmhg = "Respiratory", labs_spo2 = "Respiratory",
  labs_bicarbonate_arterial_mmol_l = "Respiratory", labs_base_excess_arterial_mmol_l = "Respiratory",
  labs_systolic_blood_pressure_mmhg = "Hemodynamic", labs_diastolic_blood_pressure_mmhg = "Hemodynamic",
  labs_ekg_rate = "Hemodynamic", labs_respirations = "Hemodynamic", labs_temperature_f = "Hemodynamic",
  labs_glucose_level_mg_dl = "Metabolic", labs_sodium_level_mmol_l = "Metabolic", labs_potassium_level_mmol_l = "Metabolic",
  labs_chloride_mmol_l = "Metabolic", labs_calcium_level_total_mg_dl = "Metabolic", labs_phosphorus_level_mg_dl = "Metabolic",
  labs_anion_gap_mmol_l = "Metabolic", labs_bicarbonate_level_mmol_l = "Metabolic",
  labs_lactate_mmol_l = "Perfusion/Stress", labs_nt_pro_bnp_ng_l = "Perfusion/Stress", labs_c_reactive_protein_mg_l = "Perfusion/Stress"
)

compute_domain_summary <- function(stage_label) {
  ids <- clinical_data %>% dplyr::filter(.data[[stage_col]] == stage_label) %>% dplyr::pull(record_id)
  data_sub <- clinical_data %>% dplyr::filter(record_id %in% ids)
  radar_df <- purrr::map_dfr(names(domain_map), function(v) {
    adm <- data_sub[[paste0("admission_", v)]]; shk <- data_sub[[paste0("shock_onset_", v)]]
    complete <- !is.na(adm) & !is.na(shk)
    if (sum(complete) < 5 || sd(adm[complete], na.rm = TRUE) == 0) return(NULL)
    z <- (median(shk[complete]) - median(adm[complete])) / sd(adm[complete], na.rm = TRUE)
    data.frame(Variable = v, Domain = domain_map[[v]], abs_z = abs(z))
  })
  radar_df %>% dplyr::group_by(Domain) %>% dplyr::summarise(mean_abs_z = mean(abs_z), n_vars = dplyr::n(), .groups = "drop") %>% dplyr::arrange(Domain) %>% dplyr::mutate(Stage = stage_label)
}

domain_summary_C <- compute_domain_summary("C"); domain_summary_D <- compute_domain_summary("D")
all_domains <- sort(union(domain_summary_C$Domain, domain_summary_D$Domain))
n_domains <- length(all_domains)
angle_lookup <- setNames(seq(pi/2, pi/2 - 2*pi, length.out = n_domains + 1)[1:n_domains], all_domains)
prep_domain_summary <- function(df) {
  df %>% dplyr::right_join(data.frame(Domain = all_domains), by = "Domain") %>% dplyr::mutate(mean_abs_z = tidyr::replace_na(mean_abs_z, 0), n_vars = tidyr::replace_na(n_vars, 0), angle = angle_lookup[Domain]) %>% dplyr::arrange(angle)
}
domain_summary_C <- prep_domain_summary(domain_summary_C) %>% dplyr::mutate(Stage = "Stage C"); domain_summary_D <- prep_domain_summary(domain_summary_D) %>% dplyr::mutate(Stage = "Stage D")
build_polygon <- function(df) df %>% dplyr::mutate(x = mean_abs_z * cos(angle), y = mean_abs_z * sin(angle)) %>% dplyr::bind_rows(dplyr::slice(., 1))
polygon_C <- build_polygon(domain_summary_C); polygon_D <- build_polygon(domain_summary_D); polygon_both <- dplyr::bind_rows(polygon_C, polygon_D)
r_max <- max(c(domain_summary_C$mean_abs_z, domain_summary_D$mean_abs_z), na.rm = TRUE) * 1.4
grid_rings <- pretty(c(0, r_max), n = 4); grid_rings <- grid_rings[grid_rings > 0]
circle_pts <- purrr::map_dfr(grid_rings, function(r) {th <- seq(0, 2*pi, length.out = 200); data.frame(x = r * cos(th), y = r * sin(th), r = r)})
spoke_pts <- domain_summary_C %>% dplyr::mutate(x = r_max * cos(angle), y = r_max * sin(angle), x0 = 0, y0 = 0)
label_pts <- domain_summary_C %>% dplyr::mutate(x = (r_max * 1.12) * cos(angle), y = (r_max * 1.12) * sin(angle), label = Domain)
ring_labels <- data.frame(x = 0.05, y = grid_rings, label = grid_rings)

p_radar <- ggplot() + geom_path(data = circle_pts, aes(x = x, y = y, group = r), color = "grey85", linewidth = 0.3) + geom_segment(data = spoke_pts, aes(x = x0, y = y0, xend = x, yend = y), color = "grey85", linewidth = 0.3) +
  geom_polygon(data = polygon_both, aes(x = x, y = y, fill = Stage, color = Stage), alpha = 0.2, linewidth = 0.9) + geom_point(data = dplyr::filter(polygon_both, !duplicated(paste(Stage, Domain))), aes(x = x, y = y, color = Stage), size = 2.2) +
  scale_fill_manual(values = stage_colors) + scale_color_manual(values = stage_colors) + 
  geom_text(data = label_pts, aes(x = x, y = y, label = label), size = 2.4, fontface = "bold", lineheight = 0.8) + geom_text(data = ring_labels, aes(x = x, y = y, label = label), size = 2, color = "grey50", hjust = 0) +
  coord_equal(clip = "off") + theme_pub() + labs(title = "Supp. 2B. Systemic derangement, Admission to Shock Onset", subtitle = "Mean |standardized change| per organ-system domain (z relative to admission SD)") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), 
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top", legend.justification = "center", plot.margin = margin(15, 45, 15, 45))
p_radar
ggsave(file.path(output_dir, "SuppFigure_Admission_Shock_radar_byStage.png"), p_radar, width = panel_width, height = 6, dpi = dpi_out, bg = "white")

