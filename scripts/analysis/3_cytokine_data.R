setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(tidyr)
library(broom)
library(forcats)

# 1. Load and clean data
data <- read_excel("cs_project/cytokine_data/CScytokineData.xlsx")
data <- data %>% mutate(`SCAI Stage` = toupper(`SCAI Stage`), Sex = tolower(Sex), Race = tolower(Race)) %>%
  rename_with(~ str_replace_all(., "-", "")) %>% rename_with(~ str_remove_all(., " \\(pg/ml\\)")) %>% rename_with(~ str_replace_all(., " ", ""))
colnames(data)
str(data)
summary(data)

data_cd <- data %>% filter(SCAIStage %in% c("C", "D")) %>% mutate(SCAIStage = factor(SCAIStage, levels = c("C", "D")))
# table(data_cd$SCAIStage)
# table(data_cd$Sex)
# table(data_cd$Race)
# range(data_cd$Age, na.rm = TRUE)

meta_cols <- c("StudyID", "Age", "Sex", "Race", "SCAIStage", "SRspecimenID")
cytokines <- setdiff(colnames(data_cd), meta_cols)
cytokines

qc_table <- data.frame(Cytokine = cytokines, MissingN = sapply(data_cd[cytokines], function(x) sum(is.na(x))), UniqueN  = sapply(data_cd[cytokines], function(x) length(unique(na.omit(x)))), 
                       NonZeroN = sapply(data_cd[cytokines], function(x) sum(x > 0, na.rm = TRUE))) %>% arrange(UniqueN, MissingN, NonZeroN)

cytokines_use <- qc_table %>% filter(UniqueN > 3) %>% pull(Cytokine) # Refined filter: keep cytokines with enough variability
cytokines_use
data_cd_log <- data_cd %>% mutate(across(all_of(cytokines_use), ~ log1p(.))) # Log-transform cytokines

# 2. D vs C: Univariate comparison (Wilcoxon test)
results <- map_dfr(cytokines_use, function(cyt) {
  test <- wilcox.test(data_cd_log[[cyt]] ~ data_cd_log$SCAIStage, exact = FALSE)
  
  median_C_log <- median(data_cd_log[[cyt]][data_cd_log$SCAIStage == "C"], na.rm = TRUE); median_D_log <- median(data_cd_log[[cyt]][data_cd_log$SCAIStage == "D"], na.rm = TRUE)
  median_C_raw <- median(data_cd[[cyt]][data_cd$SCAIStage == "C"], na.rm = TRUE); median_D_raw <- median(data_cd[[cyt]][data_cd$SCAIStage == "D"], na.rm = TRUE)
  tibble(Cytokine = cyt, p_value = test$p.value, median_C_log = median_C_log, median_D_log = median_D_log, median_C_raw = median_C_raw, median_D_raw = median_D_raw)
}) %>% mutate(p_adj = p.adjust(p_value, method = "BH"), direction = case_when(median_D_log > median_C_log ~ "Higher in D", median_D_log < median_C_log ~ "Higher in C", TRUE ~ "No difference"),
              fold_change_D_vs_C = (median_D_raw + 1) / (median_C_raw + 1)) %>% arrange(p_adj, p_value)
results

report_table <- results %>% mutate(median_C_raw = round(median_C_raw, 2), median_D_raw = round(median_D_raw, 2), median_C_log = round(median_C_log, 2), median_D_log = round(median_D_log, 2),
                                   fold_change_D_vs_C = round(fold_change_D_vs_C, 2), p_value = signif(p_value, 3), p_adj = signif(p_adj, 3)) %>% arrange(p_adj, p_value)
report_table
write.csv(report_table, "cs_project/cytokine_data/results/cytokine_univariate_results_C_vs_D.csv", row.names = FALSE)

# 3. Top cytokines for plotting
top_hits <- results %>% arrange(p_adj, p_value) %>% slice_head(n = min(6, nrow(results))) %>% pull(Cytokine)
top_hits

plot_data <- data_cd_log %>% dplyr::select(SCAIStage, all_of(top_hits)) %>% pivot_longer(cols = all_of(top_hits), names_to = "Cytokine", values_to = "LogValue")
plot_data$Cytokine <- factor(plot_data$Cytokine, levels = top_hits)

plot_annot <- results %>% filter(Cytokine %in% top_hits) %>% mutate(label = paste0("adj. p = ", signif(p_adj, 3), "\nraw p = ", signif(p_value, 3)))
y_positions <- plot_data %>% group_by(Cytokine) %>% summarise(y = max(LogValue, na.rm = TRUE) * 1.08)
plot_annot <- plot_annot %>% left_join(y_positions, by = "Cytokine")

p_box <- ggplot(plot_data, aes(x = SCAIStage, y = LogValue, fill = SCAIStage)) + geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(color = SCAIStage), width = 0.12, alpha = 0.8, size = 2, show.legend = FALSE) + facet_wrap(~ Cytokine, scales = "free_y", ncol = 3) +
  geom_text(data = plot_annot, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
  scale_fill_manual(values = c("C" = "steelblue", "D" = "darkred")) + scale_color_manual(values = c("C" = "steelblue", "D" = "darkred")) +
  labs(title = "Top cytokine differences between SCAI Stage C and D", x = "SCAI Stage", y = "log1p(concentration)") + theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "white", colour = "black"), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"), axis.text = element_text(colour = "black"), legend.position = "none")
p_box
ggsave(filename = "cs_project/cytokine_data/results/top_cytokines_boxplots.png", plot = p_box, width = 11, height = 7.5, dpi = 600)
ggsave(filename = "cs_project/cytokine_data/results/top_cytokines_boxplots.pdf", plot = p_box, width = 11, height = 7.5)

# 4. volcano-style plot
volcano_data <- results %>% mutate(log2_fold_change = log2(fold_change_D_vs_C), neg_log10_p = -log10(p_value), 
                                   direction_color = case_when(log2_fold_change > 0 ~ "Up in D", log2_fold_change < 0 ~ "Up in C", TRUE ~ "No change"), 
                                   label_flag = rank(p_adj, ties.method = "first") <= 10)

p_volcano <- ggplot(volcano_data, aes(x = log2_fold_change, y = neg_log10_p)) + geom_point(aes(color = direction_color), size = 3, alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  geom_text(data = subset(volcano_data, label_flag), aes(label = Cytokine), vjust = -0.7, size = 3.5) +
  scale_color_manual(values = c("Up in D" = "darkred", "Up in C" = "steelblue", "No change" = "grey50")) +
  labs(title = "Univariate cytokine differences: Stage D vs Stage C", x = "log2 fold change (D vs C)", y = "-log10(p-value)", color = "Direction") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
                                   axis.title = element_text(face = "bold"), axis.text = element_text(colour = "black"), legend.position = "right")
p_volcano
ggsave(filename = "cs_project/cytokine_data/results/volcano_cytokines.png", plot = p_volcano, width = 8, height = 6, dpi = 600)
ggsave(filename = "cs_project/cytokine_data/results/volcano_cytokines.pdf", plot = p_volcano, width = 8, height = 6)


# # 5. PCA
# library(FactoMineR); library(factoextra)
# 
# pca_matrix <- data_cd_log %>% dplyr::select(all_of(cytokines_use)) %>% scale()
# pca_res <- PCA(pca_matrix, graph = FALSE)
# p_pca <- fviz_pca_ind(pca_res, habillage = data_cd_log$SCAIStage, addEllipses = TRUE, ellipse.level = 0.95, 
#                       palette = c("steelblue", "darkred"), title = "PCA of cytokine profiles (Stage C vs D)")
# p_pca
# ggsave("cs_project/cytokine_data/results/PCA.png", p_pca, width = 7, height = 6, dpi = 600)

# p_loadings <- fviz_pca_var(pca_res, col.var = "contrib", gradient.cols = c("grey70", "orange", "red"), repel = TRUE) + ggtitle("Cytokine contributions to PCA")
# p_loadings
# ggsave("cs_project/cytokine_data/results/PCA_loadings.png", p_loadings, width = 7, height = 6, dpi = 600)

#========================
# 6. Predictive models
#========================
library(pROC)

y <- ifelse(data_cd_log$SCAIStage == "D", 1, 0)

model_IL6 <- glm(SCAIStage ~ ST2 + IL6, data = data_cd_log, family = binomial)
model_IL8 <- glm(SCAIStage ~ ST2 + IL8, data = data_cd_log, family = binomial)
model_S100A8 <- glm(SCAIStage ~ ST2 + S100A8, data = data_cd_log, family = binomial)
model_Trop <- glm(SCAIStage ~ ST2 + TroponinI, data = data_cd_log, family = binomial)
model_VEGFA <- glm(SCAIStage ~ ST2 + VEGFA, data = data_cd_log, family = binomial)

roc_st2 <- roc(y, data_cd_log$ST2) # Baseline
roc_IL6 <- roc(y, predict(model_IL6, type = "response"))
roc_IL8 <- roc(y, predict(model_IL8, type = "response"))
roc_S100A8 <- roc(y, predict(model_S100A8, type = "response"))
roc_Trop <- roc(y, predict(model_Trop, type = "response"))
roc_VEGFA <- roc(y, predict(model_VEGFA, type = "response"))

roc.test(roc_st2, roc_VEGFA)
auc_table <- data.frame(Model = c("ST2", "ST2 + IL6", "ST2 + IL8", "ST2 + S100A8", "ST2 + TroponinI", "ST2 + VEGFA"),
                        AUC = c(auc(roc_st2), auc(roc_IL6), auc(roc_IL8), auc(roc_S100A8), auc(roc_Trop), auc(roc_VEGFA)))
auc_table
write.csv(auc_table, "cs_project/cytokine_data/results/auc_comparison_models.csv", row.names = FALSE)
  
roc_list <- list(roc_st2, roc_IL6, roc_IL8, roc_S100A8, roc_Trop, roc_VEGFA)
names(roc_list) <- c(paste0("ST2 (AUC=", round(auc(roc_st2), 3), ")"), paste0("ST2 + IL6 (AUC=", round(auc(roc_IL6), 3), ")"),
                     paste0("ST2 + IL8 (AUC=", round(auc(roc_IL8), 3), ")"), paste0("ST2 + S100A8 (AUC=", round(auc(roc_S100A8), 3), ")"),
                     paste0("ST2 + TroponinI (AUC=", round(auc(roc_Trop), 3), ")"), paste0("ST2 + VEGFA (AUC=", round(auc(roc_VEGFA), 3), ")"))
p_roc_all <- ggroc(roc_list, size = 1.2) + scale_color_manual(values = c("darkred", "black", "steelblue", "purple", "orange", "darkgreen")) + theme_bw() +
  labs(title = "ROC comparison of ST2-based models", x = "False Positive Rate", y = "True Positive Rate", color = "Model")
p_roc_all
ggsave("cs_project/cytokine_data/results/ROC_all_models.png", p_roc_all, width = 7, height = 6, dpi = 600)
ggsave("cs_project/cytokine_data/results/ROC_all_models.pdf", p_roc_all, width = 7, height = 6)


# 7. Confounding interactive analysis
table(data_cd_log$Sex)
table(data_cd_log$SCAIStage, data_cd_log$Sex)
table(data_cd_log$Race)

cytokines <- cytokines_use
results_list <- list()
for (cyt in cytokines) {
  model_sex <- glm(as.formula(paste(cyt, "~ SCAIStage * Sex")), data = data_cd_log)
  sex_p <- summary(model_sex)$coefficients["SCAIStageD:Sexm", "Pr(>|t|)"]
  
  model_age <- glm(as.formula(paste(cyt, "~ SCAIStage * Age")), data = data_cd_log)
  age_p <- summary(model_age)$coefficients["SCAIStageD:Age", "Pr(>|t|)"]
  
  results_list[[cyt]] <- data.frame(Cytokine = cyt, Sex_interaction_p = sex_p, Age_interaction_p = age_p)
}

results_df <- bind_rows(results_list)
results_df <- results_df %>% 
  mutate(Sex_signal = case_when(Sex_interaction_p < 0.05 ~ "Significant", Sex_interaction_p < 0.15 ~ "Trend", TRUE ~ "None"),
         Age_signal = case_when(Age_interaction_p < 0.05 ~ "Significant", Age_interaction_p < 0.15 ~ "Trend", TRUE ~ "None")) %>% arrange(Sex_interaction_p)
results_df
write.csv(results_df, "cs_project/cytokine_data/results/cytokine_confounding interactive_results_C_vs_D.csv", row.names = FALSE)

# “ST2 levels increased from Stage C to Stage D in both sexes, with a trend toward a greater increase in females compared with males.”
p_st2_sex <- ggplot(data_cd_log, aes(x = SCAIStage, y = ST2, color = Sex)) + 
  geom_boxplot(aes(group = interaction(SCAIStage, Sex)), position = position_dodge(width = 0.75), width = 0.45, fill = "lightgrey", linewidth = 1) +
  stat_summary(aes(group = interaction(SCAIStage, Sex)), fun = mean, geom = "point", position = position_dodge(width = 0.75), size = 3, shape = 18, color = "black") +
  scale_color_manual(values = c("f" = "brown", "m" = "darkgreen")) + theme_bw() +
  labs(title = "ST2", x = "SCAI Stage", y = "log1p(ST2)", color = "Sex")
p_st2_sex
ggsave("cs_project/cytokine_data/results/ST2_stage_sex.png", p_st2_sex, width = 6, height = 5, dpi = 600)
ggsave("cs_project/cytokine_data/results/ST2_stage_sex.pdf", p_st2_sex, width = 6, height = 5)

# “There was a trend toward interaction between age and SCAI stage for S100A8, suggesting that the relationship between age and S100A8 may differ by disease severity”
p_s100a8_age <- ggplot(data_cd_log, aes(x = Age, y = S100A8, color = SCAIStage)) + geom_point(size = 2) + geom_smooth(method = "lm", se = TRUE) + theme_bw() +
  scale_color_manual(values = c("C" = "steelblue", "D" = "darkred")) + labs(title = "S100A8", x = "Age", y = "S100A8")
p_s100a8_age
ggsave("cs_project/cytokine_data/results/S100A8_age_stage.png", p_s100a8_age, width = 6, height = 5, dpi = 600)
ggsave("cs_project/cytokine_data/results/S100A8_age_stage.pdf", p_s100a8_age, width = 6, height = 5)
