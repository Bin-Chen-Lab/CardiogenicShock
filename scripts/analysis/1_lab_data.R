setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))

library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(dplyr)
library(openxlsx)

# ============================================================
# C vs D analysis
# ============================================================
input_file <- "cs_project/clinical_data/CS_full_clinical_data_set.csv"
output_dir <- "cs_project/clinical_data/clinical_results"

if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
fill_right <- function(x) {
  last <- NA_character_
  out <- x
  for (i in seq_along(x)) {
    if (!is.na(x[i]) && trimws(x[i]) != "") {last <- trimws(x[i])}
    out[i] <- last
  }
  out
}
clean_text <- function(x) {x <- trimws(as.character(x)); x[x == ""] <- NA; x}
clean_name <- function(x) {x <- tolower(x); x <- gsub("\\(choice\\s*=\\s*", "__", x); x <- gsub("[()]", "", x); x <- gsub("[/&,-]", " ", x); x <- gsub("[^a-z0-9]+", "_", x); x <- gsub("_+", "_", x); x <- gsub("^_|_$", "", x); x}
format_p <- function(x) {ifelse(is.na(x), NA, ifelse(x < 0.0001, "<0.0001", formatC(x, format = "f", digits = 4)))}
format_ci <- function(low, high) {ifelse(is.na(low) | is.na(high), NA, sprintf("[%.2f, %.2f]", low, high))}

raw <- read.csv(input_file, header = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "N/A"))

phase_row <- fill_right(as.character(raw[1, ]))
field_row <- as.character(raw[2, ])
phase_row <- clean_text(phase_row)
field_row <- clean_text(field_row)
phase_row[is.na(phase_row)] <- "unassigned"
missing_fields <- which(is.na(field_row))
field_row[missing_fields] <- paste0("unnamed_", missing_fields)
combined_names <- paste(phase_row, field_row, sep = "__")
combined_names <- clean_name(combined_names)
combined_names <- make.unique(combined_names, sep = "__dup_")

clinical_data <- raw[-c(1, 2), ]
colnames(clinical_data) <- combined_names
rownames(clinical_data) <- NULL
clinical_data[] <- lapply(clinical_data, function(x) {x <- trimws(as.character(x)); x[x == ""] <- NA; x})

if ("unassigned_record_id" %in% names(clinical_data)) {names(clinical_data)[names(clinical_data) == "unassigned_record_id"] <- "record_id"}
if ("unassigned_share_specimen_id" %in% names(clinical_data)) {names(clinical_data)[names(clinical_data) == "unassigned_share_specimen_id"] <- "specimen_id"}
names(clinical_data) <- gsub("^shock_onset_data_", "shock_onset_", names(clinical_data))
names(clinical_data) <- gsub("^vitals_and_labs_at_admission_", "admission_labs_", names(clinical_data))
names(clinical_data) <- gsub("^vitals_and_labs_at_shock_onset_", "shock_onset_labs_", names(clinical_data))
names(clinical_data) <- gsub("^echo_at_shock_onset_", "shock_onset_echo_", names(clinical_data))
names(clinical_data) <- gsub("^clinical_outcomes_and_complications_inpatient_", "inpatient_", names(clinical_data))
names(clinical_data) <- gsub("^inhospital_cardiovascular_procedures_", "procedures_", names(clinical_data))
names(clinical_data) <- gsub("^rhc_at_shock_onset_", "shock_onset_rhc_", names(clinical_data))
names(clinical_data) <- gsub("^pac_hemodynamics_at_shock_onset_", "shock_onset_pac_", names(clinical_data))
names(clinical_data) <- gsub("^pac_hemodynamics_at_time_of_mcs_", "mcs_pac_", names(clinical_data))
names(clinical_data) <- gsub("^pac_hemodynamics_at_shock_recovery_", "recovery_pac_", names(clinical_data))
names(clinical_data) <- gsub("^post_discharge_outcomes_", "post_discharge_", names(clinical_data))

clinical_data[] <- lapply(clinical_data, function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "Unknown", "unknown", "Not available")] <- NA
  x[x %in% c("Checked", "checked", "Yes", "yes", "YES", "Y", "True", "TRUE")] <- "1"
  x[x %in% c("Unchecked", "unchecked", "No", "no", "NO", "N", "False", "FALSE")] <- "0"
  x
})
clinical_data <- clinical_data[, colSums(!is.na(clinical_data)) > 0, drop = FALSE]
for (var in names(clinical_data)) {
  if (var %in% c("record_id", "specimen_id")) next
  x <- clinical_data[[var]]
  non_missing <- !is.na(x)
  if (sum(non_missing) == 0) next
  suppressWarnings(x_numeric <- as.numeric(x[non_missing]))
  numeric_fraction <- mean(!is.na(x_numeric))
  if (numeric_fraction >= 0.90) {suppressWarnings(clinical_data[[var]] <- as.numeric(x))}
}
clinical_data <- clinical_data %>% mutate(scai_stage = toupper(consent_scai_stage_at_time_of_informed_consent), scai_stage = str_trim(scai_stage)) %>% filter(scai_stage %in% c("C", "D"))
# write.csv(clinical_data, file.path(output_dir, "CS_full_clinical_data_cleaned.csv"), row.names = FALSE)

get_phase <- function(x) {case_when(str_detect(x, "^admission") ~ "Admission", str_detect(x, "^shock_onset") ~ "Shock onset", str_detect(x, "^recovery") ~ "Recovery", str_detect(x, "^inpatient|^discharge|^post_discharge") ~ "Outcome", TRUE ~ "Other")}
get_category <- function(x) {
  case_when(str_detect(x, "_labs_") ~ "Laboratory", str_detect(x, "_rhc_") ~ "Right heart catheterization", str_detect(x, "_pac_") ~ "Pulmonary artery catheter", str_detect(x, "_echo_") ~ "Echocardiography",
    str_detect(x, regex("pressor|inotropic|norepinephrine|epinephrine|dobutamine|milrinone|vasopressin|mcs", ignore_case = TRUE)) ~ "Hemodynamic support", 
    str_detect(x, regex("etiologies|type_of_shock|category|physiology", ignore_case = TRUE)) ~ "Shock etiology", str_detect(x, regex("cardiac_arrest", ignore_case = TRUE)) ~ "Cardiac arrest",
    str_detect(x, regex("death|discharge|stroke|renal_failure|crrt|pneumonia|bacteremia|bleeding|ischemia|dialysis", ignore_case = TRUE)) ~ "Clinical outcomes",
    str_detect(x, regex("age|bmi|height|weight|body_surface", ignore_case = TRUE)) ~ "Demographics", TRUE ~ "Other")
}
clean_variable_label <- function(x) {
  x %>% str_remove("^admission_labs_") %>% str_remove("^admission_") %>% str_remove("^shock_onset_labs_") %>% str_remove("^shock_onset_rhc_") %>% str_remove("^shock_onset_pac_") %>%
    str_remove("^shock_onset_echo_") %>% str_remove("^shock_onset_") %>% str_remove("^recovery_pac_") %>% str_remove("^recovery_") %>% str_remove("^inpatient_") %>% str_remove("^discharge_") %>%
    str_remove("^post_discharge_") %>% str_replace_all("_mmhg", " mmHg") %>% str_replace_all("_mg_dl", " mg/dL") %>% str_replace_all("_mmol_l", " mmol/L") %>% str_replace_all("_g_dl", " g/dL") %>%
    str_replace_all("_iu_l", " IU/L") %>% str_replace_all("_ng_l", " ng/L") %>% str_replace_all("_ng_ml", " ng/mL") %>% str_replace_all("_l_min", " L/min") %>% str_replace_all("_", " ") %>%
    str_squish() %>% str_to_sentence()
}
classify_variable <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return("skip")
  if (is.numeric(x)) {values <- sort(unique(x)); if (all(values %in% c(0, 1))) {return("binary")}; return("numeric")}
  return("categorical")
}
exclude_patterns <- regex("date|time|other|specify|free_text|diagnosis|line_type|rhythm|medications_at_time", ignore_case = TRUE)
candidate_vars <- names(clinical_data) %>% setdiff(c("record_id", "specimen_id", "scai_stage")) %>% .[str_detect(., "^admission_|^shock_onset_|^recovery_|^inpatient_|^discharge_|^post_discharge_")] %>% .[!str_detect(., exclude_patterns)]
variable_dictionary <- tibble(variable_original = candidate_vars, phase = map_chr(candidate_vars, get_phase), category = map_chr(candidate_vars, get_category),
                              variable = map_chr(candidate_vars, clean_variable_label), variable_type = map_chr(candidate_vars, ~ classify_variable(clinical_data[[.x]]))) %>% filter(variable_type != "skip")
# write.csv(variable_dictionary, file.path(output_dir, "CS_variable_dictionary_for_analysis.csv"), row.names = FALSE)

analyze_numeric <- function(data, var) {
  df <- data %>% dplyr::select(scai_stage, value = all_of(var)) %>% mutate(value = as.numeric(value))
  C <- df$value[df$scai_stage == "C"]; D <- df$value[df$scai_stage == "D"]
  C <- C[!is.na(C)]; D <- D[!is.na(D)]
  if (length(C) < 3 || length(D) < 3) return(NULL)
  test <- tryCatch(wilcox.test(D, C, conf.int = TRUE, exact = FALSE), error = function(e) NULL)
  tibble(variable_original = var, n_C = length(C), n_D = length(D), missing_C = sum(data$scai_stage == "C") - length(C), missing_D = sum(data$scai_stage == "D") - length(D),
         mean_C = mean(C), sd_C = sd(C), median_C = median(C), q1_C = quantile(C, 0.25), q3_C = quantile(C, 0.75), mean_D = mean(D), sd_D = sd(D), median_D = median(D),  q1_D = quantile(D, 0.25), q3_D = quantile(D, 0.75),
         event_C = NA_real_, event_rate_C = NA_real_, event_D = NA_real_, event_rate_D = NA_real_, estimate = median(D) - median(C), estimate_type = "Median difference, D-C",
         ci_low = ifelse(is.null(test) || is.null(test$conf.int), NA, test$conf.int[1]), ci_high = ifelse(is.null(test) || is.null(test$conf.int), NA, test$conf.int[2]),
         statistical_test = "Wilcoxon rank-sum", p_value_numeric = ifelse(is.null(test), NA, test$p.value))
}
analyze_binary <- function(data, var) {
  df <- data %>% dplyr::select(scai_stage, value = all_of(var)) %>% mutate(value = as.numeric(value)) %>% filter(!is.na(value))
  if (!all(unique(df$value) %in% c(0, 1))) return(NULL)
  tab <- table(df$scai_stage, df$value)
  if (!all(c("C", "D") %in% rownames(tab))) return(NULL); if (!all(c("0", "1") %in% colnames(tab))) return(NULL)
  test <- tryCatch(fisher.test(tab), error = function(e) NULL)
  c_total <- sum(tab["C", ]); d_total <- sum(tab["D", ])
  c_events <- tab["C", "1"]; d_events <- tab["D", "1"]
  tibble(variable_original = var, n_C = c_total, n_D = d_total, missing_C = sum(data$scai_stage == "C") - c_total, missing_D = sum(data$scai_stage == "D") - d_total,
         mean_C = NA_real_, sd_C = NA_real_, median_C = NA_real_, q1_C = NA_real_, q3_C = NA_real_, mean_D = NA_real_, sd_D = NA_real_, median_D = NA_real_, q1_D = NA_real_, q3_D = NA_real_,
         event_C = c_events, event_rate_C = 100 * c_events / c_total, event_D = d_events, event_rate_D = 100 * d_events / d_total, estimate = ifelse(is.null(test), NA, as.numeric(test$estimate)), estimate_type = "Odds ratio, D vs C",
         ci_low = ifelse(is.null(test), NA, test$conf.int[1]), ci_high = ifelse(is.null(test), NA, test$conf.int[2]), statistical_test = "Fisher exact", p_value_numeric = ifelse(is.null(test), NA, test$p.value))
}
analyze_categorical <- function(data, var) {
  df <- data %>% dplyr::select(scai_stage, value = all_of(var)) %>% filter(!is.na(value))
  if (nrow(df) < 10) return(NULL)
  tab <- table(df$value, df$scai_stage)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  test <- tryCatch(fisher.test(tab), error = function(e) NULL)
  tibble(variable_original = var, n_C = sum(df$scai_stage == "C"), n_D = sum(df$scai_stage == "D"), missing_C = sum(data$scai_stage == "C") - sum(df$scai_stage == "C"), missing_D = sum(data$scai_stage == "D") - sum(df$scai_stage == "D"),
         mean_C = NA_real_, sd_C = NA_real_, median_C = NA_real_, q1_C = NA_real_, q3_C = NA_real_, mean_D = NA_real_, sd_D = NA_real_, median_D = NA_real_, q1_D = NA_real_, q3_D = NA_real_,
         event_C = NA_real_, event_rate_C = NA_real_, event_D = NA_real_, event_rate_D = NA_real_, estimate = NA_real_, estimate_type = "Global categorical comparison", ci_low = NA_real_, ci_high = NA_real_,
         statistical_test = "Fisher exact", p_value_numeric = ifelse(is.null(test), NA, test$p.value))
}

results <- list()
for (i in seq_len(nrow(variable_dictionary))) {
  current_var <- variable_dictionary$variable_original[i]
  current_type <- variable_dictionary$variable_type[i]
  result <- tryCatch({if (current_type == "numeric") {analyze_numeric(clinical_data, current_var)} 
    else if (current_type == "binary") {analyze_binary(clinical_data, current_var)} 
    else if (current_type == "categorical") {analyze_categorical(clinical_data, current_var)} 
    else {NULL}}, error = function(e) {message("Failed: ", current_var, " | ", e$message); NULL})
  if (!is.null(result)) {results[[length(results) + 1]] <- result}
}
analysis_results <- bind_rows(results)

phase_levels <- c("Admission", "Shock onset", "Recovery", "Outcome", "Other")
publication_table <- analysis_results %>% left_join(variable_dictionary, by = "variable_original") %>% group_by(phase) %>% mutate(p_adjusted_numeric = p.adjust(p_value_numeric, method = "BH")) %>%
  ungroup() %>% mutate(phase = factor(phase, levels = phase_levels), missingness = sprintf("C: %d; D: %d", missing_C, missing_D), c_mean_sd = ifelse(variable_type == "numeric", sprintf("%.2f ± %.2f", mean_C, sd_C), NA),
                       d_mean_sd = ifelse(variable_type == "numeric", sprintf("%.2f ± %.2f", mean_D, sd_D), NA), c_median_iqr = ifelse(variable_type == "numeric", sprintf("%.2f (%.2f–%.2f)", median_C, q1_C, q3_C), NA),
                       d_median_iqr = ifelse(variable_type == "numeric", sprintf("%.2f (%.2f–%.2f)", median_D, q1_D, q3_D), NA), c_events = ifelse(variable_type == "binary", sprintf("%d/%d (%.1f%%)", event_C, n_C, event_rate_C), NA),
                       d_events = ifelse(variable_type == "binary", sprintf("%d/%d (%.1f%%)", event_D, n_D, event_rate_D), NA), ci_95 = format_ci(ci_low, ci_high), p_value = format_p(p_value_numeric), p_adjusted = format_p(p_adjusted_numeric)) %>% arrange(phase, category, p_value_numeric)
numeric_table <- publication_table %>% filter(variable_type == "numeric") %>% transmute(Phase = phase, Category = category, Variable = variable, Original_Variable = variable_original, N_C = n_C, N_D = n_D,
                                                                                        Missingness = missingness, `C Mean ± SD` = c_mean_sd, `D Mean ± SD` = d_mean_sd, `C Median (IQR)` = c_median_iqr, `D Median (IQR)` = d_median_iqr,
                                                                                        Estimate_Type = estimate_type, Estimate = estimate, `95% CI` = ci_95, Statistical_Test = statistical_test, P_Value = p_value,
                                                                                        Adjusted_P_Value = p_adjusted, P_Value_Numeric = p_value_numeric, Adjusted_P_Value_Numeric = p_adjusted_numeric)

binary_table <- publication_table %>% filter(variable_type == "binary") %>%
  transmute(Phase = phase, Category = category, Variable = variable, Original_Variable = variable_original, N_C = n_C, N_D = n_D, Missingness = missingness, `C Events n/N (%)` = c_events, `D Events n/N (%)` = d_events,
            Estimate_Type = estimate_type, Estimate = estimate, `95% CI` = ci_95, Statistical_Test = statistical_test, P_Value = p_value, Adjusted_P_Value = p_adjusted, P_Value_Numeric = p_value_numeric, Adjusted_P_Value_Numeric = p_adjusted_numeric)
categorical_table <- publication_table %>% filter(variable_type == "categorical") %>%
  transmute(Phase = phase, Category = category, Variable = variable, Original_Variable = variable_original, N_C = n_C, N_D = n_D, Missingness = missingness, Statistical_Test = statistical_test, P_Value = p_value,
            Adjusted_P_Value = p_adjusted, P_Value_Numeric = p_value_numeric, Adjusted_P_Value_Numeric = p_adjusted_numeric)
all_variables_table <- publication_table %>% 
  transmute(Phase = phase, Category = category, Variable = variable, Original_Variable = variable_original, Variable_Type = variable_type, N_C = n_C, N_D = n_D, Missingness = missingness,
            `C Mean ± SD` = c_mean_sd, `D Mean ± SD` = d_mean_sd, `C Median (IQR)` = c_median_iqr, `D Median (IQR)` = d_median_iqr, `C Events n/N (%)` = c_events, `D Events n/N (%)` = d_events,
            Estimate_Type = estimate_type, Estimate = estimate, `95% CI` = ci_95, Statistical_Test = statistical_test, P_Value = p_value, Adjusted_P_Value = p_adjusted, P_Value_Numeric = p_value_numeric, Adjusted_P_Value_Numeric = p_adjusted_numeric)
# write.csv(all_variables_table, file.path(output_dir, "CS_all_variables_publication_table.csv"), row.names = FALSE)
# write.xlsx(all_variables_table, file.path(output_dir, "CS_all_variables_publication_table.xlsx"), overwrite = TRUE)
# write.xlsx(list("Numeric Variables" = numeric_table, "Binary Variables" = binary_table, "Categorical Variables" = categorical_table), file.path(output_dir, "CS_publication_tables_by_type.xlsx"), overwrite = TRUE)

# ============================================================
# Admission vs Shock analysis: Optional
# ============================================================
# clinical_data <- read_csv('cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv', show_col_types = FALSE)
# head(clinical_data)
# 
# sig_star <- function(p) {dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}
# adm_cols <- grep("^admission_", colnames(clinical_data), value = TRUE)
# shock_cols <- grep("^shock_onset_", colnames(clinical_data), value = TRUE)
# common_suffix <- intersect(sub("^admission_", "", adm_cols), sub("^shock_onset_", "", shock_cols))
# is_numeric_var <- function(suffix) is.numeric(clinical_data[[paste0("admission_", suffix)]])
# numeric_suffix <- Filter(is_numeric_var, common_suffix)
# length(numeric_suffix)
# 
# admshock_results <- purrr::map_dfr(numeric_suffix, function(v) {
#   adm <- clinical_data[[paste0("admission_", v)]]
#   shk <- clinical_data[[paste0("shock_onset_", v)]]
#   complete <- !is.na(adm) & !is.na(shk)
#   n <- sum(complete)
#   if (n < 5) return(NULL)
#   test <- wilcox.test(shk[complete], adm[complete], paired = TRUE)
#   data.frame(Variable_ = v, N_paired = n, Admission_median = median(adm[complete]), Shock_median = median(shk[complete]), Median_diff = median(shk[complete] - adm[complete]), P_Value = test$p.value)
# })
# 
# admshock_results$Adjusted_P_Value <- p.adjust(admshock_results$P_Value, method = "BH")
# admshock_results <- admshock_results[order(admshock_results$P_Value), ]
# admshock_results$sig <- sig_star(admshock_results$Adjusted_P_Value)
# print(admshock_results, row.names = FALSE)
# write.csv(admshock_results, file.path("cs_project/clinical_data/clinical_results/CS_Admission_vs_Shock_paired_results.csv"), row.names = FALSE)
# ============================================================
clinical_data <- read_csv('cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv', show_col_types = FALSE)
sig_star <- function(p) {dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}
stage_col <- "consent_scai_stage_at_time_of_informed_consent"
adm_cols <- grep("^admission_", colnames(clinical_data), value = TRUE)
shock_cols <- grep("^shock_onset_", colnames(clinical_data), value = TRUE)
common_suffix <- intersect(sub("^admission_", "", adm_cols), sub("^shock_onset_", "", shock_cols))
is_numeric_var <- function(suffix) is.numeric(clinical_data[[paste0("admission_", suffix)]])
numeric_suffix <- Filter(is_numeric_var, common_suffix)
length(numeric_suffix)
run_admshock_test <- function(data, vars) {
  purrr::map_dfr(vars, function(v) {
    adm <- data[[paste0("admission_", v)]]
    shk <- data[[paste0("shock_onset_", v)]]
    complete <- !is.na(adm) & !is.na(shk)
    n <- sum(complete)
    if (n < 2) return(NULL)
    test <- tryCatch(wilcox.test(shk[complete], adm[complete], paired = TRUE), error = function(e) NULL)
    if (is.null(test)) return(NULL)
    data.frame(Variable_ = v, N_paired = n, Admission_median = median(adm[complete]), Shock_median = median(shk[complete]), Median_diff = median(shk[complete] - adm[complete]), P_Value = test$p.value)
  })
}
admshock_C <- clinical_data %>% dplyr::filter(.data[[stage_col]] == "C") %>% run_admshock_test(numeric_suffix)
admshock_C$Adjusted_P_Value <- p.adjust(admshock_C$P_Value, method = "BH")
admshock_C <- admshock_C[order(admshock_C$P_Value), ]
admshock_C$sig <- sig_star(admshock_C$Adjusted_P_Value)
admshock_C$Stage <- "Stage C"
admshock_D <- clinical_data %>% dplyr::filter(.data[[stage_col]] == "D") %>% run_admshock_test(numeric_suffix)
admshock_D$Adjusted_P_Value <- p.adjust(admshock_D$P_Value, method = "BH")
admshock_D <- admshock_D[order(admshock_D$P_Value), ]
admshock_D$sig <- sig_star(admshock_D$Adjusted_P_Value)
admshock_D$Stage <- "Stage D"
print(admshock_C, row.names = FALSE)
print(admshock_D, row.names = FALSE)
admshock_by_stage <- dplyr::bind_rows(admshock_C, admshock_D)
write.csv(admshock_by_stage, file.path("cs_project/clinical_data/clinical_results/CS_Admission_vs_Shock_paired_results_byStage.csv"), row.names = FALSE)
saveRDS(list(admshock_C = admshock_C, admshock_D = admshock_D, admshock_by_stage = admshock_by_stage, clinical_data = clinical_data), file = "cs_project/clinical_data/clinical_results/CS_Admission_vs_Shock_byStage.rds")

# ============================================================
# Table 1: Demographics, Stage C vs Stage D
# ============================================================
clinical_data <- read_csv('cs_project/clinical_data/clinical_results/CS_full_clinical_data_cleaned.csv', show_col_types = FALSE) %>% dplyr::filter(scai_stage %in% c("C", "D"))
continuous_vars <- c(Age = "admission_patient_age_at_admission", BMI = "admission_bmi", "Height (cm)" = "admission_height_cm_at_admission", "Weight (kg)" = "admission_weight_kg_at_admission", "BSA" = "admission_body_surface_area_bsa")
categorical_vars <- c(Sex = "demographics_gender", Race = "demographics_race", Ethnicity = "demographics_ethnicity")
outcome_vars <- c("In-hospital death" = "inpatient_death", "Renal failure" = "inpatient_renal_failure",
                  "New CRRT/dialysis" = "inpatient_crrt", "Pneumonia" = "inpatient_pneumonia", "Bacteremia" = "inpatient_bacteremia")

fmt_num <- function(x) sprintf("%.1f", x)
n_c_total <- sum(clinical_data$scai_stage == "C"); n_d_total <- sum(clinical_data$scai_stage == "D"); n_total <- nrow(clinical_data)

continuous_summary <- purrr::imap_dfr(continuous_vars, function(col, label) {
  all_vals <- clinical_data[[col]]
  c_vals <- clinical_data[[col]][clinical_data$scai_stage == "C"]; d_vals <- clinical_data[[col]][clinical_data$scai_stage == "D"]
  test <- wilcox.test(d_vals, c_vals)
  data.frame(Variable = label, Overall = paste0(fmt_num(median(all_vals, na.rm = TRUE)), " (", fmt_num(quantile(all_vals, 0.25, na.rm = TRUE)), "\u2013", fmt_num(quantile(all_vals, 0.75, na.rm = TRUE)), ")"),
             Stage_C = paste0(fmt_num(median(c_vals, na.rm = TRUE)), " (", fmt_num(quantile(c_vals, 0.25, na.rm = TRUE)), "\u2013", fmt_num(quantile(c_vals, 0.75, na.rm = TRUE)), ")"),
             Stage_D = paste0(fmt_num(median(d_vals, na.rm = TRUE)), " (", fmt_num(quantile(d_vals, 0.25, na.rm = TRUE)), "\u2013", fmt_num(quantile(d_vals, 0.75, na.rm = TRUE)), ")"),
             P_Value = test$p.value)
})
categorical_summary <- purrr::imap_dfr(categorical_vars, function(col, label) {
  tbl_overall <- table(clinical_data[[col]])
  tbl <- table(clinical_data[[col]], clinical_data$scai_stage); test <- fisher.test(tbl)
  purrr::map_dfr(rownames(tbl), function(level) {
    n_all <- tbl_overall[level]; n_c <- tbl[level, "C"]; n_d <- tbl[level, "D"]
    data.frame(Variable = paste0(label, ": ", level), Overall = paste0(n_all, "/", n_total, " (", fmt_num(100 * n_all / n_total), "%)"),
               Stage_C = paste0(n_c, "/", n_c_total, " (", fmt_num(100 * n_c / n_c_total), "%)"), Stage_D = paste0(n_d, "/", n_d_total, " (", fmt_num(100 * n_d / n_d_total), "%)"), P_Value = test$p.value)})})
outcome_summary <- purrr::imap_dfr(outcome_vars, function(col, label) {
  vals_bin <- ifelse(!is.na(clinical_data[[col]]) & clinical_data[[col]] == 1, 1, 0)
  n_all <- sum(vals_bin)
  tbl <- table(vals_bin, clinical_data$scai_stage)
  test <- tryCatch(fisher.test(tbl), error = function(e) list(p.value = NA))
  n_c <- if ("1" %in% rownames(tbl)) tbl["1", "C"] else 0
  n_d <- if ("1" %in% rownames(tbl)) tbl["1", "D"] else 0
  data.frame(Variable = label, Overall = paste0(n_all, "/", n_total, " (", fmt_num(100 * n_all / n_total), "%)"), 
             Stage_C = paste0(n_c, "/", n_c_total, " (", fmt_num(100 * n_c / n_c_total), "%)"), Stage_D = paste0(n_d, "/", n_d_total, " (", fmt_num(100 * n_d / n_d_total), "%)"), P_Value = test$p.value)})
discharge_overall <- table(clinical_data$discharge_what_was_the_patients_outcome_at_discharge)
discharge_tbl <- table(clinical_data$discharge_what_was_the_patients_outcome_at_discharge, clinical_data$scai_stage)
discharge_test <- fisher.test(discharge_tbl)
discharge_summary <- purrr::map_dfr(rownames(discharge_tbl), function(level) {
  n_all <- discharge_overall[level]; n_c <- discharge_tbl[level, "C"]; n_d <- discharge_tbl[level, "D"]
  data.frame(Variable = paste0("Discharge outcome: ", level), Overall = paste0(n_all, "/", n_total, " (", fmt_num(100 * n_all / n_total), "%)"),
             Stage_C = paste0(n_c, "/", n_c_total, " (", fmt_num(100 * n_c / n_c_total), "%)"), Stage_D = paste0(n_d, "/", n_d_total, " (", fmt_num(100 * n_d / n_d_total), "%)"), P_Value = discharge_test$p.value)})

demo <- dplyr::bind_rows(continuous_summary, categorical_summary, outcome_summary, discharge_summary)
demo$P_Value_display <- ifelse(demo$P_Value < 0.001, "<0.001", sprintf("%.3f", demo$P_Value))
colnames(demo)[colnames(demo) == "Overall"] <- paste0("Overall (n=", n_total, ")")
colnames(demo)[colnames(demo) %in% c("Stage_C", "Stage_D")] <- c(paste0("Stage C (n=", n_c_total, ")"), paste0("Stage D (n=", n_d_total, ")"))
write.csv(demo, file.path("cs_project/clinical_data/clinical_results/Table1_Demographics.csv"), row.names = FALSE)
# Of 42 enrolled patients, bulk RNA-sequencing was performed on 40 and plasma cytokine profiling on 41, reflecting specimen availability. 
# Demographic characteristics of the full enrolled cohort are shown in Table 1.



# # Sampe mismatch check: Optional
# library(dplyr)
# normalize_sr_id <- function(sr) {
#   digits <- gsub("[^0-9]", "", sr)
#   year <- substr(digits, 1, 2)
#   serial <- as.character(as.integer(substr(digits, 3, nchar(digits))))
#   paste0("SR", year, serial)
# }
# cytokine_data <- readxl::read_excel("cs_project/cytokine_data/CScytokineData.xlsx") %>% dplyr::rename(SR_specimen_ID = `SR specimen ID`) %>% dplyr::mutate(SR_norm = normalize_sr_id(SR_specimen_ID))
# metadata_bulk <- read.csv("cs_project/bulk_data/metadata_bulk.csv") %>% dplyr::mutate(SR_norm = normalize_sr_id(SR.specimen.ID))
# only_in_cytokine <- setdiff(cytokine_data$SR_norm, metadata_bulk$SR_norm)
# only_in_bulk <- setdiff(metadata_bulk$SR_norm, cytokine_data$SR_norm)
# cat("Only in cytokine data:", only_in_cytokine, "\n")
# cat("Only in bulk metadata:", only_in_bulk, "\n")
# cat("Common to both:", length(intersect(cytokine_data$SR_norm, metadata_bulk$SR_norm)), "\n")
# merged <- inner_join(cytokine_data, metadata_bulk, by = "SR_norm")
# dim(merged)