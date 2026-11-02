setwd(Sys.getenv("CS_ROOT", "/mnt/home/koiralas"))
# ============================================================
# Bulk RNA-seq analysis: SCAI D vs C
# ============================================================

library(stringr)
library(tibble)
library(ggplot2)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(dplyr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(writexl)
library(tidyverse)
library(msigdbr)
library(AnnotationDbi)
library(KEGGREST)
library(reactome.db)
library(GSVA)


# 1. Load bulk RNA-seq data
load("cs_project/bulk_data/count_data.RData")
ls()
dim(log2.fpkm.matrix)
dim(log2.read.count.matrix)
dim(log2.tpm.matrix)
head(rownames(log2.tpm.matrix))
head(colnames(log2.tpm.matrix))
head(rownames(log2.read.count.matrix))
head(colnames(log2.read.count.matrix))

metadata <- read.csv("cs_project/bulk_data/metadata_bulk.csv")
head(metadata)

# 2. Choose expression matrix for exploratory analysis
expr <- log2.tpm.matrix
rna_samples <- colnames(expr)
# head(rna_samples)

# 3. Clean RNA metadata
metadata <- metadata %>% mutate(SR.specimen.ID = str_trim(SR.specimen.ID), scai_stage_bulk = toupper(str_trim(SCAI.Stage)))
# table(metadata$scai_stage_bulk, useNA = "ifany")
# sum(colnames(expr) %in% metadata$SR.specimen.ID)
# setdiff(colnames(expr), metadata$SR.specimen.ID)
# setdiff(metadata$SR.specimen.ID, colnames(expr))

rna_meta <- metadata %>% filter(scai_stage_bulk %in% c("C", "D"))
expr_cd <- expr[, rna_meta$SR.specimen.ID]
table(rna_meta$scai_stage_bulk)
dim(expr_cd)

# 4. Exploratory PCA using log2 TPM
gene_var <- apply(expr_cd, 1, var, na.rm = TRUE)
expr_filt <- expr_cd[gene_var > quantile(gene_var, 0.25), ]
pca <- prcomp(t(expr_filt), scale. = TRUE)
pca_df <- data.frame(sample = rna_meta$SR.specimen.ID, scai = rna_meta$scai_stage_bulk, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = scai)) + geom_point(size = 3) + theme_minimal() + labs(title = "RNA PCA: SCAI C vs D")
# p_pca

# 5. Reconstruct approximate raw counts from log2(count + 1)
# range(log2.read.count.matrix)
counts <- round(2^log2.read.count.matrix - 1)
# range(counts)
counts[counts < 0] <- 0
counts[is.na(counts)] <- 0
# range(counts)

summary(as.vector(counts))
all(counts == round(counts))

counts <- counts[, rna_meta$SR.specimen.ID]
stopifnot(all(colnames(counts) == rna_meta$SR.specimen.ID))

# 6. DESeq2 differential expression: D vs C
coldata <- rna_meta %>% dplyr::select(SR.specimen.ID, scai_stage_bulk) %>% column_to_rownames("SR.specimen.ID")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ scai_stage_bulk)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("scai_stage_bulk", "D", "C"))
res <- as.data.frame(res)
res$gene <- rownames(res)
resLFC <- lfcShrink(dds, coef = "scai_stage_bulk_D_vs_C", type = "apeglm")
resLFC <- as.data.frame(resLFC)
resLFC$gene <- rownames(resLFC)
# head(resLFC, 20)
# summary(resLFC$padj)

# 7. Add gene annotation
mapping <- read.csv("cs_project/gene_info_biomart.csv")
gene_map <- mapping %>% dplyr::select(gene_id, gene_name, gene_description) %>% distinct(gene_id, .keep_all = TRUE)
resLFC_annot <- resLFC %>% left_join(gene_map, by = c("gene" = "gene_id")) %>% relocate(gene, gene_name, gene_description)
# head(resLFC_annot)
# hist(resLFC_annot$log2FoldChange)
# hist(resLFC_annot$lfcSE)

dir.create("cs_project/bulk_data/results", recursive = TRUE, showWarnings = FALSE)
write.csv(resLFC_annot, "cs_project/bulk_data/results/dge_results_annotated.csv", row.names = FALSE)

# 8. Define DE genes: LFC cutoff - 1 or 0.5?
# de_sig <- resLFC_annot %>% filter(!is.na(padj)) %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1)
# head(de_sig)
# range(de_sig$log2FoldChange)
# hist(de_sig$log2FoldChange)
# hist(de_sig$lfcSE)
# write.csv(de_sig, "cs_project/bulk_data/results/dge_results_significant_1.csv", row.names = FALSE)

de_sig <- resLFC_annot %>% filter(!is.na(padj)) %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 0.5)
head(de_sig)
range(de_sig$log2FoldChange)
hist(de_sig$log2FoldChange)
hist(de_sig$lfcSE)
write.csv(de_sig, "cs_project/bulk_data/results/dge_results_significant_0.5.csv", row.names = FALSE)

# 9. Clinically anchored gene-axis analysis
# Based on clinical findings:
# - hyperglycemia / metabolic stress
# - hypoalbuminemia / capillary leak
# - vasoplegia / low SVR
# - pressor escalation
# - inflammation / innate immunity
# - endothelial dysfunction
# - hypoxia / mitochondrial stress
# - coagulation / complement
clinical_gene_axes <- list(
  metabolic_stress_hyperglycemia = c("SLC2A1", "SLC2A3", "HK2", "PFKFB3", "LDHA", "PDK1", "HIF1A", "DDIT3", "XBP1", "ATF4", "FOS", "JUN", "NR4A1", "TXNIP", "G6PD", "ENO1", "ALDOA", "PKM"),
  inflammation_capillary_leak = c("IL6", "IL1B", "TNF", "NFKBIA", "NFKB1", "RELA", "CXCL8", "CCL2", "CCL3", "CCL4", "SOCS3", "S100A8", "S100A9", "PTGS2", "MMP9", "LCN2", "OSM", "IL1RN", "STAT3"),
  endothelial_activation_permeability = c("VWF", "ICAM1", "VCAM1", "SELE", "SELP", "PECAM1", "CDH5", "CLDN5", "KDR", "FLT1", "TEK", "ANGPT1", "ANGPT2", "VEGFA", "SERPINE1", "THBD", "ENG", "ESAM"),
  vasoplegia_vascular_tone = c("NOS2", "NOS3", "ARG1", "ARG2", "GCH1", "EDN1", "ADM", "CALCRL", "RAMP2", "RAMP3", "PTGS2", "ADRA1A", "ADRA1B", "ADRA1D", "ADRB1", "ADRB2", "AVPR1A", "AVPR1B", "AGTR1", "AGTR2"),
  neutrophil_activation = c("S100A8", "S100A9", "MPO", "ELANE", "MMP8", "MMP9", "CEACAM8", "FCGR3B", "CXCR1", "CXCR2", "OLFM4", "LCN2", "DEFA1", "DEFA3", "AZU1", "CTSG"),
  monocyte_macrophage_activation = c("CD14", "FCN1", "LST1", "VCAN", "CCR2", "CCL2", "CD68", "MSR1", "ITGAM", "TYROBP", "FCGR1A", "FCGR2A", "LYZ", "S100A12", "TREM1"),
  hypoxia_mitochondrial_stress = c("HIF1A", "EGLN1", "EGLN3", "VEGFA", "BNIP3", "NDUFA1", "NDUFB8", "COX4I1", "ATP5F1A", "SOD2", "TXNIP", "UCP2", "PGK1", "ENO1", "CA9", "ADM"),
  coagulation_complement = c("C3", "C5", "C5AR1", "CFD", "CFB", "F3", "SERPINE1", "THBD", "PLAU", "PLAUR", "ITGA2B", "PF4", "PPBP", "GP1BA", "GP9", "CLU", "C1QA", "C1QB", "C1QC"))

axis_gene_table <- bind_rows(lapply(names(clinical_gene_axes), function(axis) {
  data.frame(clinical_axis = axis, gene_name = clinical_gene_axes[[axis]], stringsAsFactors = FALSE)
  })
)

# 10. Search DE results within clinical axes
de_annot <- de_sig %>% mutate(direction = case_when(log2FoldChange > 0.5 ~ "Higher in D", log2FoldChange < 0.5 ~ "Lower in D"), fdr_sig = ifelse(!is.na(padj) & padj < 0.05, TRUE, FALSE))
# table(de_annot$direction)
targeted_de_hits <- de_annot %>% inner_join(axis_gene_table, by = "gene_name") %>% arrange(clinical_axis, padj)
write.csv(targeted_de_hits, "cs_project/bulk_data/results/targeted_DE_genes_by_clinical_axis.csv", row.names = FALSE)

# # 11. Summarize clinical axes
# axis_summary <- targeted_de_hits %>% group_by(clinical_axis) %>% 
#   summarise(n_genes_detected = n(), n_fdr_p05 = sum(fdr_sig, na.rm = TRUE), 
#             n_up_in_D = sum(log2FoldChange > 0, na.rm = TRUE), n_down_in_D = sum(log2FoldChange < 0, na.rm = TRUE),
#             median_log2FC = median(log2FoldChange, na.rm = TRUE), top_nominal_genes = paste(head(gene_name[order(pvalue)], 10), collapse = ", "), .groups = "drop") %>%
#   arrange(desc(abs(median_log2FC)))
# axis_summary
# write.csv(axis_summary, "cs_project/bulk_data/results/clinical_axis_DE_summary.csv", row.names = FALSE)

# 12. Heatmap of clinically anchored genes
vsd <- vst(dds, blind = TRUE)
expr_vst <- assay(vsd)

target_genes_ensembl <- targeted_de_hits %>% pull(gene) %>% unique()
target_genes_ensembl <- intersect(target_genes_ensembl, rownames(expr_vst))

if (length(target_genes_ensembl) >= 2) {
  expr_target <- expr_vst[target_genes_ensembl, ]
  gene_symbols <- targeted_de_hits %>% dplyr::select(gene, gene_name) %>% distinct()
  rownames(expr_target) <- gene_symbols$gene_name[match(rownames(expr_target), gene_symbols$gene)]
  annotation_col <- data.frame(SCAI = rna_meta$scai_stage_bulk)
  rownames(annotation_col) <- rna_meta$SR.specimen.ID
  annotation_col <- annotation_col[colnames(expr_target), , drop = FALSE]
  ord <- order(annotation_col$SCAI)
  expr_target <- expr_target[, ord]
  annotation_col <- annotation_col[ord, , drop = FALSE]
  scai_counts <- table(annotation_col$SCAI)
  gaps_col <- cumsum(scai_counts)[-length(scai_counts)]
  
  pheatmap(expr_target, scale = "row", annotation_col = annotation_col, cluster_cols = FALSE, gaps_col = gaps_col, show_colnames = TRUE, fontsize_row = 8,
           main = "Clinically anchored DE genes: SCAI D vs C", filename = "cs_project/bulk_data/results/clinically_anchored_DE_heatmap.png", width = 7, height = 3)
}

# 13. Volcano plot with targeted axis genes highlighted
volcano_df <- resLFC_annot %>% 
  mutate(neg_log10_p = -log10(padj), direction = case_when(padj < 0.05 & log2FoldChange > 0.5  ~ "Up in D", padj < 0.05 & log2FoldChange < -0.5 ~ "Down in D", TRUE ~ "Not significant"),
         is_target_gene = gene_name %in% targeted_de_hits$gene_name, label_gene = ifelse(is_target_gene, gene_name, NA_character_))
p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_p)) + geom_point(aes(color = direction), alpha = 0.5, size = 1.5) +
  geom_point(data = volcano_df %>% filter(is_target_gene), color = "black", size = 2.5, alpha = 1) +
  geom_text_repel(data = volcano_df %>% filter(is_target_gene), aes(label = label_gene), color = "black", fontface = "bold", size = 3.5, box.padding = 1, point.padding = 0.5,
                  max.overlaps = Inf, segment.color = "black", segment.size = 0.4, arrow = arrow(length = unit(0.018, "npc"), type = "closed")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Down in D" = "steelblue", "Up in D" = "darkred", "Not significant" = "lightgrey")) + theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")) + labs(x = "log2 fold change: D vs C", y = "-log10 adjusted p-value")
p_volcano
ggsave("cs_project/bulk_data/results/volcano_clinical_axis_genes.png", plot = p_volcano, width = 7, height = 6, dpi = 300)

# 14. Pathway enrichment
head(resLFC_annot)
head(de_annot)
table(de_annot$direction)
range(de_annot$log2FoldChange); mean(de_annot$log2FoldChange)
range(de_annot$padj); mean(de_annot$padj)

genes <- list(Up = de_annot %>% filter(log2FoldChange > 0.5, padj < 0.05, gene_name != "") %>% pull(gene_name) %>% unique(),
              Down = de_annot %>% filter(log2FoldChange < -0.5, padj < 0.05, gene_name != "") %>% pull(gene_name) %>% unique())
universe <- resLFC_annot %>% filter(baseMean > 10, !is.na(padj)) %>% pull(gene_name) %>% unique()
to_entrez <- function(x) {bitr(unique(x), "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db) %>% pull(ENTREZID) %>% unique()}

genes_id <- lapply(genes, to_entrez)
universe_id <- to_entrez(universe)
go_bp <- compareCluster(genes_id, fun = "enrichGO", universe = universe_id, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
reactome <- compareCluster(genes_id, fun = "enrichPathway", universe = universe_id, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 1, readable = TRUE)
kegg <- compareCluster(genes_id, fun = "enrichKEGG", universe = universe_id, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1)

ora_df <- bind_rows(as.data.frame(go_bp) %>% mutate(Database = "GO_BP"), as.data.frame(reactome) %>% mutate(Database = "Reactome"), as.data.frame(kegg) %>% mutate(Database = "KEGG")) %>%
  mutate(Direction = Cluster, Sig_FDR_0.05 = p.adjust < 0.05, Sig_nominal_0.05 = pvalue < 0.05, logFDR = -log10(p.adjust)) %>%
  dplyr::select(Database, Direction, Sig_FDR_0.05, Sig_nominal_0.05, ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, qvalue, logFDR, geneID) %>%
  arrange(Database, Direction, p.adjust)

plot_df <- ora_df %>% mutate(plot_p = ifelse(Sig_FDR_0.05, p.adjust, pvalue), Evidence = ifelse(Sig_FDR_0.05, "FDR < 0.05", "Nominal only"),
                                 score = -log10(plot_p), signed_score = ifelse(Direction == "Down", -score, score), label = str_wrap(Description, 42)) %>%
  filter(Sig_FDR_0.05 | Sig_nominal_0.05) %>% group_by(Database, Direction) %>% slice_min(plot_p, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(label2 = paste0(Database, " | ", label), label2 = fct_reorder(label2, signed_score)) %>%
  mutate(Group = case_when(Direction == "Up"   & Evidence == "FDR < 0.05" ~ "Up | FDR sig", Direction == "Up"   & Evidence == "Nominal only" ~ "Up | nominal",
                           Direction == "Down" & Evidence == "FDR < 0.05" ~ "Down | FDR sig", TRUE ~ "Down | nominal"))

p_ora <- ggplot(plot_df, aes(x = signed_score, y = label2, fill = Group)) + geom_col(width = 0.78) + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  facet_grid(Database ~ ., scales = "free_y", space = "free_y") + scale_x_continuous(labels = abs) +
  scale_fill_manual(values = c("Up | FDR sig" = "darkred", "Up | nominal" = "#D55E00", "Down | FDR sig" = "steelblue", "Down | nominal" = "#92C5DE")) +
  labs(title = "Pathway enrichment ORA", subtitle = "Right = Upregulated genes; Left = Downregulated genes", x = expression(-log[10]~"(adjusted p or nominal p)"), y = NULL, fill = NULL) +
  theme_classic(base_size = 12) + theme(legend.position = "top", strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold"),
                                        axis.text.y = element_text(size = 8.5), plot.title = element_text(face = "bold"))
p_ora
ggsave("cs_project/bulk_data/results/pathway_enrichment_ORA.pdf", p_ora, width = 11, height = 8)
ggsave("cs_project/bulk_data/results/pathway_enrichment_ORA.png", p_ora, width = 11, height = 8, dpi = 300)
write_xlsx(list(pathway_enrichment_all = ora_df), "cs_project/bulk_data/results/pathway_Enrichment_ORA.xlsx")

# 14. Extended Pathway enrichment using ssGSEA
pathway_file <- "cs_project/bulk_data/results/pathway_Enrichment_ORA.xlsx"
pathway <- readxl::read_excel(pathway_file, sheet = "pathway_enrichment_all")
cs_theme_lookup <- tibble::tribble(
  ~Description,                                                                              ~Theme,
  "Hemostasis",                                                                               "Coagulation & platelet",
  "Platelet activation, signaling and aggregation",                                           "Coagulation & platelet",
  "Platelet degranulation",                                                                   "Coagulation & platelet",
  "Formation of Fibrin Clot (Clotting Cascade)",                                               "Coagulation & platelet",
  "Complement and coagulation cascades",                                                       "Coagulation & platelet",
  "platelet activation",                                                                       "Coagulation & platelet",
  "blood coagulation",                                                                         "Coagulation & platelet",
  "hemostasis",                                                                                "Coagulation & platelet",
  "complement activation",                                                                     "Complement",
  "regulation of complement activation",                                                       "Complement",
  "complement receptor mediated signaling pathway",                                            "Complement",
  "pattern recognition receptor signaling pathway",                                            "PRR & innate immune",
  "toll-like receptor signaling pathway",                                                      "PRR & innate immune",
  "innate immune response-activating signaling pathway",                                       "PRR & innate immune",
  "regulation of innate immune response",                                                      "PRR & innate immune",
  "Neutrophil extracellular trap formation",                                                   "Neutrophil & NETs",
  "neutrophil mediated cytotoxicity",                                                          "Neutrophil & NETs",
  "neutrophil chemotaxis",                                                                     "Neutrophil & NETs",
  "humoral immune response",                                                                   "Immunoglobulin & humoral",
  "immunoglobulin production",                                                                 "Immunoglobulin & humoral",
  "immunoglobulin mediated immune response",                                                   "Immunoglobulin & humoral",
  "interferon-mediated signaling pathway",                                                     "Interferon",
  "response to type II interferon",                                                            "Interferon",
  "mitochondrion organization",                                                                "Mitochondrial & apoptosis",
  "regulation of mitochondrial membrane permeability",                                         "Mitochondrial & apoptosis",
  "apoptotic mitochondrial changes",                                                            "Mitochondrial & apoptosis",
  "release of cytochrome c from mitochondria",                                                 "Mitochondrial & apoptosis",
  "Cell surface interactions at the vascular wall",                                            "Vascular, endothelial & ECM",
  "ECM-receptor interaction",                                                                  "Vascular, endothelial & ECM",
  "leukocyte adhesion to vascular endothelial cell",                                           "Vascular, endothelial & ECM",
  "vascular associated smooth muscle cell proliferation",                                      "Vascular, endothelial & ECM",
  "cellular response to vascular endothelial growth factor stimulus",                          "Vascular, endothelial & ECM",
  "wound healing",                                                                             "Wound healing",
  "response to wounding",                                                                      "Wound healing",
  "regulation of wound healing",                                                               "Wound healing",
  "RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function",      "Megakaryocyte & TF programs",
  "megakaryocyte differentiation",                                                             "Megakaryocyte & TF programs"
)

target_pathways <- pathway %>% dplyr::inner_join(cs_theme_lookup, by = "Description") %>% dplyr::distinct(ID, Description, Database, Theme)
msig_gobp     <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msig_kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
msig_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
msig_all <- dplyr::bind_rows(msig_gobp, msig_kegg, msig_reactome) %>% dplyr::distinct(gs_exact_source, gene_symbol, .keep_all = TRUE)
matched <- target_pathways %>% dplyr::inner_join(msig_all, by = c("ID" = "gs_exact_source"))
known_gene_sets <- split(matched$gene_symbol, matched$Description)
unmatched <- target_pathways %>% dplyr::anti_join(matched, by = "ID")
get_go_genes <- function(go_id) {AnnotationDbi::select(org.Hs.eg.db, keys = go_id, keytype = "GOALL", columns = "SYMBOL") %>% dplyr::pull(SYMBOL) %>% unique()}
get_kegg_genes <- function(kegg_id) {
  genes_raw <- KEGGREST::keggGet(kegg_id)[[1]]$GENE
  entrez <- genes_raw[seq(1, length(genes_raw), by = 2)]
  AnnotationDbi::select(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = "SYMBOL") %>% dplyr::pull(SYMBOL) %>% unique()
}
get_reactome_genes_by_name <- function(description) {
  pathname <- paste0("Homo sapiens: ", description)
  entrez <- AnnotationDbi::select(reactome.db, keys = pathname, keytype = "PATHNAME", columns = "ENTREZID") %>% dplyr::pull(ENTREZID) %>% unique()
  AnnotationDbi::select(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = "SYMBOL") %>% dplyr::pull(SYMBOL) %>% unique()
}
fallback_sets <- unmatched %>% purrr::pmap(function(ID, Description, Database, Theme) {
  genes <- switch(Database, "GO_BP" = get_go_genes(ID), "KEGG" = get_kegg_genes(ID), "Reactome" = get_reactome_genes_by_name(Description), character(0))
  setNames(list(genes), Description)
  }) %>% purrr::flatten()
known_gene_sets <- c(known_gene_sets, fallback_sets)
# gene_set_sizes <- tibble::tibble(Description = names(known_gene_sets), n_genes = purrr::map_int(known_gene_sets, length))
# print(gene_set_sizes)
# stopifnot(all(target_pathways$Description %in% names(known_gene_sets)))
# stopifnot(all(gene_set_sizes$n_genes > 0))

gene_map <- gene_map %>% dplyr::filter(!is.na(gene_name), gene_name != "")
expr_symbol <- as.data.frame(expr_cd) %>% tibble::rownames_to_column("gene_id") %>% dplyr::mutate(gene_id = stringr::str_remove(gene_id, "\\..*")) %>% 
  dplyr::inner_join(gene_map %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% dplyr::mutate(mean_expr = rowMeans(dplyr::across(where(is.numeric)))) %>% dplyr::group_by(gene_name) %>%
  dplyr::slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::select(-gene_id, -mean_expr) %>% tibble::column_to_rownames("gene_name") %>% as.matrix()
# dim(expr_symbol)

gsva_par <- GSVA::ssgseaParam(expr = expr_symbol, geneSets = known_gene_sets, minSize = 5, maxSize = 1000)
ssgsea_known <- GSVA::gsva(gsva_par)
scores_long_known <- as.data.frame(ssgsea_known) %>% tibble::rownames_to_column("pathway") %>% tidyr::pivot_longer(-pathway, names_to = "SR.specimen.ID", values_to = "es") %>%
  dplyr::left_join(rna_meta %>% dplyr::select(SR.specimen.ID, scai_stage_bulk), by = "SR.specimen.ID") %>% dplyr::left_join(target_pathways %>% dplyr::select(Description, Theme), by = c("pathway" = "Description"))
stat_res_known <- scores_long_known %>% dplyr::group_by(pathway) %>% dplyr::summarise(p_wilcox = wilcox.test(es ~ scai_stage_bulk)$p.value, mean_C = mean(es[scai_stage_bulk == "C"]), mean_D = mean(es[scai_stage_bulk == "D"]), .groups = "drop") %>%
  dplyr::mutate(p_adj = p.adjust(p_wilcox, method = "BH"), Direction = dplyr::if_else(mean_D > mean_C, "Higher in Stage D", "Higher in Stage C")) %>% dplyr::arrange(p_adj)
print(stat_res_known)
# saveRDS(list(stat_res = stat_res_known, scores_long = scores_long_known), file = "cs_project/bulk_data/results/ssgsea_known_pathways_tpm_results.rds")
# saveRDS(list(stat_res = stat_res_known, scores_long = scores_long_known, target_pathways = target_pathways, cs_theme_lookup = cs_theme_lookup), file = "cs_project/bulk_data/results/ssgsea_known_pathways_tpm_results.rds")

# # 15a. Custom pathway set enrichment - using DE genes
# library(tidyverse)
# library(clusterProfiler)
# library(enrichplot)
# 
# de_file        <- "cs_project/bulk_data/results/dge_results_annotated.csv"
# pathway_gmt_dir <- "cs_project/bulk_data/pathway sets/"
# 
# de_ranked <- read.csv(de_file) %>% dplyr::filter(!is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0, !is.na(gene_name), gene_name != "") %>% dplyr::mutate(stat_calc = log2FoldChange / lfcSE) %>%
#   dplyr::group_by(gene_name) %>% dplyr::slice_max(order_by = abs(stat_calc), n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::arrange(desc(stat_calc))
# gene_list <- setNames(de_ranked$stat_calc, de_ranked$gene_name)
# gene_list <- sort(gene_list, decreasing = TRUE)
# term2gene <- list.files(pathway_gmt_dir, pattern = "\\.gmt$", full.names = TRUE) %>% lapply(clusterProfiler::read.gmt) %>% dplyr::bind_rows()
# term2gene <- term2gene %>% dplyr::mutate(term = stringr::str_replace(term, "\\s*\\(?([A-Z]{2,3}):(\\d+)\\)?\\s*$", " (\\1:\\2)"), term = stringr::str_squish(term))
# 
# set.seed(42)
# gsea_res <- clusterProfiler::GSEA(gene_list, TERM2GENE = term2gene, minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, eps = 0)
# gsea_df <- as.data.frame(gsea_res) %>% dplyr::mutate(Sig_FDR_0.05 = p.adjust < 0.05, Direction = dplyr::if_else(NES > 0, "Higher in Stage D", "Higher in Stage C"), logFDR = -log10(p.adjust)) %>% dplyr::arrange(p.adjust)
# print(gsea_df %>% dplyr::select(ID, NES, pvalue, p.adjust, Direction, Sig_FDR_0.05))
# write.csv(gsea_df, file.path("cs_project/bulk_data/results/gsea_custom_metabolic.csv"), row.names = FALSE)
# saveRDS(gsea_res, file.path("cs_project/bulk_data/results/gsea_res_custom_metabolic.rds"))

# 15b. Custom pathway set enrichment - II using TPM
library(tidyverse)
library(GSVA)
library(clusterProfiler)

pathway_gmt_dir <- "cs_project/bulk_data/pathway sets/"
term2gene <- list.files(pathway_gmt_dir, pattern = "\\.gmt$", full.names = TRUE) %>% lapply(clusterProfiler::read.gmt) %>% dplyr::bind_rows() %>% 
  dplyr::mutate(term = stringr::str_replace(term, "\\s*\\(?([A-Z]{2,3}):(\\d+)\\)?\\s*$", " (\\1:\\2)"), term = stringr::str_squish(term))
gene_sets <- split(term2gene$gene, term2gene$term)
gsva_par <- GSVA::ssgseaParam(expr = expr_symbol, geneSets = gene_sets, minSize = 5, maxSize = 1000)
ssgsea_scores <- GSVA::gsva(gsva_par)   # pathways x samples matrix
scores_long <- as.data.frame(ssgsea_scores) %>% tibble::rownames_to_column("pathway") %>% tidyr::pivot_longer(-pathway, names_to = "SR.specimen.ID", values_to = "es") %>%
  dplyr::left_join(rna_meta %>% dplyr::select(SR.specimen.ID, scai_stage_bulk), by = "SR.specimen.ID")
stat_res <- scores_long %>% dplyr::group_by(pathway) %>% dplyr::summarise(p_wilcox = wilcox.test(es ~ scai_stage_bulk)$p.value, mean_C = mean(es[scai_stage_bulk == "C"]), mean_D = mean(es[scai_stage_bulk == "D"]), .groups = "drop") %>%
  dplyr::mutate(p_adj = p.adjust(p_wilcox, method = "BH"), Direction = dplyr::if_else(mean_D > mean_C, "Higher in Stage D", "Higher in Stage C")) %>% dplyr::arrange(p_adj)
saveRDS(list(stat_res = stat_res, scores_long = scores_long), file = "cs_project/bulk_data/results/ssgsea_custom_pathway_tpm_results.rds")

# 16. Co-Expression Analysis
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

load("cs_project/bulk_data/count_data.RData")
metadata <- read.csv("cs_project/bulk_data/metadata_bulk.csv")
expr_raw <- log2.tpm.matrix
colnames(expr_raw) <- sub("\\..*", "", colnames(expr_raw))
metadata <- metadata[match(colnames(expr_raw), metadata$SR.specimen.ID), ]
keep_samples <- metadata$SCAI.Stage %in% c("C", "D")
expr_raw <- expr_raw[, keep_samples]; metadata <- metadata[keep_samples, ]

expr <- expr_raw[rowSums(expr_raw > 1) >= 0.8 * ncol(expr_raw), ]
libSize_scaled <- scale(colSums(log2.read.count.matrix[, colnames(expr)]))[, 1]
expr_corrected <- t(apply(expr, 1, function(gene) resid(lm(gene ~ libSize_scaled))))
dimnames(expr_corrected) <- dimnames(expr)
expr_var <- expr_corrected[names(sort(apply(expr_corrected, 1, var), decreasing = TRUE))[1:5000], ]

# Module detection
softPower <- 14
net <- blockwiseModules(t(expr_var), power = softPower, networkType = "signed", TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.05,
                        numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0)
table(net$colors)
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
sample_order <- colnames(expr_var)
rownames(MEs) <- sample_order

# Module-Trait/Stage relationship
traitData <- data.frame(stage = ifelse(metadata$SCAI.Stage[match(sample_order, metadata$SR.specimen.ID)] == "D", 1, 0))
rownames(traitData) <- sample_order

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(MEs))
stage_results <- data.frame(module = rownames(moduleTraitCor), cor = moduleTraitCor[, "stage"], pval = moduleTraitP[, "stage"])
stage_results <- stage_results[order(stage_results$pval), ]
stage_results$FDR <- p.adjust(stage_results$pval, method = "BH")
print(stage_results, row.names = FALSE)

top_module <- stage_results$module[1]
top_label  <- as.numeric(sub("ME", "", top_module))
top_module_genes <- rownames(expr_var)[net$colors == top_label]
length(top_module_genes)
wgcna_results <- list(net = net, moduleColors = moduleColors, MEs = MEs, expr_var = expr_var, softPower = softPower, stage_results = stage_results, 
                      top_module = top_module, top_module_genes = top_module_genes, sample_order = sample_order, metadata = metadata)
saveRDS(wgcna_results, file.path("cs_project/bulk_data/results/wgcna_results.rds"))



# # Export RSEM-derived expression matrices to plain-text files for GEO submission: Optional
# rdata_file <- "cs_project/bulk_data/count_data.RData"
# meta_file  <- "cs_project/bulk_data/metadata_bulk.csv"
# out_dir    <- "cs_project/geo_upload"
# 
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---- Load matrices ----------------------------------------------------------
# load(rdata_file)   # loads: log2.fpkm.matrix, log2.read.count.matrix, log2.tpm.matrix
# 
# # ---- Select the samples being submitted (SCAI Stage C or D only) ------------
# meta <- read.csv(meta_file, stringsAsFactors = FALSE)
# meta <- meta[meta$SCAI.Stage %in% c("C", "D"), ]
# keep_ids <- meta$SR.specimen.ID
# cat("Samples to export:", length(keep_ids), "\n")   # expect 39
# 
# mats <- list(
#   log2_read_count_matrix = log2.read.count.matrix,
#   log2_fpkm_matrix       = log2.fpkm.matrix,
#   log2_tpm_matrix        = log2.tpm.matrix
# )
# 
# # ---- Sanity  ----------------------------------------------------------
# for (nm in names(mats)) {
#   missing <- setdiff(keep_ids, colnames(mats[[nm]]))
#   if (length(missing) > 0) {
#     stop(nm, ": these samples are missing from the matrix: ",
#          paste(missing, collapse = ", "))
#   }
#   if (anyDuplicated(colnames(mats[[nm]]))) {
#     stop(nm, ": duplicated column names found")
#   }
# }
# 
# # ---- Subset, order columns to match the sample sheet, and write -------------
# for (nm in names(mats)) {
#   m <- mats[[nm]][, keep_ids, drop = FALSE]
#   
#   out <- data.frame(gene_id = rownames(m), m,
#                     check.names = FALSE, stringsAsFactors = FALSE)
#   
#   out_file <- file.path(out_dir, paste0(nm, ".txt"))
#   write.table(out, file = out_file, sep = "\t",
#               quote = FALSE, row.names = FALSE)
#   
#   cat(sprintf("Wrote %s  (%d genes x %d samples)\n",
#               out_file, nrow(m), ncol(m)))
# }
#   
# # ---- Confirm the Stage E sample was excluded --------------------------------
# for (nm in names(mats)) {
#   chk <- read.delim(file.path(out_dir, paste0(nm, ".txt")),
#                     nrows = 1, check.names = FALSE)
#   stopifnot(!("SR2518012" %in% colnames(chk)))
# }