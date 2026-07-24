# Day-specific marker analysis of cardiomyocytes
# Correlating expression with GAA repeat lengths in patient samples

# Activating project-specific renv environment
renv::load("scripts/R")

# Activate helper scripts
source("scripts/R/helper_scripts.R")

library(DESeq2)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(rstatix)
library(tidyverse)
library(yaml)

# Load Processed Data and Config
processed_data <- readRDS("data/processed_data_batch_effect_corrected_qced_samples.rds")
counts <- processed_data$filtered_counts_batch_effect_corrected
counts_samples_removed <- processed_data$filtered_counts_batch_effect_corrected_samples_removed
metadata <- processed_data$metadata
metadata_samples_removed <- processed_data$metadata_samples_removed

# Cardiac maturation and identity markers
markers <- c(
  # Pluripotent state
  "OCT4", "SOX2", "NANOG", "POU5F1",
  # Primitive Mesoderm
  "TBXT",
  # Cardiac Progenitors
  "ISL1", "KDR", "MESP1",
  # Early Cardiomyocytes
  "NKX2-5", "TNNT2", "ACTN2", "MYH6",
  # Maturing Cardiomyocytes
  "MYL2", "MYL7", # Day 12-15+ and Day 20-40
  "GJA1", # GJA1 encoding Connexin 43 # Day 12-15+
  "GATA4", "MEF2C", # Day 12-15+
  # Contractile & Sarcomere
  "TNNT2", "MYH6", "MYH7", "ACTC1", "ACTN2", "MYL2", "MYL7",
  # Calcium Handling
  "RYR2", "PLN", "ATP2A2", "CACNA1C", "SLC8A1",
  # Transcription Factors
  "NKX2-5", "GATA4", "MEF2C", "TBX5",
  # Natriuretic Peptides (Stress/development)
  "NPPA", "NPPB"
)

markers_present <- markers[markers %in% rownames(counts_samples_removed)] |> unique()
# OCT4 absent
# Brachyury, a transcription factor 
# encoded by Gene TBXT is filtered out due to low counts
cat("Defined", length(unique(markers)), "markers. Present in dataset:", length(unique(markers_present)), "\n")
# Total 25 markers out of 27 defined are present in the dataset

marker_exp <- counts_samples_removed[markers_present, , drop = FALSE]
marker_df <- as.data.frame(marker_exp) |> 
  rownames_to_column("Gene") |> 
  pivot_longer(-Gene, names_to = "sample_id", values_to = "Expression") |> 
  left_join(metadata_samples_removed, by = "sample_id")

# Comparing expression level of marker genes between cardiac phenotypes within and between day 15 and day 52
# Using Wilcoxon test for comparison. Reasons:
# t-test is sensitive to outliers
# Sample size is less to follow the CLT
# Gene expression data rarely follows a bell-curve normal distribution
stat_df <- marker_df |> 
  mutate(group_comb = paste(day, cardiac_phenotype, "_"))

stats_within_days <- stat_df |> 
  group_by(Gene, day) |> 
  wilcox_test(Expression ~ cardiac_phenotype) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance() |> 
  mutate(comparison_type = paste0("Within Day Phenotype Comparison Between ", day))

stats_across_days <- stat_df %>%
  group_by(Gene, cardiac_phenotype) %>%
  wilcox_test(Expression ~ day) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  mutate(comparison_type = paste("Across Days Timepoint Comparison for ", cardiac_phenotype))

combined_stats <- bind_rows(stats_within_days, stats_across_days) %>%
  select(Gene, comparison_type, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)

comparison_file <- file.path(marker_analysis_dir, "marker_genes_significance_tests.xlsx")
write.xlsx(combined_stats, file = comparison_file)

# Boxplot of marker genes between cardiac phenotypes stratified by day
p_box <- ggplot(marker_df, aes(x = day, y = Expression, fill = cardiac_phenotype)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85, width = 0.6, position = position_dodge(0.7)) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = color_palettes$cardiac_phenotype) +
  labs(
    title = "Expression of Cardiomyocyte Maturation & Identity Markers",
    subtitle = "Batch-effect corrected counts across Day 15 and Day 52",
    x = "Timepoint (Day)",
    y = "Expression",
    fill = "Cardiac Phenotype"
  ) +
  theme_custom() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9))

boxplots_file <- file.path(marker_analysis_dir, "maturation_markers_boxplot.png")
ggsave(
  filename = boxplots_file,
  plot = p_box,
  width = fig_width + 5,
  height = fig_height + 4,
  dpi = fig_dpi
)

# Expression Heatmap for each sample
# Order samples by Day, then Condition, then Cardiac Phenotype to highlight maturation
ordered_metadata <- metadata |> 
  arrange(day, condition, cardiac_phenotype)

heatmap_data <- counts[markers_present, ordered_metadata$sample_id, drop = FALSE]

# Row-wise Z-score normalization
heatmap_data_z <- t(scale(t(heatmap_data)))

annotation_df <- data.frame(
  Condition = ordered_metadata$condition,
  Timepoint = ordered_metadata$day,
  Phenotype = ordered_metadata$cardiac_phenotype,
  Exclusion = ordered_metadata$`Samples to include`,
  row.names = ordered_metadata$sample_id
)

annotation_colors <- list(
  Condition = color_palettes$condition,
  Timepoint = color_palettes$day,
  Phenotype = color_palettes$cardiac_phenotype,
  Exclusion = c("Yes" = "#1b4332", "No" = "#d62828") 
)

heatmap_file <- file.path(marker_analysis_dir, "maturation_markers_heatmap.png")
png(heatmap_file, width = (fig_width + 10) * fig_dpi, height = (fig_height + 1) * fig_dpi, res = fig_dpi)
pheatmap(
  heatmap_data_z,
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE, # Maintain chronological ordering in ordered_metadata
  cluster_rows = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = colorRampPalette(c("#457b9d", "#f1faee", "#e63946"))(100),
  main = "Cardiac Maturation Marker Expression (Z-score)"
)
dev.off()

# Correlating patient expression with GAA repeat length (LA and UA)
patient_meta_day15 <- metadata_samples_removed |> dplyr::filter(condition == "Patient", day == "D15")
patient_counts_day15 <- counts_samples_removed[markers_present, patient_meta_day15$sample_id, drop = FALSE]

patient_meta_day52 <- metadata_samples_removed |> dplyr::filter(condition == "Patient", day == "D52")
patient_counts_day52 <- counts_samples_removed[markers_present, patient_meta_day52$sample_id, drop = FALSE]

cat("Analyzing correlations for", nrow(patient_meta_day15), "Day 15 patient samples and", nrow(patient_meta_day52), "Day 52 patient samples\n")
# 37 Day 15 patient samples
# 16 Day 52 patient samples
# 53 patient samples

cor_results <- list()
for (gene in markers_present) {
  gene_expr_day15 <- patient_counts_day15[gene, ]
  gene_expr_day52 <- patient_counts_day52[gene, ]

  # Pearson correlation
  # Day 15
  p_la_15 <- cor(gene_expr_day15, patient_meta_day15$LA, method = "pearson", use = "complete.obs")
  p_ua_15 <- cor(gene_expr_day15, patient_meta_day15$UA, method = "pearson", use = "complete.obs")

  p_val_la_15 <- cor.test(gene_expr_day15, patient_meta_day15$LA, method = "pearson")$p.value
  p_val_ua_15 <- cor.test(gene_expr_day15, patient_meta_day15$UA, method = "pearson")$p.value

  # Day 52
  p_la_52 <- cor(gene_expr_day52, patient_meta_day52$LA, method = "pearson", use = "complete.obs")
  p_ua_52 <- cor(gene_expr_day52, patient_meta_day52$UA, method = "pearson", use = "complete.obs")

  p_val_la_52 <- cor.test(gene_expr_day52, patient_meta_day52$LA, method = "pearson")$p.value
  p_val_ua_52 <- cor.test(gene_expr_day52, patient_meta_day52$UA, method = "pearson")$p.value

  cor_results[[gene]] <- data.frame(
    Gene = gene,
    Pearson_LA_D15 = p_la_15,
    Pearson_UA_D15 = p_ua_15,
    Pvalue_LA_D15 = p_val_la_15,
    Pvalue_UA_D15 = p_val_ua_15,
    Pearson_LA_D52 = p_la_52,
    Pearson_UA_D52 = p_ua_52,
    Pvalue_LA_D52 = p_val_la_52,
    Pvalue_UA_D52 = p_val_ua_52
  )
}

cor_df <- bind_rows(cor_results)
cor_csv_file <- file.path(marker_analysis_dir, "marker_GAA_correlations.csv")
write.csv(cor_df, cor_csv_file, row.names = FALSE)

# Plot correlation heatmap
cor_matrix <- cor_df |> 
  dplyr::select(Pearson_LA_D15, Pearson_UA_D15, Pearson_LA_D52, Pearson_UA_D52) |> 
  as.matrix()
rownames(cor_matrix) <- cor_df$Gene

cor_heatmap_file <- file.path(marker_analysis_dir, "marker_GAA_correlation_heatmap.png")
png(cor_heatmap_file, width = (fig_width - 1) * fig_dpi, height = fig_height * fig_dpi, res = fig_dpi)
pheatmap(
  cor_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  angle_col = 0,
  fontsize_col = 8,
  color = colorRampPalette(c("#1d3557", "#f1faee", "#e63946"))(100),
  main = "Correlation of Markers with GAA Repeat Lengths (Patients)"
)
dev.off()

write.csv(installed.packages(), "environment/R_packages.csv")
yaml::write_yaml(devtools::session_info(), "environment/sessionInfo.yaml")