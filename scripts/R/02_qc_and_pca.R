# Filtering counts, running DESeq2 VST normalization, performing PCA,
# and generating QC plots (PCA, sample-to-sample correlation heatmap)

# Activating project-specific renv environment
renv::load("scripts/R")

# Activate helper scripts
source("scripts/R/helper_scripts.R")

# Required packages
library(DESeq2)
library(ggplot2)
library(limma)
library(matrixStats)
library(openxlsx)
library(pheatmap)
library(readxl)
library(tidyverse)
library(yaml)

# Load processed data from 01_data_preprocessing and Config
processed_data <- readRDS("data/processed_data.rds")
counts <- processed_data$counts
metadata <- processed_data$metadata

cat("Counts dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata dimensions:", nrow(metadata), "samples x", ncol(metadata), "variables\n")

path_config <- yaml::read_yaml("config/paths.yaml")
qc_dir <- path_config$qc
pca_blinded_to_design <- path_config$pca$blind_to_design
pca_design_batch_effect_removed <- path_config$pca$design_and_batch_effect_removed

# Filtering out lowly-expressed genes
# Keep genes with at least 10 counts in at least 10 samples
keep <- rowSums(counts >= 10) >= 10
filtered_counts <- counts[keep, ]
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(filtered_counts), "\n")
# Genes before filtering: 38341
# Genes after filtering: 19413

# Creating DESeqDataSet and Variance Stabilizing Transformation Normalization
# 1. Raw VST PCA (blind = TRUE)
# To assess overall sample quality, identify outliers, and 
# check whether any unexpected technical issues or sample swaps are present
dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = metadata,
    design = ~1
)

vsd <- vst(dds, blind = TRUE)
norm_counts <- assay(vsd)

# PCA colours by the different variables
cols <- c("condition", "day", "sequencing_batch", "sex", "cardiac_phenotype")
all_plots <- lapply(cols, function(col_name) {
    clean_title_name <- tools::toTitleCase(gsub("_", " ", col_name))
    plot_and_save_pca(
        expr_matrix = filtered_counts, 
        metadata = metadata,
        color_by = col_name,
        shape_by = NULL,
        palettes = color_palettes[[col_name]],
        title = paste("PCA Coloured by", clean_title_name),
        pca_output_dir = pca_blinded_to_design,
        filename = paste0("pca_", col_name)
    )
})

# PCA colour by phenotype and shape by batch
plot_and_save_pca(
    expr_matrix = filtered_counts, 
    metadata = metadata,
    color_by = "cardiac_phenotype",
    shape_by = "sequencing_batch",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Batch",
    pca_output_dir = pca_blinded_to_design,
    filename = "pca_cardiac_phenotype_batch"
)

# PCA colour by phenotype and shape by day
plot_and_save_pca(
    expr_matrix = filtered_counts, 
    metadata = metadata,
    color_by = "cardiac_phenotype",
    shape_by = "day",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Day",
    pca_output_dir = pca_blinded_to_design,
    filename = "pca_cardiac_phenotype_day"
)

# 2. Batch-corrected VST PCA (blind = FALSE) with design = ~ Batch + Sex + Day + Condition
# To visualize the biological structure of the data after accounting for technical batch effects
# Since condition and cardiac_phenotype is correlated, 
# we can drop condition in design as cardiac_phenotype has control as well as the two types of patients
check_design_collinearity(
    formula = ~ sequencing_batch + sex + age + day + cardiac_phenotype,
    data = as.data.frame(metadata)
)

# Setting Control as baseline
metadata$cardiac_phenotype <- relevel(metadata$cardiac_phenotype, ref = "Control")

dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = metadata,
    design = ~ sequencing_batch + sex + age + day + cardiac_phenotype
)

vsd <- vst(dds, blind = FALSE)

mat <- removeBatchEffect(
    assay(vsd),
    batch = metadata$sequencing_batch,
    design = model.matrix(~ sex + age + day + cardiac_phenotype, metadata)
)

# PCA colour by phenotype and shape by batch
plot_and_save_pca(
    expr_matrix = mat, 
    metadata = metadata,
    color_by = "cardiac_phenotype",
    shape_by = "sequencing_batch",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Batch",
    pca_output_dir = pca_design_batch_effect_removed,
    filename = "pca_cardiac_phenotype_batch"
)

# PCA colour by phenotype and shape by day
plot_and_save_pca(
    expr_matrix = mat, 
    metadata = metadata,
    color_by = "cardiac_phenotype",
    shape_by = "day",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Day",
    pca_output_dir = pca_design_batch_effect_removed,
    filename = "pca_cardiac_phenotype_day"
)

# Sample-to-Sample Correlation Heatmap using batch-effect corrected counts
sample_cor <- cor(mat)

annotation_df <- data.frame(
    Condition = metadata$condition,
    Timepoint = metadata$day,
    Batch = metadata$sequencing_batch,
    Phenotype = metadata$cardiac_phenotype,
    row.names = metadata$sample_id
)

annotation_colors <- list(
    Phenotype = color_palettes$cardiac_phenotype,
    Timepoint = color_palettes$day,
    Batch = color_palettes$sequencing_batch
)

heatmap_file <- file.path(qc_dir, "sample_correlation_heatmap_batch_effect_corrected.png")
png(
    heatmap_file, 
    width = fig_width * fig_dpi, 
    height = (fig_height + 1) * fig_dpi,
    res = fig_dpi
)
pheatmap(
    sample_cor,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = colorRampPalette(c("#457b9d", "#f1faee", "#e63946"))(100),
    main = "Sample-to-Sample Pearson Correlation Heatmap"
)
dev.off()

# Removing samples with either low mapping or low properly-paired or low exonic or high rRNA
mapping_metrics <- read_excel(
    file.path(qc_dir, "mapping_metrics/mapping_metrics_for_all_batches.xlsx"),
    sheet = "QC_compiled"
)
metadata_samples_to_keep_col_add <- metadata |> 
  rownames_to_column(var = "row_names_temp") |> 
  left_join(
    mapping_metrics |> select(Sample, `Samples to exclude`), 
    join_by("sample_id" == "Sample")) |> 
  column_to_rownames(var = "row_names_temp")

samples_to_exclude <- mapping_metrics$Sample[mapping_metrics$`Samples to exclude` == "No"]
# 12 samples to exclude
filtered_counts_samples_removed <- filtered_counts[, !colnames(filtered_counts) %in% samples_to_exclude]
metadata_samples_removed <- metadata[!metadata$sample_id %in% samples_to_exclude, ]

# QC on the filtered samples
dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts_samples_removed,
    colData = metadata_samples_removed,
    design = ~ sequencing_batch + sex + age + day + cardiac_phenotype
)

vsd <- vst(dds, blind = FALSE)

mat_samples_removed <- removeBatchEffect(
    assay(vsd),
    batch = metadata_samples_removed$sequencing_batch,
    design = model.matrix(~ sex + age + day + cardiac_phenotype, metadata_samples_removed)
)

# PCA colour by phenotype and shape by batch
plot_and_save_pca(
    expr_matrix = mat_samples_removed, 
    metadata = metadata_samples_removed,
    color_by = "cardiac_phenotype",
    shape_by = "sequencing_batch",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Batch",
    pca_output_dir = pca_design_batch_effect_removed,
    filename = "pca_cardiac_phenotype_batch_samples_removed"
)

# PCA colour by phenotype and shape by day
plot_and_save_pca(
    expr_matrix = mat_samples_removed, 
    metadata = metadata_samples_removed,
    color_by = "cardiac_phenotype",
    shape_by = "day",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Day",
    pca_output_dir = pca_design_batch_effect_removed,
    filename = "pca_cardiac_phenotype_day_samples_removed"
)

# Sample-to-Sample Correlation Heatmap using batch-effect corrected counts
sample_cor <- cor(mat_samples_removed)

annotation_df <- data.frame(
    Condition = metadata_samples_removed$condition,
    Timepoint = metadata_samples_removed$day,
    Batch = metadata_samples_removed$sequencing_batch,
    Phenotype = metadata_samples_removed$cardiac_phenotype,
    row.names = metadata_samples_removed$sample_id
)

annotation_colors <- list(
    Phenotype = color_palettes$cardiac_phenotype,
    Timepoint = color_palettes$day,
    Batch = color_palettes$sequencing_batch
)

heatmap_file <- file.path(qc_dir, "sample_correlation_heatmap_batch_effect_corrected_samples_removed.png")
png(
    heatmap_file, 
    width = fig_width * fig_dpi, 
    height = (fig_height + 1) * fig_dpi,
    res = fig_dpi
)
pheatmap(
    sample_cor,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = colorRampPalette(c("#457b9d", "#f1faee", "#e63946"))(100),
    main = "Sample-to-Sample Pearson Correlation Heatmap"
)
dev.off()

# Saving the processed data along with batch-effect corrected counts
processed_data$metadata <- metadata_samples_to_keep_col_add
processed_data$filtered_counts <- filtered_counts
processed_data$filtered_counts_batch_effect_corrected <- mat
processed_data$filtered_counts_batch_effect_corrected_samples_removed <- mat_samples_removed
processed_data$metadata_samples_removed <- metadata_samples_removed

saveRDS(processed_data, "data/processed_data_batch_effect_corrected_qced_samples.rds")

long_names <- names(processed_data)
short_names <- substr(long_names, 1, 28)
final_names <- make.unique(short_names, sep = ".")
names(processed_data) <- substr(final_names, 1, 31)
write.xlsx(processed_data, file = "data/processed_data_batch_effect_corrected_qced_samples.xlsx", rowNames = TRUE)

write.csv(installed.packages(), "environment/R_packages.csv")
yaml::write_yaml(devtools::session_info(), "environment/sessionInfo.yaml")
