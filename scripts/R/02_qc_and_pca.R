# Filtering counts, running DESeq2 VST normalization, performing PCA,
# and generating QC plots (PCA, sample-to-sample correlation heatmap)

# Activating project-specific renv environment
renv::load("scripts/R")

# Activate helper scripts
source("scripts/R/helper_scripts.R")

# Required packages
library(DESeq2)
library(dplyr)
library(ggplot2)
library(limma)
library(matrixStats)
library(pheatmap)
library(yaml)

# Load processed data from 01_data_preprocessing and Config
processed_data <- readRDS("data/processed_data.rds")
counts <- processed_data$counts
metadata <- processed_data$metadata

cat("Counts dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata dimensions:", nrow(metadata), "samples x", ncol(metadata), "variables\n")

path_config <- yaml::read_yaml("config/paths.yaml")
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
        filename = paste0("pca_", col_name, ".png")
    )
})

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

plot_and_save_pca(
    expr_matrix = mat, 
    metadata = metadata,
    color_by = "cardiac_phenotype",
    shape_by = "day",
    palettes = color_palettes$cardiac_phenotype,
    title = "PCA Coloured by Cardiac Phenotype shape Day",
    pca_output_dir = pca_design_batch_effect_removed,
    filename = "pca_cardiac_phenotype_day.png"
)
