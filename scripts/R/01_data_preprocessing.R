# Preprocessing the feature count matrices and metadata
# merging them, and saving the merged dataset as an RDS object

# Activate project-specific renv environment
renv::load("scripts/R")

# Load required packages
library(AnnotationDbi)
library(dplyr)
library(org.Hs.eg.db)
library(readxl)
library(tidyr)

# Load Path Configuration
paths_config <- yaml::read_yaml("config/paths.yaml")
metadata_path <- file.path(paths_config$metadata, "metadata_clean.xlsx")
counts_dir <- paths_config$count_data

# Load Metadata
metadata <- read_excel(metadata_path, sheet = "Metadata_Clean")

# Load Feature Count Matrices
# Day15_Batch1 = 44
# Day15_Batch2 = 31
# Day52_Batch3 = 49
# Total samples (Patients + Controls) = 124
count_files <- list.files(counts_dir, full.names = TRUE, pattern = "\\.txt$")
print(basename(count_files))

counts_list <- list()
for (f in count_files) {
  cat("Reading file:", basename(f), "\n")
  df <- read.delim(f, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  counts_list[[basename(f)]] <- df
}

# Ensembl ID to Gene Symbol Mapping
# Retain only Gene Symbols which includes Genes as well as lncRNA and miRNA Genes in GTF
# Total Ensembl IDs = 78691
mapped_counts_list <- list()

for (name in names(counts_list)) {
  df <- counts_list[[name]]

  # Strip version from Ensembl ID (e.g. ENSG00000310526.1 -> ENSG00000310526)
  ens_ids <- gsub("\\..*$", "", df$Geneid)

  # Mapping using org.Hs.eg.db
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ens_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  # Adding symbols back to dataframe
  df$Symbol <- symbols

  # Remove rows with NA Symbols
  mapped_df <- df |> filter(!is.na(Symbol) & Symbol != "")
  cat("Removed", nrow(df) - nrow(mapped_df), "unmapped Ensembl IDs. Remaining rows:", nrow(mapped_df), "\n")

  # Removing the original Geneid column and aggregating counts by Symbol (sum duplicates)
  agg_df <- mapped_df |> 
    select(-Geneid) |> 
    group_by(Symbol) |> 
    summarise(across(everything(), sum)) |> 
    ungroup()

  cat("Aggregated to unique Symbols. Rows:", nrow(agg_df), "\n")
  mapped_counts_list[[name]] <- agg_df
}
# Number of rows with Gene Symbols = 38512
# After aggregation (due to more than 1 row with same Gene Symbol) = 38341

# Merge all Count Matrices
# Key = Symbol
merged_counts <- mapped_counts_list[[1]]
for (i in 2:length(mapped_counts_list)) {
  merged_counts <- full_join(merged_counts, mapped_counts_list[[i]], by = "Symbol")
}

# Align with metadata
# Check if all metadata sample IDs are in the counts
missing_samples <- metadata$sample_id[!metadata$sample_id %in% colnames(merged_counts)]
if (length(missing_samples) > 0) {
  stop("Error: The following metadata samples are missing from counts: ", paste(missing_samples, collapse = ", "))
}

# Order and keep sample columns present in metadata
aligned_counts_df <- merged_counts |> 
  select(Symbol, all_of(metadata$sample_id))
all(colnames(count_matrix) == metadata$sample_id)

# Converting to matrix and set rownames
count_matrix <- as.matrix(aligned_counts_df[, -1])
rownames(count_matrix) <- aligned_counts_df$Symbol
# Dimension = 38341 x 101

# Saving the processed data
processed_data <- list(
  counts = count_matrix,
  metadata = metadata
)
saveRDS(processed_data, "data/processed_data.rds")

write.csv(installed.packages(), "environment/R_packages.csv")
yaml::write_yaml(devtools::session_info(), "environment/sessionInfo.yaml")
