# Load config files
library(yaml)

plotting_config <- yaml::read_yaml("config/plotting.yaml")
fig_width <- plotting_config$figure_width
fig_height <- plotting_config$figure_height
fig_dpi <- plotting_config$dpi
pca_ntop <- plotting_config$PCA$ntop

# Convert metadata colums
# Categorical variables = Factor
# Continuous variables = Numeric
prepare_metadata_types <- function(metadata, factor_cols = NULL, numeric_cols = NULL) {

  metadata <- as.data.frame(metadata)
  
  # 1. Convert specified columns to factors
  if (!is.null(factor_cols)) {
    for (col in factor_cols) {
      if (col %in% colnames(metadata)) {
        metadata[[col]] <- as.factor(metadata[[col]])
      } else {
        warning(paste("Column", col, "not found in metadata. Skipping factor conversion."))
      }
    }
  }
  
  # 2. Convert specified columns to numeric
  if (!is.null(numeric_cols)) {
    for (col in numeric_cols) {
      if (col %in% colnames(metadata)) {
        metadata[[col]] <- as.numeric(as.character(metadata[[col]]))
      } else {
        warning(paste("Column", col, "not found in metadata. Skipping numeric conversion."))
      }
    }
  }
  
  return(metadata)
}

# Check Design Matrix
# Check whether the model matrix is full rank or not 
# (i.e. there are no collinearities between covariates in the design formula)
check_design_collinearity <- function(formula, data) {
  design_matrix <- model.matrix(formula, data = data)
  # Check rank
  is_full_rank <- qr(design_matrix)$rank == ncol(design_matrix) # QR decomposition to check rank
  cat("Design matrix full rank: ", is_full_rank, "\n")
  if (!is_full_rank) {
    cat("\nDetected linear dependencies:\n")
    lm_obj <- lm(rep(1, nrow(design_matrix)) ~ design_matrix) # Fit linear model to identify dependencies
    print(alias(lm_obj)) # Display linear dependencies
  } else {
    cat("No linear dependencies detected.\n")
  }
  invisible(design_matrix)
}

# Define color palettes
color_palettes <- list(
  condition = c("Control" = "#2ec4b6", "Patient" = "#e71d36"),
  day = c("D15" = "#ff9f1c", "D52" = "#4361ee"),
  sequencing_batch = c("Batch1" = "#7209b7", "Batch2" = "#f72585", "Batch3" = "#4cc9f0"),
  sex = c("M" = "#0077b6", "F" = "#ffb703"),
  cardiac_phenotype = c("Control" = "#2ec4b6", "CMP" = "#d90429", "NCMP" = "#ffb703")
)

# Custom theme
theme_premium <- function() {
  theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = "gray90", linewidth = 0.5),
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5)
    )
}

# Plot: PCA
plot_and_save_pca <- function(expr_matrix, metadata, top_n = pca_ntop, color_by, shape_by, palettes, title, pca_output_dir, filename) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(matrixStats)
    })

    expr_matrix <- as.matrix(expr_matrix)
    metadata <- metadata[colnames(expr_matrix), , drop = FALSE]
    topVarGenes <- order(matrixStats::rowVars(expr_matrix), decreasing = TRUE)[seq_len(top_n)]
    
    pca_res <- prcomp(t(expr_matrix[topVarGenes, ]))
    percent_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
    
    pca_data <- data.frame(
        Sample = colnames(expr_matrix),
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2],
        Color = metadata[[color_by]]
    )
    
    if (!is.null(shape_by)) {
        pca_data$Shape <- metadata[[shape_by]]
    }
    
    pca_data <- na.omit(pca_data)
  
    counts <- table(pca_data$Color)
    legend_labels <- paste0(names(palettes), " (n = ", counts[names(palettes)], ")")
    names(legend_labels) <- names(palettes)
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
         scale_color_manual(values = palettes, labels = legend_labels) +
         labs(
             title = title,
             subtitle = paste0("Top ", top_n, " most variable genes (VST normalized)"),
             x = paste0("PC1: ", round(percent_var[1] * 100, 1), "% variance"),
             y = paste0("PC2: ", round(percent_var[2] * 100, 1), "% variance"),
             color = tools::toTitleCase(color_by)
         ) +
         theme_premium()
  
    if (!is.null(shape_by)) {
      shape_counts <- table(pca_data$Shape)
      unique_shapes <- names(shape_counts)

      shape_labels <- paste0(unique_shapes, " (n = ", shape_counts[unique_shapes], ")")
      names(shape_labels) <- unique_shapes

      p <- p + 
        geom_point(aes(color = Color, shape = Shape), size = 3, alpha = 0.8) +
        scale_shape_manual(values = 15:(15 + length(unique_shapes) - 1), labels = shape_labels) +
        labs(shape = tools::toTitleCase(shape_by))
    } else {
      p <- p + 
        geom_point(aes(color = Color), size = 3, alpha = 0.8)
    }
      
    ggsave(
      filename = file.path(pca_output_dir, filename),
      plot = p,
      width = fig_width,
      height = fig_height,
      dpi = fig_dpi
    )
    cat("Saved PCA plot to:", file.path(pca_output_dir, filename), "\n")
    
    return(p)
}