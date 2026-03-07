#' @title Single-Cell RNA-Seq QC, Filtering, and Pre-Integration Processing
#' @description This script performs comprehensive quality control, filtering, normalization, variable feature finding, scaling, and PCA on a Seurat object.
#' @details It is designed for single-cell RNA-seq data analysis workflows, producing both pre- and post-filtering statistics and visualizations. The script includes setup for HPC parallel processing and memory management for large datasets.
#'
#' The main steps include:
#' \itemize{
#'   \item Setting up `future` for parallel processing on HPC.
#'   \item Loading an aggregated Seurat object.
#'   \item Generating pre-filtering QC plots and summary statistics, including 'lab' information.
#'   \item Filtering cells based on `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.
#'   \item Generating post-filtering QC plots and summary statistics.
#'   \item Normalizing data, finding variable features, scaling data (HVG only), and running PCA.
#'   \item Saving the processed Seurat object as an RDS checkpoint before integration steps.
#' }
#'
#' @section Output Files:
#' \itemize{
#'   \item `01_QC_Pre_Scatter_Uniform.png`: Scatter plots showing cell counts vs mitochondrial percentage and cell counts vs gene features before filtering.
#'   \item `01_QC_Pre_Filtering_Stats.csv`: CSV file summarizing QC metrics (number of cells, median features, counts, MT%) per batch and lab before filtering.
#'   \item `02_QC_Post_Scatter_Uniform.png`: Scatter plots showing cell counts vs mitochondrial percentage and cell counts vs gene features after filtering.
#'   \item `02_QC_Post_Filtering_Stats.csv`: CSV file summarizing QC metrics per batch and lab after filtering.
#'   \item `seurat_obj_pre_integration.rds`: The final Seurat object saved as an RDS file, ready for integration.
#' }
#' @seealso
#' `Seurat` for single-cell object manipulation, `ggplot2` for visualization, `patchwork` for plot arrangement, `future` for parallel processing, `dplyr` for data wrangling.
#' @keywords single-cell RNA-seq QC filtering Seurat HPC
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import future
#' @import dplyr

library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

#' @section 1. HPC & PARALLEL SETUP:
#' Configures the `future` package for parallel processing on High-Performance Computing (HPC) environments using a multicore strategy.
#' Sets the maximum memory size for global variables to prevent allocation issues with large Seurat objects.
#'
#' We use 'multicore' for Linux HPCs.
plan("multicore", workers = 4)
#' Set max memory to 200GB (or 300GB if needed)
options(future.globals.maxSize = 300 * 1024^3)
print("--- HPC MODE: 4 CORES ENABLED ---")

#' @section 2. LOAD DATA:
#' Defines the file path for the input Seurat object and the column used to identify batches.
#' Loads the Seurat object from the specified RDS file.
file_path <- "aggregated_cellbender_filtered.rds"
batch_column <- "orig.ident"

seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

#' @section 3. PRE-FILTERING (UNIFORM BLACK PLOTS):
#' Calculates the percentage of mitochondrial genes for each cell and generates pre-filtering quality control plots.
#' Plots include cell counts versus mitochondrial percentage and cell counts versus the number of detected genes.
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

print("--- GENERATING PRE-QC PLOTS ---")

#' Plotting using ggplot2 for uniform black color and high contrast.
p1 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=percent.mt)) +
      geom_point(color="black", size=0.5) +
      theme_classic() +
      ggtitle("Pre-Filter: Counts vs Mito%") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

p2 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point(color="black", size=0.5) +
      theme_classic() +
      ggtitle("Pre-Filter: Counts vs Genes") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

ggsave("01_QC_Pre_Scatter_Uniform.png", plot = p1 + p2, width = 16, height = 8, dpi = 300)

#' @subsection CSV SUMMARY (PRE) - MODIFIED TO INCLUDE 'lab':
#' Generates a CSV summary of pre-filtering QC metrics, grouped by batch and 'lab' information.
print("--- GENERATING PRE-FILTER CSV (WITH LAB INFO) ---")

#' We group by the batch column AND 'lab'.
pre_stats <- seurat_obj@meta.data %>%
  group_by(.data[[batch_column]], lab) %>%
  summarise(
    Num_Cells = n(),
    Median_nFeature = median(nFeature_RNA),
    Median_nCount = median(nCount_RNA),
    Median_MT = median(percent.mt),
    .groups = 'drop'
  )

#' Create TOTAL row with a placeholder for 'lab'.
total_pre <- data.frame(
  batch_column = "TOTAL",
  lab = "ALL",  # Placeholder for the total row
  Num_Cells = ncol(seurat_obj),
  Median_nFeature = median(seurat_obj$nFeature_RNA),
  Median_nCount = median(seurat_obj$nCount_RNA),
  Median_MT = median(seurat_obj$percent.mt)
)
#' Rename the first column to match the batch column name (orig.ident).
colnames(total_pre)[1] <- batch_column

#' Combine and Save pre-filtering statistics.
write.csv(rbind(pre_stats, total_pre), "01_QC_Pre_Filtering_Stats.csv", row.names = FALSE)
print("--- PRE-FILTER METRICS SAVED ---")


#' @section 4. FILTERING:
#' Filters the Seurat object based on defined thresholds for `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.
#' Recovers memory after filtering.
print("--- FILTERING DATA ---")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 500 &
                              nFeature_RNA < 9000 &
                              percent.mt < 20)
#' Free up memory
gc()
print(paste("Cells remaining:", ncol(seurat_obj)))

#' @subsection POST-FILTERING PLOTS:
#' Generates post-filtering quality control plots, similar to the pre-filtering plots, to visualize the effect of filtering.
print("--- GENERATING POST-QC PLOTS ---")

p3 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=percent.mt)) +
      geom_point(color="black", size=0.5) +
      theme_classic() +
      ggtitle("Post-Filter: Counts vs Mito%") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

p4 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point(color="black", size=0.5) +
      theme_classic() +
      ggtitle("Post-Filter: Counts vs Genes") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

ggsave("02_QC_Post_Scatter_Uniform.png", plot = p3 + p4, width = 16, height = 8, dpi = 300)

#' @subsection CSV SUMMARY (POST) - MODIFIED TO INCLUDE 'lab':
#' Generates a CSV summary of post-filtering QC metrics, grouped by batch and 'lab' information.
print("--- GENERATING POST-FILTER CSV (WITH LAB INFO) ---")

post_stats <- seurat_obj@meta.data %>%
  group_by(.data[[batch_column]], lab) %>%
  summarise(
    Num_Cells = n(),
    Median_nFeature = median(nFeature_RNA),
    Median_nCount = median(nCount_RNA),
    Median_MT = median(percent.mt),
    .groups = 'drop'
  )

total_post <- data.frame(
  batch_column = "TOTAL",
  lab = "ALL",
  Num_Cells = ncol(seurat_obj),
  Median_nFeature = median(seurat_obj$nFeature_RNA),
  Median_nCount = median(seurat_obj$nCount_RNA),
  Median_MT = median(seurat_obj$percent.mt)
)
colnames(total_post)[1] <- batch_column

write.csv(rbind(post_stats, total_post), "02_QC_Post_Filtering_Stats.csv", row.names = FALSE)
print("--- POST-FILTER METRICS SAVED ---")


#' @section 5. NORMALIZE & FIND FEATURES (NO SPLIT/JOIN):
#' Normalizes the gene expression data using a global scaling normalization method and identifies highly variable features.
#' This step prepares the data for dimensionality reduction and clustering.
print("--- STARTING NORMALIZATION ---")

#' KEY FIX: We do NOT JoinLayers (causes overflow).
#' We do NOT split (already split).
#' We just run NormalizeData, which automatically handles the existing split layers.
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)


#' @section 6. SCALE & PCA:
#' Scales the highly variable features to ensure equal variance across genes and performs Principal Component Analysis (PCA) for dimensionality reduction.
print("--- SCALING DATA (HVG ONLY) ---")
#' Scaling works on the split layers automatically in V5.
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

print("--- RUNNING PCA ---")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30, verbose = FALSE)

print("--- SAVING CHECKPOINT RDS (NO COMPRESSION) ---")
#' Keep compress = FALSE to prevent crash.
saveRDS(seurat_obj, "seurat_obj_pre_integration.rds", compress = FALSE)

print("--- SCRIPT COMPLETED ---")
