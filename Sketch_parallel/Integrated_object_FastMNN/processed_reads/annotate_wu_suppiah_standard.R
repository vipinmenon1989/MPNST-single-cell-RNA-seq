#' @title Generate UMAP comparison plots for Seurat object resolutions and annotations
#'
#' @description This script loads a pre-processed Seurat object, filters it based on a specific annotation,
#' and then generates a PDF containing side-by-side UMAP plots. The left UMAP displays
#' clustering results for various resolutions, while the right UMAP displays a specified
#' cell type annotation.
#'
#' @details The script iterates through a predefined set of clustering resolutions.
#' For each resolution, it generates two UMAP plots: one showing the clustering
#' based on that resolution for all cells, and another showing a specific annotation
#' for a filtered subset of cells (excluding 'unknown' or NA annotations).
#' The plots are combined using `patchwork` and saved into a single PDF file.
#'
#' @return A PDF file named "UMAP_Comparison_Clean_<annot>.pdf" is created in the specified working directory.
#'
#' @author V. Menon (adapted)
#' @date 2023-10-27
#' @keywords UMAP Seurat Clustering Annotation Plotting Visualization
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @seealso \code{\link[Seurat]{DimPlot}}, \code{\link[patchwork]{wrap_plots}}
#'
library(Seurat)
library(ggplot2)
library(patchwork)

#' Set working directory to locate input Seurat object and save the output PDF.
setwd("/local/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/Integrated_object_FastMNN/processed_reads/")

#' Load the pre-processed Seurat object containing integrated data, UMAP reduction, and metadata.
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

#' Configuration parameters for resolutions and annotation column.
#' Define a vector of clustering resolutions to iterate through for UMAP plotting.
resolutions <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")
#' Specify the name of the metadata column containing cell type annotations.
annot <- "wu_annot"

#' Filter the Seurat object for a clean annotation UMAP.
#' Identify cells where the specified annotation is NOT 'unknown', 'Unknown', or NA.
cells_to_keep <- which(!seurat_obj[[annot, drop = TRUE]] %in% c("unknown", "Unknown", NA))
#' Create a new Seurat object containing only the filtered cells for annotation UMAP.
seurat_annot_filtered <- seurat_obj[, cells_to_keep]

#' Initialize PDF device to save the UMAP comparison plots.
pdf(paste0("UMAP_Comparison_Clean_", annot, ".pdf"), width = 18, height = 8)

#' Loop through each defined resolution to generate and combine UMAP plots.
for (res in resolutions) {
  #' Construct the full column name for the current resolution's clustering.
  res_col <- paste0("integrated_snn_res.", res)
  
  #' Check if the resolution column exists in the Seurat object's metadata.
  if (res_col %in% colnames(seurat_obj@meta.data)) {
    
    #' Generate the left UMAP plot showing clusters for the current resolution.
    #' Displays all cells with clustering labels, with legend suppressed.
    p1 <- DimPlot(seurat_obj, reduction = "umap.mnn", group.by = res_col, 
                  label = TRUE, raster = TRUE) + 
          ggtitle(paste("Resolution:", res)) + 
          theme(legend.position = "none")

    #' Generate the right UMAP plot showing specified annotations for filtered cells.
    #' Displays filtered cells with annotation labels.
    p2 <- DimPlot(seurat_annot_filtered, reduction = "umap.mnn", group.by = annot, 
                  label = TRUE, raster = TRUE) + 
          ggtitle(paste("Annotation:", annot))

    #' Print the combined UMAP plots side-by-side to the PDF.
    print(p1 + p2)
  }
}
#' Close the PDF device, saving all generated plots to the file.
dev.off()
