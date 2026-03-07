#' @file run_seurat_marker_analysis.R
#' @title Seurat Marker Analysis and Heatmap Generation Script for Large Datasets
#' @description This script performs marker gene identification and heatmap generation for a Seurat object across multiple clustering resolutions.
#' It implements a "subset-first" strategy to manage memory with large datasets (e.g., 1 million cells).
#' @details The script iterates through a predefined range of clustering resolutions. For each resolution, it subsets a maximum of 500 cells per cluster to create a temporary, memory-efficient object. Marker genes are then identified on this subsetted object. Top 5 positive markers per cluster are extracted and used to generate a heatmap, further downsampling cells to 50 per cluster for visualization. Error handling is included for each resolution.
#' @section Memory Management Strategy:
#' To handle large Seurat objects efficiently and prevent out-of-memory errors, the script employs a critical "subset-first" strategy:
#' \enumerate{
#'   \item Global memory limit: `options(future.globals.maxSize)` is set to a high value.
#'   \item Marker finding subset: For `FindAllMarkers`, cells are subsetted to a maximum of 500 per cluster. This is statistically sufficient and dramatically reduces memory usage.
#'   \item Heatmap visualization subset: For `DoHeatmap`, cells are further subsetted to a maximum of 50 per cluster for a cleaner visual representation.
#' }
#' This approach ensures that memory-intensive operations like `JoinLayers` and `FindAllMarkers` are performed on manageable object sizes.
#' @author Your Name/Organization (Placeholder)
#' @date 2023-10-27 (Placeholder)

# 1. SET MEMORY LIMITS
#' @description Sets the maximum global size for future operations.
#' @details This is crucial for handling large Seurat objects and preventing out-of-memory errors during parallel computations or large data manipulations.
options(future.globals.maxSize = 500 * 1024^3)

#' @description Loads necessary R packages for Seurat analysis, plotting, and data manipulation.
#' @import Seurat Provides single-cell RNA sequencing analysis tools.
#' @import ggplot2 Provides a system for declaratively creating graphics.
#' @import dplyr Provides a grammar of data manipulation.
#' @import data.table Provides an enhanced version of R's data.frame.
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

# 2. LOAD DATA
print("--- LOADING OBJECT ---")
#' @description Loads the pre-clustered Seurat object from a specified RDS file.
#' @param "seurat_obj_clustered_1M.rds" The file path to the input Seurat object.
#' @return A Seurat object named `seurat_obj` loaded into the global environment.
seurat_obj <- readRDS("seurat_obj_clustered_1M.rds")

# CRITICAL: Do NOT run JoinLayers(seurat_obj) here. It will crash.
#' @description Sets the default assay for the loaded Seurat object to "RNA".
#' @details This ensures that downstream functions operating on assay data (e.g., `FindAllMarkers`)
#' correctly target the RNA counts and data slots.
DefaultAssay(seurat_obj) <- "RNA"

#' @description Creates output directories for marker gene tables and heatmaps.
#' @details These directories are created if they do not already exist. `showWarnings = FALSE`
#' suppresses warnings if the directories already exist.
dir.create("marker_tables_final", showWarnings = FALSE)
dir.create("heatmap_outputs_final", showWarnings = FALSE)

# 3. DEFINE RESOLUTIONS
#' @description Defines the range of clustering resolutions to be processed.
#' @details The script will iterate through these resolutions, which correspond to metadata columns
#' in the format `integrated_snn_res.X`.
res_range <- seq(0.1, 1, by = 0.1)

print("--- STARTING SUBSET-FIRST PIPELINE ---")

#' @description Main loop to iterate through each defined clustering resolution.
#' @param res An individual resolution value from `res_range`.
#' @section Processing Steps within Loop:
#' \enumerate{
#'   \item Construct the metadata column name for the current resolution.
#'   \item Validate the existence of the resolution-specific metadata column.
#'   \item Set the active identity of the Seurat object to the current resolution's clusters.
#'   \item Perform smart subsetting: Extract a maximum of 500 cells per cluster to create a temporary, memory-efficient object (`sub_obj`).
#'   \item Join layers on the subsetted object (`sub_obj`).
#'   \item Find all positive marker genes for each cluster in `sub_obj`.
#'   \item Save the identified marker genes to a CSV file.
#'   \item Prepare data for heatmap: Extract the top 5 positive marker genes per cluster.
#'   \item Generate and save a heatmap: Further subset cells to 50 per cluster for visualization and save as PDF.
#'   \item Clean up temporary objects (`sub_obj`, `heatmap_obj`, `markers`) and perform garbage collection (`gc()`) to free RAM.
#' }
for (res in res_range) {
  
  #' @description Constructs the metadata column name corresponding to the current resolution.
  res_col <- paste0("integrated_snn_res.", res)
  print(paste(">>> PROCESSING RESOLUTION:", res))
  
  # A. VALIDATE METADATA
  #' @description Validates if the resolution-specific metadata column exists in the Seurat object.
  #' @details If the column is not found (e.g., if a resolution was not run during clustering),
  #' a message is printed, and the current resolution is skipped, moving to the next in the range.
  if (!res_col %in% colnames(seurat_obj@meta.data)) {
    print(paste("   Skipping:", res_col, "- not found."))
    next
  }
  
  # Set Identity to current resolution
  #' @description Sets the active identity of the Seurat object to the clusters defined by the current resolution.
  #' @param seurat_obj The Seurat object.
  #' @param value The metadata column containing cluster assignments for the current resolution.
  Idents(seurat_obj) <- seurat_obj@meta.data[[res_col]]
  
  # B. SMART SUBSET (THE FIX)
  # We extract max 500 cells per cluster.
  # This reduces the matrix size from 1,000,000 to ~15,000 (if 30 clusters).
  # This is statistically sufficient for marker detection and safe for memory.
  print("   Subsetting 500 cells/cluster to bypass integer limit...")
  
  #' @description Selects a maximum of 500 cells per cluster to create a lightweight temporary object.
  #' @details This step is critical for handling large datasets by significantly reducing the memory footprint
  #' for subsequent operations like `JoinLayers` and `FindAllMarkers`, preventing crashes while
  #' maintaining statistical power for marker detection.
  #' @param levels(Idents(seurat_obj)) The unique cluster identities in the current resolution.
  #' @return A character vector `cells_to_keep` containing the barcodes of the selected cells.
  cells_to_keep <- unlist(lapply(levels(Idents(seurat_obj)), function(x) {
    cells_in_cluster <- WhichCells(seurat_obj, idents = x)
    if (length(cells_in_cluster) > 500) {
      sample(cells_in_cluster, 500)
    } else {
      cells_in_cluster
    }
  }))
  
  # Create a lightweight temporary object
  #' @description Creates a temporary, subsetted Seurat object from the original object using selected cells.
  #' @param seurat_obj The original Seurat object.
  #' @param cells A character vector of cell barcodes (`cells_to_keep`) to include in the subset.
  #' @return A new Seurat object `sub_obj` containing only the specified cells.
  sub_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  # C. JOIN LAYERS ON SUBSET
  # This works because sub_obj is small.
  # It unifies the split layers into one 'counts' matrix for FindAllMarkers.
  #' @description Joins the split assay layers of the subsetted Seurat object into a single 'counts' matrix.
  #' @details This operation is necessary for `FindAllMarkers` to function correctly when assays contain multiple layers.
  #' It is performed on the small `sub_obj` to avoid memory issues that would occur on the full `seurat_obj`.
  #' @param sub_obj The subsetted Seurat object with potentially multiple assay layers.
  #' @return The `sub_obj` with assay layers unified.
  sub_obj <- JoinLayers(sub_obj)
  
  # D. FIND MARKERS
  print("   Calculating markers on subset...")
  
  #' @description Attempts to find all positive marker genes for each cluster in the subsetted object.
  #' @details This block is wrapped in a `tryCatch` to gracefully handle potential errors during marker finding
  #' or subsequent heatmap generation for a specific resolution, preventing the entire script from crashing.
  tryCatch({
    # No need for max.cells.per.ident here because sub_obj is already downsampled
    #' @description Identifies differentially expressed genes that are positive markers for each cluster.
    #' @param sub_obj The subsetted and layer-joined Seurat object.
    #' @param only.pos If TRUE, returns only positive markers (upregulated in the current cluster).
    #' @param min.pct Only test genes that are detected in a minimum fraction of cells in either of the two groups.
    #' @param logfc.threshold Only test genes that show a minimum log2 fold change between the two groups.
    #' @param verbose Suppresses verbose output from the function.
    #' @return A data frame `markers` containing the identified marker genes, their statistics, and cluster assignments.
    markers <- FindAllMarkers(
      sub_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = FALSE
    )
    
    #' @description Checks if any marker genes were found by `FindAllMarkers`.
    #' @details If no markers are found (the `markers` data frame is empty), a warning is printed,
    #' and the current resolution processing is skipped.
    if (nrow(markers) == 0) {
      print("   WARNING: No markers found. Skipping.")
      next
    }
    
    # Save Markers
    #' @description Saves the identified marker genes data frame to a CSV file.
    #' @param markers The data frame of marker genes.
    #' @param file The output file path, dynamically named to include the current resolution.
    #' @param row.names A logical value indicating whether row names should be written to the file. Set to `FALSE`.
    write.csv(
      markers, 
      file = paste0("marker_tables_final/markers_res_", res, ".csv"),
      row.names = FALSE
    )
    
    # E. PREPARE HEATMAP DATA
    #' @description Ensures the 'cluster' column exists in the markers data frame before proceeding.
    if (!"cluster" %in% colnames(markers)) next
    
    #' @description Extracts the top 5 positive marker genes (by `avg_log2FC`) for each cluster.
    #' @details This step prepares a curated list of genes to be visualized in the heatmap,
    #' focusing on the most representative markers. Unique genes are selected across all clusters.
    #' @param markers The data frame of marker genes.
    #' @return A character vector `top5_genes` containing unique gene names.
    top5_genes <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 5, order_by = avg_log2FC) %>%
      pull(gene) %>%
      unique()
      
    #' @description Checks if any top genes were identified for heatmap plotting.
    #' @details If no top genes are found (e.g., if there are no markers or an error in processing),
    #' the heatmap generation for the current resolution is skipped.
    if (length(top5_genes) == 0) next
    
    # F. GENERATE HEATMAP
    # We downsample again to 50 cells/cluster for the visual plot (cleaner)
    print("   Generating Heatmap...")
    #' @description Subsets cells again to a maximum of 50 per cluster specifically for heatmap visualization.
    #' @details This further reduction in cell count makes the heatmap cleaner, more interpretable,
    #' and reduces computation time for plotting without sacrificing the representation of gene expression patterns.
    #' @param levels(Idents(sub_obj)) The unique cluster identities in the subsetted object.
    #' @return A character vector `heatmap_cells` containing cell barcodes for the heatmap object.
    heatmap_cells <- unlist(lapply(levels(Idents(sub_obj)), function(x) {
      cells <- WhichCells(sub_obj, idents = x)
      if (length(cells) > 50) sample(cells, 50) else cells
    }))
    
    #' @description Creates a new Seurat object (`heatmap_obj`) containing the further subsetted cells.
    #' Scales the expression data for the `top5_genes` within this object for heatmap visualization.
    #' @param sub_obj The previously subsetted Seurat object (500 cells/cluster).
    #' @param cells The cell barcodes to include (`heatmap_cells`).
    #' @param features The specific genes to scale (`top5_genes`).
    #' @param verbose Suppresses verbose output from `ScaleData`.
    #' @return A Seurat object `heatmap_obj` with scaled data for the specified genes.
    heatmap_obj <- subset(sub_obj, cells = heatmap_cells)
    heatmap_obj <- ScaleData(heatmap_obj, features = top5_genes, verbose = FALSE)
    
    #' @description Generates a heatmap of the top marker genes for the current resolution.
    #' @param heatmap_obj The Seurat object prepared for heatmap, with scaled data.
    #' @param features The genes to plot on the heatmap (`top5_genes`).
    #' @param size Font size for cell labels.
    #' @param angle Angle for feature labels (gene names).
    #' @return A `ggplot2` object `p` representing the heatmap with a title and adjusted y-axis text size.
    p <- DoHeatmap(
      heatmap_obj, 
      features = top5_genes, 
      size = 3, 
      angle = 90
    ) + 
    ggtitle(paste("Resolution:", res)) +
    theme(axis.text.y = element_text(size = 5))

    #' @description Calculates an appropriate height for the output heatmap PDF.
    #' @details The height is determined dynamically, being at least 8 inches, and increasing
    #' based on the number of `top5_genes` to ensure all gene labels are visible.
    #' @param top5_genes The list of genes being plotted on the heatmap.
    #' @return A numeric value `calc_height` in inches for the PDF height.
    calc_height <- max(8, length(top5_genes) * 0.15)
    
    #' @description Saves the generated heatmap as a PDF file.
    #' @param filename The output file path, including the resolution in the filename.
    #' @param plot The `ggplot2` object `p` containing the heatmap.
    #' @param width Width of the PDF in inches.
    #' @param height Height of the PDF in inches (dynamically calculated `calc_height`).
    #' @param limitsize A logical value indicating whether to enforce the specified `width` and `height` limits. Set to `FALSE`.
    ggsave(
      filename = paste0("heatmap_outputs_final/Heatmap_res_", res, ".pdf"),
      plot = p, 
      width = 15, 
      height = calc_height,
      limitsize = FALSE
    )
    
    print(paste("   SUCCESS: Saved Res", res))
    
  }, error = function(e) {
    #' @description Error handler function for the `tryCatch` block.
    #' @details If an error occurs during marker finding or heatmap generation for a specific resolution,
    #' this function prints an error message containing the resolution and the error details.
    #' This allows the script to continue processing other resolutions.
    #' @param e The error object containing details about the error.
    print(paste("   ERROR at Res", res, ":", e$message))
  })
  
  # Cleanup to free RAM
  #' @description Removes temporary Seurat objects and marker data frames from memory.
  #' @details This step is crucial for long-running loops, especially when processing large datasets,
  #' to prevent memory accumulation and potential out-of-memory errors by explicitly freeing up RAM
    #' after each resolution's processing is complete. `gc()` forces garbage collection.
  rm(sub_obj, heatmap_obj, markers)
  gc()
}

print("--- PIPELINE FINISHED ---")
