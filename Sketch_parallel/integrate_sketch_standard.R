#' @file sketch_rpca_integration.R
#' @title Seurat Sketch RPCA Integration Script
#' @description This script performs RPCA-based layer integration on a Seurat sketch object.
#' @details The script initializes parallel processing via `future`, loads a Seurat object,
#'          applies `IntegrateLayers` using `RPCAIntegration`, and saves the resulting object.
#'          It is designed for processing large single-cell datasets efficiently.
#' @author [Your Name/Organization]
#' @keywords Seurat, RPCA, integration, single-cell, sketch
#' @import Seurat
#' @import future

library(Seurat)
library(future)

# Use the 8 cores
plan("multicore", workers = 1)
# Increased to 80GB to be safe, though DietSeurat should make the object tiny.
options(future.globals.maxSize = 100 * 1024^3)

print("--- SUB-PROCESS: LOADING SKETCH ---")

#' @title Initial Seurat Sketch Object
#' @description A Seurat object loaded from "temp_sketch_only.rds", representing a
#'              subsetted dataset ready for integration.
#' @format Seurat object
#' @source temp_sketch_only.rds
sketch_subset <- readRDS("temp_sketch_only.rds")

print("--- SUB-PROCESS: STARTING PARALLEL RPCA ---")

#' @title RPCA Integrated Seurat Sketch Object
#' @description The `sketch_subset` object after performing RPCA-based layer integration.
#' @details The `IntegrateLayers` function is applied to align different data layers
#'          within the object. It uses `RPCAIntegration` method, leveraging the
#'          "pca" reduction and creating a new "integrated.rpca" reduction.
#'          Integration is performed across dimensions 1 to 30, with an anchor
#'          number of 20 (`k.anchor = 20`).
#' @seealso \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{RPCAIntegration}}
sketch_subset <- IntegrateLayers(
  object = sketch_subset,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.anchor = 20,
  verbose = TRUE
)

print("--- SUB-PROCESS: SAVING RESULTS ---")

#' @title Save Integrated Sketch Object
#' @description Saves the final RPCA-integrated Seurat sketch object to disk.
#' @details The object is saved as an RDS file named "temp_sketch_integrated.rds"
#'          using R's `saveRDS` function, specified with `version = 3` and
#'          `compress = FALSE` for compatibility and performance.
#' @output temp_sketch_integrated.rds
saveRDS(sketch_subset, "temp_sketch_integrated.rds", version = 3, compress = FALSE)
