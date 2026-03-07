#' @title Subsetting a Seurat Object by Lab for Memory Efficiency
#' @description This script provides a memory-efficient method to subset a large Seurat object
#'   into smaller objects based on unique 'lab' identifiers found in its metadata.
#'   Each lab-specific subset is saved as an individual RDS file.
#' @details The process involves:
#'   1. Temporarily loading the full Seurat object to extract all unique 'lab' identifiers.
#'   2. Creating a dedicated directory ('lab_subsets') for the output files.
#'   3. Iterating through each identified lab:
#'      a. Identifying cells belonging to the current lab.
#'      b. Creating a subset of the main object for these cells.
#'      c. Saving the subsetted object to an RDS file in the 'lab_subsets' directory.
#'      d. Immediately removing the subset object from memory and running garbage collection
#'         to free up resources before processing the next lab.
#'   Finally, the main Seurat object is removed from memory.
#' @return This script does not return an R object. Its primary outcome is the creation
#'   of multiple '.rds' files in the 'lab_subsets' directory, each containing a Seurat object
#'   corresponding to a specific lab.
#' @seealso \code{\link[Seurat]{readRDS}}, \code{\link[Seurat]{subset}}
#' @keywords Seurat, subsetting, memory, large datasets, lab, processing

library(Seurat)

# 1. READ METADATA ONLY (To get Lab Names and Cell IDs)
print("--- IDENTIFYING LABS ---")
# We load the object briefly just to get the metadata
temp_obj <- readRDS("seurat_obj_mnn_integrated_1M.rds")
all_labs <- unique(as.character(temp_obj$lab))
all_labs <- all_labs[!is.na(all_labs)]

# 2. THE MEMORY-EFFICIENT LOOP
if (!dir.exists("lab_subsets")) dir.create("lab_subsets")

for (current_lab in all_labs) {
    print(paste0("--- EXTRACTING LAB: ", current_lab, " ---"))
    
    # Identify barcodes for this lab while the object is in memory
    cells_to_keep <- colnames(temp_obj)[temp_obj$lab == current_lab]
    
    # Subset
    lab_subset <- subset(temp_obj, cells = cells_to_keep)
    
    # Save
    saveRDS(lab_subset, paste0("lab_subsets/lab_", current_lab, ".rds"))
    
    # Clean up the subset immediately
    rm(lab_subset)
    gc()
}

# Now clean up the master object
rm(temp_obj)
gc()
print("--- ALL LABS SAVED ---")
