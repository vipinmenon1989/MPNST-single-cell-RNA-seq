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