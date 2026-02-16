library(Seurat)
library(future)

# Use the 8 cores
plan("multicore", workers = 1)
# Increased to 80GB to be safe, though DietSeurat should make the object tiny.
options(future.globals.maxSize = 100 * 1024^3) 

print("--- SUB-PROCESS: LOADING SKETCH ---")
sketch_subset <- readRDS("temp_sketch_only.rds")

print("--- SUB-PROCESS: STARTING PARALLEL RPCA ---")

# FIX: dims = 1:30 matches your PCA. 
# k.anchor = 20 is explicit (fixes the "multiple arguments" crash).
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
saveRDS(sketch_subset, "temp_sketch_integrated.rds", version = 3, compress = FALSE)