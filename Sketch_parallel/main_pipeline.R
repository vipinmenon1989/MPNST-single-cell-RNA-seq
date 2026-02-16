library(Seurat)
library(ggplot2)
library(future)
library(dplyr)

# 1. SETUP
plan("sequential")
options(future.globals.maxSize = 450 * 1024^3) 

# 2. LOAD & SKETCH
print("--- LOADING DATA ---")
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

DefaultAssay(seurat_obj) <- "RNA"
print("--- SKETCHING ---")
seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

# 3. PREPARE SKETCH ASSAY
DefaultAssay(seurat_obj) <- "sketch"
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# 4. THE HAND-OFF (NUCLEAR OPTION)
print("--- EXPORTING SUBSET FOR ISOLATED INTEGRATION ---")
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))

# CRITICAL FIX: Strip the massive 'RNA' assay. Keep ONLY 'sketch' and 'pca'.
# This reduces the file size from ~250GB to ~2GB.
print("--- SLIMMING DOWN OBJECT (DietSeurat) ---")
sketch_subset <- DietSeurat(
  sketch_subset, 
  assays = "sketch", 
  dimreducs = "pca"
)

# Verify size before saving (Optional check for your logs)
print(format(object.size(sketch_subset), units = "GB"))

# Save with version 3
saveRDS(sketch_subset, "temp_sketch_only.rds", version = 3, compress = FALSE)
rm(sketch_subset)
gc()

print("--- LAUNCHING EXTERNAL R PROCESS (integrate_sketch.R) ---")
system("Rscript integrate_sketch.R")

# 5. RE-IMPORT AND MERGE
print("--- MERGING INTEGRATED RESULTS ---")
if(!file.exists("temp_sketch_integrated.rds")) {
  stop("Integration script failed to produce output.")
}
integrated_subset <- readRDS("temp_sketch_integrated.rds")
seurat_obj[["integrated.rpca"]] <- integrated_subset[["integrated.rpca"]]

# Cleanup temp files
file.remove("temp_sketch_only.rds")
file.remove("temp_sketch_integrated.rds")
rm(integrated_subset)
gc()

# 6. PROJECT & UMAP
print("--- PROJECTING TO FULL DATASET ---")
# Keep your current projection logic
seurat_obj <- ProjectData(
  object = seurat_obj,
  assay = "RNA",
  full.reduction = "pca",       
  sketched.assay = "sketch",
  sketched.reduction = "pca",    
  nmisc = 20
)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ProjectIntegration(
  object = seurat_obj, 
  sketched.assay = "sketch", 
  integration.reduction = "integrated.rpca", 
  new.reduction = "integrated.rpca.full"
)

print("--- RUNNING UMAP ---")
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction = "integrated.rpca.full", 
  dims = 1:30, 
  reduction.name = "umap.full", 
  verbose = FALSE
)

# 7. SAVE FINAL
print("--- SAVING FINAL OBJECT ---")
saveRDS(seurat_obj, "seurat_obj_integrated_sketch.rds", version = 3, compress = FALSE)
print("--- PIPELINE COMPLETE ---")