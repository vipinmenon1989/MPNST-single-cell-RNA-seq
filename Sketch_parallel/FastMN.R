library(Seurat)
library(SeuratWrappers) # Required for FastMNNIntegration
library(batchelor)      # The underlying engine for FastMNN
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC SETUP
# -------------------------------------------------------------------------
# Using 'multicore' for Linux HPC; use 'multisession' if on Windows.
plan("sequential") 
options(future.globals.maxSize = 400 * 1024^3) # 400GB limit for globals

print("--- HPC MODE: INITIALIZED ---")

# 2. LOAD DATA
# -------------------------------------------------------------------------
file_path <- "seurat_obj_pre_integration.rds" 

if (!file.exists(file_path)) {
  stop("Input file not found.")
}

seurat_obj <- readRDS(file_path)
print(paste("--- DATA LOADED: ", ncol(seurat_obj), " cells ---"))

# 3. GEOMETRIC SKETCHING
# -------------------------------------------------------------------------
# We sketch 50k cells to represent the 1M population. 
# This is the "Blunt Truth" for speed: 1M cell integration is unnecessary 
# if a 50k sketch captures the manifold.
print("--- RUNNING GEOMETRIC SKETCH (50k cells) ---")
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

# 4. MEMORY EVICTION & SUBSETTING
# -------------------------------------------------------------------------
print("--- SUBSETTING SKETCH & EVICTING FULL OBJECT ---")
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))
DefaultAssay(sketch_subset) <- "sketch"

# Save the full object to disk and WIPE from RAM
saveRDS(seurat_obj, "temp_seurat_full.rds", compress = FALSE)
rm(seurat_obj)
gc() # Deep garbage collection

# 5. PROCESS SKETCH
# -------------------------------------------------------------------------
print("--- PROCESSING SKETCH SUBSET ---")
sketch_subset <- FindVariableFeatures(sketch_subset, verbose = FALSE)
sketch_subset <- ScaleData(sketch_subset, verbose = FALSE)
sketch_subset <- RunPCA(sketch_subset, verbose = FALSE)

# 6. PARALLEL FastMNN INTEGRATION
# -------------------------------------------------------------------------
# Switching to multicore NOW because the environment is empty of the 1M cell object.
print("--- HPC MODE: MULTICORE (8 WORKERS) ---")
plan("multicore", workers = 8)

print("--- INTEGRATING LAYERS VIA FastMNN ---")
# FastMNN integrates in the PCA space; it is faster and more stable than RPCA for high-K batches.
sketch_subset <- IntegrateLayers(
  object = sketch_subset,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = TRUE
)

# 7. RE-LINKING & PROJECTION
# -------------------------------------------------------------------------
plan("sequential")
print("--- RELOADING MAIN OBJECT & PROJECTING ---")
seurat_obj <- readRDS("temp_seurat_full.rds")

# Transfer integration from sketch to full
seurat_obj[["integrated.mnn"]] <- sketch_subset[["integrated.mnn"]]
rm(sketch_subset)
gc()

seurat_obj <- ProjectIntegration(
  object = seurat_obj, 
  sketched.assay = "sketch", 
  reduction = "integrated.mnn",           # Changed from integration.reduction
  reduction.name = "integrated.mnn.full"  # Changed from new.reduction
)

# 8. FINAL DIM REDUCTION & SAVE
# -------------------------------------------------------------------------ls
print("--- FINAL UMAP FOR 1M CELLS ---")
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction = "integrated.mnn.full", 
  dims = 1:30, 
  reduction.name = "umap.mnn", 
  verbose = FALSE
)

print("--- SAVING FINAL INTEGRATED OBJECT ---")
saveRDS(seurat_obj, "seurat_obj_mnn_integrated_1M.rds", compress = FALSE)
file.remove("temp_seurat_full.rds") 

print("--- PIPELINE COMPLETED ---")
