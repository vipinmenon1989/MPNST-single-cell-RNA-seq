library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC SETUP (INITIAL)
# -------------------------------------------------------------------------
# Start sequential to load data safely. We will toggle this later.
plan("sequential") 
options(future.globals.maxSize = 100 * 1024^3)

print("--- HPC MODE: INITIALIZED ---")

# 2. LOAD DATA
# -------------------------------------------------------------------------
file_path <- "seurat_obj_pre_integration.rds" 
batch_column <- "orig.ident"

if (!file.exists(file_path)) {
  stop("Input file not found. Ensure Part 1 ran successfully.")
}

seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

# 3. SKETCH THE DATASET
# -------------------------------------------------------------------------
print("--- RUNNING GEOMETRIC SKETCH ---")
DefaultAssay(seurat_obj) <- "RNA"

# Note: Ignore 'nsketch' warnings.
seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

print("--- SKETCH COMPLETE ---")


# 4. OPTIMIZED INTEGRATION (DETACH & PARALLELIZE)
# -------------------------------------------------------------------------
print("--- OPTIMIZATION: DETACHING SKETCH FOR PARALLEL PROCESSING ---")

# A. Create a lightweight subset (Drops the 1M RNA cells from memory context)
# This prevents the "304GB export" error.
sketch_subset <- subset(seurat_obj, cells = Cells(seurat_obj[["sketch"]]))
DefaultAssay(sketch_subset) <- "sketch"

# B. Process the subset
print("--- PROCESSING SKETCH SUBSET ---")
sketch_subset <- FindVariableFeatures(sketch_subset, verbose = FALSE)
sketch_subset <- ScaleData(sketch_subset, verbose = FALSE)
sketch_subset <- RunPCA(sketch_subset, verbose = FALSE)

# C. Enable Parallelization
# We requested 8 CPUs. We use 'multicore' for Linux stability.
print("--- HPC MODE: SWITCHING TO MULTICORE (8 WORKERS) ---")
plan("multicore", workers = 8)

print("--- INTEGRATING SKETCH LAYERS (RPCA - PARALLEL) ---")
# This will now run on 8 cores and should complete in < 45 mins.
sketch_subset <- IntegrateLayers(
  object = sketch_subset, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

# D. Revert to Sequential & Reattach
plan("sequential")
print("--- HPC MODE: REVERTED TO SEQUENTIAL ---")

print("--- TRANSFERRING INTEGRATION TO MAIN OBJECT ---")
# We manually copy the calculated reduction back to the master object
seurat_obj[["integrated.rpca"]] <- sketch_subset[["integrated.rpca"]]

# Clean up to free RAM
rm(sketch_subset)
gc()

print("--- SKETCH INTEGRATION COMPLETE ---")


# 5. PROJECT INTEGRATION
# -------------------------------------------------------------------------
print("--- PROJECTING INTEGRATION TO FULL DATASET ---")
# This maps the 950k other cells onto the 50k anchors.

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- ProjectIntegration(
  object = seurat_obj, 
  sketched.assay = "sketch", 
  integration.reduction = "integrated.rpca", 
  new.reduction = "integrated.rpca.full"
)

print("--- PROJECTION COMPLETE ---")


# 6. RUN UMAP ON FULL DATASET
# -------------------------------------------------------------------------
print("--- RUNNING UMAP ON PROJECTED DATA ---")
# We use the full projected reduction 'integrated.rpca.full'
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction = "integrated.rpca.full", 
  dims = 1:30, 
  reduction.name = "umap.full", 
  verbose = FALSE
)

print("--- GENERATING INTEGRATED PLOT ---")
p2 <- DimPlot(seurat_obj, 
              reduction = "umap.full", 
              group.by = batch_column, 
              label = TRUE, 
              label.size = 2, 
              repel = TRUE, 
              raster = FALSE) + 
      NoLegend
