library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC SETUP (SEQUENTIAL MODE)
# -------------------------------------------------------------------------
# CRITICAL: We force sequential processing. 
# This bypasses the "304.77 GiB" export error completely.
plan("sequential") 

# We set the limit high just in case, but sequential mode shouldn't trigger it.
options(future.globals.maxSize = 500 * 1024^3)

print("--- HPC MODE: SEQUENTIAL (1 CORE) ---")

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
print("--- RUNNING GEOMETRIC SKETCH (SEQUENTIAL) ---")
# 1 Core is sufficient here.

DefaultAssay(seurat_obj) <- "RNA"

# Note: Ignore 'nsketch' warnings.
seurat_obj <- SketchData(
  object = seurat_obj, 
  ncells = 50000, 
  method = "LeverageScore",  
  sketched.assay = "sketch"
)

print("--- SKETCH COMPLETE ---")


# 4. INTEGRATE THE SKETCH
# -------------------------------------------------------------------------
print("--- PROCESSING SKETCH LAYERS ---")
# We switch to the sketch assay (50k cells)
DefaultAssay(seurat_obj) <- "sketch"

# Basic processing on the sketch
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

print("--- INTEGRATING SKETCH LAYERS (RPCA - SEQUENTIAL) ---")
# On 1 core, integrating 50k cells takes about 30-60 mins.
seurat_obj <- IntegrateLayers(
  object = seurat_obj, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

print("--- SKETCH INTEGRATION COMPLETE ---")


# 5. PROJECT INTEGRATION
# -------------------------------------------------------------------------
print("--- PROJECTING INTEGRATION TO FULL DATASET ---")
# This maps the 950k other cells onto the 50k anchors.
# This uses Matrix Algebra, which R handles efficiently on 1 core.

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
      NoLegend() + 
      ggtitle("Integrated (Sketch-Projected)") +
      theme(plot.title = element_text(hjust = 0.5))

ggsave("04_UMAP_Integrated_Sketch.png", plot = p2, width = 16, height = 12, dpi = 300)
print("--- INTEGRATED UMAP SAVED ---")


# 7. SAVE FINAL CHECKPOINT
# -------------------------------------------------------------------------
print("--- SAVING INTEGRATED OBJECT (NO COMPRESSION) ---")
# Always compress=FALSE for large objects on HPC
saveRDS(seurat_obj, "seurat_obj_integrated_sketch.rds", compress = FALSE)

print("--- PIPELINE COMPLETED ---")