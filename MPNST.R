library(Seurat)
library(ggplot2)
library(patchwork)
library(future)
library(dplyr)

# 1. HPC & PARALLEL SETUP
# -------------------------------------------------------------------------
# We use 'multicore' for Linux HPCs. 
plan("multicore", workers = 4) 
# Set max memory to 200GB (or 300GB if needed)
options(future.globals.maxSize = 300 * 1024^3)
print("--- HPC MODE: 4 CORES ENABLED ---")

# 2. LOAD DATA
# -------------------------------------------------------------------------
file_path <- "aggregated_cellbender_filtered.rds"
batch_column <- "orig.ident" 

seurat_obj <- readRDS(file_path)
print("--- DATA LOADED ---")

# 3. PRE-FILTERING (UNIFORM BLACK PLOTS)
# -------------------------------------------------------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

print("--- GENERATING PRE-QC PLOTS ---")

# Plotting using ggplot2 for uniform black color and high contrast
p1 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=percent.mt)) + 
      geom_point(color="black", size=0.5) + 
      theme_classic() +
      ggtitle("Pre-Filter: Counts vs Mito%") + 
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

p2 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + 
      geom_point(color="black", size=0.5) + 
      theme_classic() +
      ggtitle("Pre-Filter: Counts vs Genes") + 
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

ggsave("01_QC_Pre_Scatter_Uniform.png", plot = p1 + p2, width = 16, height = 8, dpi = 300)

# --- CSV SUMMARY (PRE) - MODIFIED TO INCLUDE 'lab' ---
print("--- GENERATING PRE-FILTER CSV (WITH LAB INFO) ---")

# We group by the batch column AND 'lab'
pre_stats <- seurat_obj@meta.data %>%
  group_by(.data[[batch_column]], lab) %>%
  summarise(
    Num_Cells = n(),
    Median_nFeature = median(nFeature_RNA),
    Median_nCount = median(nCount_RNA),
    Median_MT = median(percent.mt),
    .groups = 'drop'
  )

# Create TOTAL row with a placeholder for 'lab'
total_pre <- data.frame(
  batch_column = "TOTAL",
  lab = "ALL",  # Placeholder for the total row
  Num_Cells = ncol(seurat_obj),
  Median_nFeature = median(seurat_obj$nFeature_RNA),
  Median_nCount = median(seurat_obj$nCount_RNA),
  Median_MT = median(seurat_obj$percent.mt)
)
# Rename the first column to match the batch column name (orig.ident)
colnames(total_pre)[1] <- batch_column 

# Combine and Save
write.csv(rbind(pre_stats, total_pre), "01_QC_Pre_Filtering_Stats.csv", row.names = FALSE)
print("--- PRE-FILTER METRICS SAVED ---")


# 4. FILTERING
# -------------------------------------------------------------------------
print("--- FILTERING DATA ---")
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 500 & 
                              nFeature_RNA < 9000 & 
                              percent.mt < 20)
gc() # Free up memory
print(paste("Cells remaining:", ncol(seurat_obj)))

# --- POST-FILTERING PLOTS ---
print("--- GENERATING POST-QC PLOTS ---")

p3 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=percent.mt)) + 
      geom_point(color="black", size=0.5) + 
      theme_classic() +
      ggtitle("Post-Filter: Counts vs Mito%") + 
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

p4 <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + 
      geom_point(color="black", size=0.5) + 
      theme_classic() +
      ggtitle("Post-Filter: Counts vs Genes") + 
      theme(axis.text = element_text(size=12), axis.title = element_text(size=14))

ggsave("02_QC_Post_Scatter_Uniform.png", plot = p3 + p4, width = 16, height = 8, dpi = 300)

# --- CSV SUMMARY (POST) - MODIFIED TO INCLUDE 'lab' ---
print("--- GENERATING POST-FILTER CSV (WITH LAB INFO) ---")

post_stats <- seurat_obj@meta.data %>%
  group_by(.data[[batch_column]], lab) %>%
  summarise(
    Num_Cells = n(),
    Median_nFeature = median(nFeature_RNA),
    Median_nCount = median(nCount_RNA),
    Median_MT = median(percent.mt),
    .groups = 'drop'
  )

total_post <- data.frame(
  batch_column = "TOTAL",
  lab = "ALL",
  Num_Cells = ncol(seurat_obj),
  Median_nFeature = median(seurat_obj$nFeature_RNA),
  Median_nCount = median(seurat_obj$nCount_RNA),
  Median_MT = median(seurat_obj$percent.mt)
)
colnames(total_post)[1] <- batch_column

write.csv(rbind(post_stats, total_post), "02_QC_Post_Filtering_Stats.csv", row.names = FALSE)
print("--- POST-FILTER METRICS SAVED ---")


# 5. NORMALIZE & FIND FEATURES (NO SPLIT/JOIN)
# -------------------------------------------------------------------------
print("--- STARTING NORMALIZATION ---")

# KEY FIX: We do NOT JoinLayers (causes overflow). 
# We do NOT split (already split).
# We just run NormalizeData, which automatically handles the existing split layers.

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)


# 6. SCALE & PCA
# -------------------------------------------------------------------------
print("--- SCALING DATA (HVG ONLY) ---")
# Scaling works on the split layers automatically in V5
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

print("--- RUNNING PCA ---")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30, verbose = FALSE)

print("--- SAVING CHECKPOINT RDS (NO COMPRESSION) ---")
# Keep compress = FALSE to prevent crash
saveRDS(seurat_obj, "seurat_obj_pre_integration.rds", compress = FALSE)

print("--- SCRIPT COMPLETED ---")