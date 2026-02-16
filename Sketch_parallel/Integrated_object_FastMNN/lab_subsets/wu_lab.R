library(Seurat)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)
library(patchwork)
options(future.globals.maxSize = 300 * 1024^3)

# ==============================================================================
# 1. EXTRACT & CLEAN THE SUBSET
# ==============================================================================
message("Loading full object...")
obj_full <- readRDS("lab_Wu.rds")

message("Subsetting for Lab 'Wu'...")
# Verify 'lab' is the correct column name. If it's 'orig.ident', change it here.
# We also filter out cells where 'wu_annot' is NA to ensure clean validation later.
obj <- subset(obj_full, subset = lab == "Wu" & (wu_annot!= "unknown"))

message(paste0("Subset complete. Remaining cells: ", ncol(obj)))

# CLEAN SLATE PROTOCOL: Remove old reductions and graphs
# We want the clustering to be driven ONLY by this subset's variance
obj@reductions <- list()
obj@graphs <- list()

# ==============================================================================
# 2. RE-PROCESS (THE CRITICAL STEP)
# ==============================================================================
message("Re-processing the Wu subset...")

# We normalize and find variable features SPECIFICALLY for this subset.
# This ensures we pick up markers relevant to these cells, not the 1M dataset.
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

# ==============================================================================
# 3. MULTI-RESOLUTION CLUSTERING LOOP
# ==============================================================================
message("Clustering resolutions 0.1 - 1.0...")
resolutions <- seq(0.1, 1.0, by = 0.1)

# Explicitly name the graph to prevent Seurat v5 errors
graph_name <- "wu_subset"
obj <- FindNeighbors(obj, dims = 1:30, graph.name = c(paste0(graph_name, "_nn"), paste0(graph_name, "_snn")))
obj <- FindClusters(obj, graph.name = paste0(graph_name, "_snn"), resolution = resolutions, verbose = FALSE)

# ==============================================================================
# 4. STATISTICAL MATCHING (ARI)
# ==============================================================================
message("Calculating ARI agreement with 'wu_annot'...")

# Identify the new resolution columns (e.g., 'wu_subset_snn_res.0.5')
res_columns <- grep(paste0(graph_name, "_snn_res\\.\\d"), colnames(obj@meta.data), value = TRUE)
res_columns <- res_columns[order(as.numeric(gsub(".*_res\\.", "", res_columns)))]

ari_scores <- numeric(length(res_columns))
cleaned_resolutions <- as.numeric(gsub(".*_res\\.", "", res_columns))

for (i in seq_along(res_columns)) {
  # Calculate ARI (we already removed NAs in step 1, so this is safe)
  ari_scores[i] <- adjustedRandIndex(obj$wu_annot, obj@meta.data[[res_columns[i]]])
}

validation_df <- data.frame(Resolution = cleaned_resolutions, 
                            ColumnName = res_columns, 
                            ARI = ari_scores)

best_res_row <- validation_df[which.max(validation_df$ARI), ]
best_res <- best_res_row$Resolution
best_col <- best_res_row$ColumnName

message(paste0("--------------------------------------------------"))
message(paste0(" Winning Resolution: ", best_res))
message(paste0(" Peak ARI: ", round(best_res_row$ARI, 3)))
message(paste0(" (An ARI > 0.4 indicates strong biological alignment)"))
message(paste0("--------------------------------------------------"))

# ==============================================================================
# 5. OUTPUT PLOTS (HPC SAFE)
# ==============================================================================

# PDF 1: The Validation Metrics
pdf("Wu_Subset_Validation.pdf", width = 12, height = 6)

# Elbow Plot
p1 <- ggplot(validation_df, aes(x = Resolution, y = ARI)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Which Resolution matches Wu Annot?", y = "ARI Score") +
  theme_minimal()

# Confusion Matrix Heatmap
conf_data <- obj@meta.data %>%
  group_by(wu_annot, !!sym(best_col)) %>%
  tally() %>%
  group_by(wu_annot) %>%
  mutate(Percent = n / sum(n))

p2 <- ggplot(conf_data, aes(x = !!sym(best_col), y = wu_annot, fill = Percent)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#D55E00") +
  labs(title = paste("Best Fit Heatmap (Res", best_res, ")"), x = "Cluster", y = "Annotation") +
  theme_minimal()

print(p1 + p2)
dev.off()

# PDF 2: The Side-by-Side UMAP
pdf("Wu_Subset_UMAP.pdf", width = 14, height = 7)
p3 <- DimPlot(obj, group.by = "wu_annot", label = TRUE, raster=FALSE) + ggtitle("Ground Truth")
p4 <- DimPlot(obj, group.by = best_col, label = TRUE, raster=FALSE) + ggtitle(paste("Best Clustering (Res", best_res, ")"))
print(p3 + p4)
dev.off()
