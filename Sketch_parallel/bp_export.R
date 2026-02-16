library(BPCells)
library(Matrix)

# 1. Open your H5 matrix
mat <- open_matrix_hdf5(path = "counts.h5", group = "counts")

# 2. Stream to Matrix Market (.mtx)
# This writes directly to disk and is the most robust bridge to Python.
writeMM(obj = mat, file = "counts.mtx")

# 3. Export Metadata (Ensuring they are fresh)
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- EXPORT SUCCESSFUL: counts.mtx created ---")
