library(Seurat)
library(BPCells)

# 1. LOAD DATA
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

# 2. WRITE MATRIX TO DISK (BPCells handles v5 layers automatically)
# This creates a folder 'counts_bpcells' with the 1M cell matrix
write_matrix_dir(
  mat = seurat_obj[["RNA"]]$counts, 
  dir = "counts_bpcells",
  overwrite = TRUE
)

# 3. EXPORT METADATA & GENES
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- R EXPORT COMPLETE: Matrix stored in 'counts_bpcells' ---")
