library(Seurat)

# 1. RELOAD THE OBJECT (Just for metadata)
print("--- RELOADING SEURAT OBJECT FOR METADATA ---")
seurat_obj <- readRDS("seurat_obj_pre_integration.rds")

# 2. EXPORT METADATA & COORDINATES
# These are small files and won't trigger the 32-bit integer error
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- METADATA EXPORT COMPLETE ---")
