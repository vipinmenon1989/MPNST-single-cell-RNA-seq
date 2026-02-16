library(BPCells)

# 1. Open the MatrixDir
mat <- open_matrix_dir("/autofs/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/counts_bpcells")

# 2. Corrected HDF5 Write
# In the latest version, the argument is 'path', not 'file'
write_matrix_hdf5(
  mat = mat, 
  path = "counts.h5", 
  group = "counts",
  compress = TRUE
)

# 3. Ensure Metadata is exported for Python
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- R EXPORT SUCCESSFUL: counts.h5 created ---")
