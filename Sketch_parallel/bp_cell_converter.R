library(BPCells)

# 1. Open the MatrixDir you already created
mat <- open_matrix_dir("/autofs/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/counts_bpcells")

# 2. Write to HDF5 (The Universal Format)
# This will create a file that Python's 'h5py' or 'scanpy' can read.
# We name the group "counts" so we know where to find it in Python.
write_matrix_hdf5(
  mat = mat, 
  file = "counts.h5", 
  group = "counts",
  compress = TRUE
)

# 3. Export Metadata (We still need these as CSVs)
# Metadata is small, so standard write.csv is fine.
# Assuming you still have your seurat_obj in this session
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

print("--- EXPORT COMPLETE: counts.h5 created ---")
