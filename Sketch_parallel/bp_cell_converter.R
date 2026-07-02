library(BPCells)

# 1. Open the MatrixDir you already created (see r_metadata.R, which writes
#    this directory as a relative path so the pipeline stays portable across
#    machines/checkouts)
mat <- open_matrix_dir("counts_bpcells")

# 2. Write to HDF5 (The Universal Format)
# This will create a file that Python's 'h5py' or 'scanpy' can read.
# We name the group "counts" so we know where to find it in Python.
# NOTE: the BPCells argument is 'path', not 'file' (see bp.R) -- fixed below.
write_matrix_hdf5(
  mat = mat, 
  path = "counts.h5", 
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
