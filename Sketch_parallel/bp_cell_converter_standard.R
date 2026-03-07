#' @description Loads the BPCells package, providing functionalities for handling large-scale single-cell RNA-seq count matrices.
library(BPCells)

#' @title Open BPCells Matrix Directory
#' @description Opens an existing BPCells MatrixDir from the specified path, creating an `IterableMatrix` object.
#'   This object provides an efficient, memory-sparing way to work with large matrices stored on disk.
#' @details The path `/autofs/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/counts_bpcells`
#'   should point to a directory previously created and populated with BPCells matrix data.
#' @seealso \code{\link[BPCells]{open_matrix_dir}} for details on `IterableMatrix` objects.
mat <- open_matrix_dir("/autofs/projects-t3/lilab/vmenon/Zhu_MPNST/Sketch_parallel/counts_bpcells")

#' @title Export BPCells Matrix to HDF5 Format
#' @description Writes the opened `IterableMatrix` object to an HDF5 file.
#'   HDF5 is a universal format that allows the matrix data to be easily read by other
#'   programming languages and tools, such as Python's `h5py` or `scanpy`.
#' @param mat The `IterableMatrix` object to be exported. This is typically the result
#'   of `open_matrix_dir`.
#' @param file A character string specifying the name of the output HDF5 file (e.g., "counts.h5").
#' @param group A character string specifying the HDF5 group name under which the
#'   matrix data will be stored within the HDF5 file. This helps locate the data
#'   when reading from Python.
#' @param compress A logical value indicating whether to apply compression to the data
#'   within the HDF5 file. Setting to `TRUE` can significantly reduce file size.
#' @details The generated HDF5 file will contain the count matrix under the specified group name.
#'   For example, if `group = "counts"`, Python code can access it via `f['counts']`
#'   where `f` is an `h5py.File` object.
#' @seealso \code{\link[BPCells]{write_matrix_hdf5}} for more details on HDF5 export options.
write_matrix_hdf5(
  mat = mat,
  file = "counts.h5",
  group = "counts",
  compress = TRUE
)

#' @title Export Seurat Object Metadata and Feature/Barcode Names
#' @description Exports essential metadata, gene names (features), and cell barcode names
#'   from a Seurat object into separate CSV and TSV files. These files provide
#'   contextual information for the exported HDF5 matrix.
#' @param seurat_obj A Seurat object from which `meta.data`, `rownames`, and `colnames`
#'   will be extracted. It is assumed that this object is available in the current
#'   R session.
#' @details Metadata is exported as a standard CSV file, which is suitable for
#'   small to medium-sized tabular data. Gene names and cell barcodes are exported
#'   as simple plain text files (TSV-like, one item per line) for easy parsing.
#' @section Files Created:
#' \itemize{
#'   \item `metadata.csv`: Contains cell-level metadata (e.g., clustering, UMAP coordinates)
#'     from `seurat_obj@meta.data`. Row names are retained.
#'   \item `genes.tsv`: Contains the feature (gene) identifiers, derived from `rownames(seurat_obj)`.
#'   \item `barcodes.tsv`: Contains the cell barcode identifiers, derived from `colnames(seurat_obj)`.
#' }
#' @seealso \code{\link[utils]{write.csv}}, \code{\link[base]{write}}, \code{\link[SeuratObject]{Seurat}}
write.csv(seurat_obj@meta.data, "metadata.csv", row.names = TRUE)
write(rownames(seurat_obj), "genes.tsv")
write(colnames(seurat_obj), "barcodes.tsv")

#' @description Prints a confirmation message to the console indicating that the data
#'   export process has been successfully completed and the `counts.h5` file has been created.
print("--- EXPORT COMPLETE: counts.h5 created ---")
