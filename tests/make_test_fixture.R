#!/usr/bin/env Rscript
# tests/make_test_fixture.R
# ---------------------------------------------------------------------------
# Generates a tiny, synthetic Seurat object shaped like the real
# 'aggregated_cellbender_filtered.rds' input this pipeline expects, so CI can
# smoke-test the QC + RPCA-integration stages (MPNST.R, MPNST_RPCA_Integration.R)
# without needing the real (multi-GB, non-public) scRNA-seq dataset.
#
# This is NOT real biological data -- counts are random Poisson draws. It
# exists purely to exercise the pipeline's code paths (I/O, QC filtering,
# normalization, PCA, RPCA integration) end to end.
#
# Usage: Rscript tests/make_test_fixture.R [output_path]
# ---------------------------------------------------------------------------

suppressMessages(library(Seurat))
suppressMessages(library(Matrix))

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else "aggregated_cellbender_filtered.rds"

# n_genes must comfortably clear MPNST.R's nFeature_RNA > 500 QC filter even
# after Poisson sampling noise -- see percent-nonzero math in the comment
# below.
n_genes <- 1000
n_batches <- 3
# MPNST_RPCA_Integration.R hardcodes k.weight = 100 for RPCAIntegration, which
# requires at least k.weight cells in the smallest batch/dataset being
# integrated, or Seurat errors out during anchor weighting. Keep a healthy
# margin above 100 per batch so the smoke test exercises that code path
# instead of just failing on cell count.
cells_per_batch <- 220
n_cells <- n_batches * cells_per_batch

# Gene names: a handful of mitochondrial genes (for percent.mt) + generic genes
mt_genes <- paste0("MT-", c("ND1", "ND2", "CO1", "CO2", "ATP6", "CYB"))
other_genes <- paste0("GENE", seq_len(n_genes - length(mt_genes)))
gene_names <- c(mt_genes, other_genes)

# Simulate counts: modest depth per cell, mito genes given a low mean so
# percent.mt stays well under the pipeline's 20% filter for most cells.
counts <- matrix(
  rpois(n_genes * n_cells, lambda = 2),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(gene_names, paste0("cell_", seq_len(n_cells)))
)
counts[seq_along(mt_genes), ] <- matrix(
  rpois(length(mt_genes) * n_cells, lambda = 0.5),
  nrow = length(mt_genes)
)
counts <- Matrix::Matrix(counts, sparse = TRUE)  # dgCMatrix, version-stable coercion

batch_ids <- rep(paste0("batch", seq_len(n_batches)), each = cells_per_batch)
lab_ids <- rep(c("LabA", "LabB", "LabC"), length.out = n_cells)

meta <- data.frame(
  orig.ident = batch_ids,
  lab = lab_ids,
  row.names = colnames(counts)
)

seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  min.cells = 0,
  min.features = 0
)

# MPNST.R explicitly assumes its input arrives with RNA layers already split
# per batch (Seurat v5 layer architecture) -- it deliberately does NOT call
# split() itself (see the "KEY FIX" comment in MPNST.R section 5). Mimic
# that here so RPCA/FastMNN integration has real per-batch layers to work
# with, matching the shape of the real aggregated_cellbender_filtered.rds.
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$orig.ident)

saveRDS(seurat_obj, out_path, compress = FALSE)
cat(sprintf(
  "--- TEST FIXTURE WRITTEN: %s (%d genes x %d cells, %d batches) ---\n",
  out_path, n_genes, n_cells, n_batches
))
