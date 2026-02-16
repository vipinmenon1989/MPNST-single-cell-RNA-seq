library(Seurat)
library(ggplot2)
library(parallel)

# 1. SETUP
# -------------------------------------------------------------------------
CORES <- 4  # Adjust based on RAM
print(paste0("--- HPC MODE: PARALLEL PLOTTING (", CORES, " CORES) ---"))

# 2. LOAD & REPAIR OBJECTS (SEQUENTIAL)
# -------------------------------------------------------------------------
print("--- LOADING OBJECTS ---")
pbj <- readRDS("seurat_obj_pre_integration.rds")      
xbj <- readRDS("seurat_obj_mnn_integrated_1M.rds")    

# --- REPAIR UNINTEGRATED (pbj) ---
if ("umap.unintegrated" %in% names(pbj@reductions)) {
  red_unint <- "umap.unintegrated"
  print("Status: Found 'umap.unintegrated' in pbj.")
} else if ("umap" %in% names(pbj@reductions)) {
  red_unint <- "umap"
  print("Status: Found default 'umap' in pbj.")
} else {
  print("--- MISSING UMAP IN PBJ. CALCULATING NOW... ---")
  pbj <- RunUMAP(pbj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", verbose = FALSE)
  red_unint <- "umap.unintegrated"
  print("--- CALCULATION COMPLETE ---")
}

# --- REPAIR INTEGRATED (xbj) ---
if ("umap.mnn" %in% names(xbj@reductions)) {
  red_int <- "umap.mnn"
  print("Status: Found 'umap.mnn' in xbj.")
} else {
  print("--- MISSING UMAP IN XBJ. CALCULATING NOW... ---")
  xbj <- RunUMAP(xbj, reduction = "integrated.mnn.full", dims = 1:30, reduction.name = "umap.mnn", verbose = FALSE)
  red_int <- "umap.mnn"
  print("--- CALCULATION COMPLETE ---")
}

# 3. DEFINE JOBS
# -------------------------------------------------------------------------
features <- c("wu_annot", "suppiah_annot", "lab", "tumor_classification", "orig.ident")
jobs <- list()

for (f in features) {
  # Add Unintegrated Job
  if (f %in% colnames(pbj@meta.data)) {
    jobs[[paste0("Unint_", f)]] <- list(
      type = "Unintegrated",
      obj = pbj, 
      feature = f, 
      reduction = red_unint
    )
  }
  # Add Integrated Job
  if (f %in% colnames(xbj@meta.data)) {
    jobs[[paste0("Int_", f)]] <- list(
      type = "Integrated",
      obj = xbj, 
      feature = f, 
      reduction = red_int
    )
  }
}

print(paste0("--- QUEUED ", length(jobs), " PLOTS ---"))

# 4. EXECUTE PARALLEL PLOTTING (WITH LABEL LOGIC)
# -------------------------------------------------------------------------
generate_plot <- function(job) {
  tryCatch({
    file_name <- paste0(job$type, "_", job$feature, ".png")
    
    # Check if reduction actually exists inside the worker
    if (!job$reduction %in% names(job$obj@reductions)) {
      return(paste0("FAILED: Reduction ", job$reduction, " lost in forking."))
    }

    # LOGIC: Turn on labeling ONLY if feature is NOT "orig.ident"
    should_label <- (job$feature != "orig.ident")

    p <- DimPlot(job$obj, 
                 reduction = job$reduction, 
                 group.by = job$feature, 
                 label = should_label,      # TRUE for annotations, FALSE for batch ID
                 label.size = 4,            # Clean, readable font size
                 repel = TRUE,              # Pushes labels apart so they don't overlap
                 raster = FALSE) + 
         ggtitle(paste0(job$type, ": ", job$feature)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         NoLegend()
    
    ggsave(file_name, plot = p, width = 12, height = 10, dpi = 300)
    return(paste0("DONE: ", file_name))
    
  }, error = function(e) {
    return(paste0("ERROR on ", job$feature, ": ", e$message))
  })
}

results <- mclapply(jobs, generate_plot, mc.cores = CORES, mc.preschedule = FALSE)

# 5. FINAL STATUS
# -------------------------------------------------------------------------
print("--- RESULTS ---")
print(unlist(results))
