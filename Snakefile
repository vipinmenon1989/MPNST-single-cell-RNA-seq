"""
MPNST scRNA-seq Snakemake workflow
====================================
Wraps the R/Seurat pipeline in this repo (QC -> batch integration -> optional
FastMNN clustering/annotation) into a reproducible, resumable Snakemake DAG.

The underlying R scripts each read/write fixed filenames relative to their
own working directory (this mirrors how they're launched by the .slurm job
scripts in this repo). Snakemake rules below `cd` into the right directory
and stage inputs under the filename each script expects via a symlink, so
none of the original R scripts had to be modified to work with Snakemake.

Configure via config/config.yaml, then:

    snakemake --cores 1 -n              # preview the DAG
    snakemake --cores 1                 # run with deps already on PATH
    snakemake --cores 1 --use-conda     # run using envs/environment.yaml

NOTE: this pipeline processes real scRNA-seq Seurat objects (hundreds of
GB, 30-40+ hr runtimes per the README) and requires R + Seurat 5 + several
Bioconductor/CRAN packages (see envs/environment.yaml). It is not something
that can be smoke-tested against the full dataset in CI; CI instead runs
the DAG against a tiny synthetic fixture (see tests/make_test_fixture.R).
"""

configfile: "config/config.yaml"

INPUT_RDS = config["input_rds"]
INTEGRATION = config.get("integration", "rpca")
RUN_CLUSTERING = config.get("run_clustering", False)

VALID_INTEGRATIONS = ("rpca", "sketch", "sketch_parallel", "fastmnn")
if INTEGRATION not in VALID_INTEGRATIONS:
    raise ValueError(
        f"config['integration'] must be one of {VALID_INTEGRATIONS}, got '{INTEGRATION}'"
    )

SKETCH_DIR = "Sketch_parallel"
PROCESSED_DIR = "Sketch_parallel/Integrated_object_FastMNN/processed_reads"

# Per-integration-mode outputs that `rule all` should request.
_INTEGRATION_TARGETS = {
    "rpca": ["seurat_obj_integrated_rpca.rds"],
    "sketch": ["seurat_obj_integrated_sketch.rds"],
    "sketch_parallel": [f"{SKETCH_DIR}/seurat_obj_integrated_sketch.rds"],
    "fastmnn": [f"{SKETCH_DIR}/seurat_obj_mnn_integrated_1M.rds"],
}

_CLUSTERING_TARGETS = [
    f"{PROCESSED_DIR}/seurat_obj_clustered_1M.rds",
    f"{PROCESSED_DIR}/resolutions_comparison_final.pdf",
]


def _all_targets():
    targets = ["seurat_obj_pre_integration.rds"] + _INTEGRATION_TARGETS[INTEGRATION]
    if RUN_CLUSTERING:
        if INTEGRATION != "fastmnn":
            raise ValueError(
                "config['run_clustering'] = true requires config['integration'] = 'fastmnn' "
                "(the clustering stage consumes seurat_obj_mnn_integrated_1M.rds)."
            )
        targets += _CLUSTERING_TARGETS
    return targets


rule all:
    input:
        _all_targets(),


# ---------------------------------------------------------------------------
# Stage 1: QC & preprocessing (repo root)
# ---------------------------------------------------------------------------
rule qc_preprocessing:
    """MPNST.R: load raw counts, QC filter, normalize, scale, PCA."""
    input:
        rds=INPUT_RDS,
    output:
        rds="seurat_obj_pre_integration.rds",
        pre_scatter="01_QC_Pre_Scatter_Uniform.png",
        pre_stats="01_QC_Pre_Filtering_Stats.csv",
        post_scatter="02_QC_Post_Scatter_Uniform.png",
        post_stats="02_QC_Post_Filtering_Stats.csv",
    log:
        "logs/qc_preprocessing.log",
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p logs
        if [ "{input.rds}" != "aggregated_cellbender_filtered.rds" ]; then
            ln -sf "$(realpath {input.rds})" aggregated_cellbender_filtered.rds
        fi
        Rscript MPNST.R > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Stage 2: batch integration (choice of strategy, see config['integration'])
# ---------------------------------------------------------------------------
rule integrate_rpca:
    """MPNST_RPCA_Integration.R: canonical RPCA integration (repo root)."""
    input:
        rds="seurat_obj_pre_integration.rds",
    output:
        rds="seurat_obj_integrated_rpca.rds",
        unintegrated="03_UMAP_Unintegrated.png",
        integrated="04_UMAP_Integrated_RPCA.png",
        comparison="05_Integration_Comparison.png",
    log:
        "logs/integrate_rpca.log",
    conda:
        "envs/environment.yaml"
    shell:
        "mkdir -p logs && Rscript MPNST_RPCA_Integration.R > {log} 2>&1"


rule integrate_sketch:
    """Sketch.R: single-process 50k-cell sketch + RPCA (repo root)."""
    input:
        rds="seurat_obj_pre_integration.rds",
    output:
        rds="seurat_obj_integrated_sketch.rds",
        plot="04_UMAP_Integrated_Sketch.png",
    log:
        "logs/integrate_sketch.log",
    conda:
        "envs/environment.yaml"
    shell:
        "mkdir -p logs && Rscript Sketch.R > {log} 2>&1"


rule integrate_sketch_parallel:
    """Sketch_parallel/main_pipeline.R: subprocess hand-off sketch + RPCA."""
    input:
        rds="seurat_obj_pre_integration.rds",
    output:
        rds=f"{SKETCH_DIR}/seurat_obj_integrated_sketch.rds",
    log:
        "logs/integrate_sketch_parallel.log",
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p logs {SKETCH_DIR}/logs
        ln -sf "$(realpath {input.rds})" {SKETCH_DIR}/seurat_obj_pre_integration.rds
        cd {SKETCH_DIR} && Rscript main_pipeline.R > ../logs/integrate_sketch_parallel.log 2>&1
        """


rule integrate_fastmnn:
    """Sketch_parallel/FastMN.R: sketch + FastMNN integration for 1M cells."""
    input:
        rds="seurat_obj_pre_integration.rds",
    output:
        rds=f"{SKETCH_DIR}/seurat_obj_mnn_integrated_1M.rds",
    log:
        "logs/integrate_fastmnn.log",
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p logs {SKETCH_DIR}/logs
        ln -sf "$(realpath {input.rds})" {SKETCH_DIR}/seurat_obj_pre_integration.rds
        cd {SKETCH_DIR} && Rscript FastMN.R > ../logs/integrate_fastmnn.log 2>&1
        """


# ---------------------------------------------------------------------------
# Stage 3 (optional, fastmnn only): clustering + marker detection
# ---------------------------------------------------------------------------
rule cluster_fastmnn:
    """processing.R: multi-resolution clustering + marker tables."""
    input:
        rds=f"{SKETCH_DIR}/seurat_obj_mnn_integrated_1M.rds",
    output:
        rds=f"{PROCESSED_DIR}/seurat_obj_clustered_1M.rds",
        grid=f"{PROCESSED_DIR}/umap_outputs/UMAP_Grid_All_Res.pdf",
    log:
        "logs/cluster_fastmnn.log",
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p logs {PROCESSED_DIR}
        ln -sf "$(realpath {input.rds})" {PROCESSED_DIR}/seurat_obj_mnn_integrated_1M.rds
        cd {PROCESSED_DIR} && Rscript processing.R > ../../../../logs/cluster_fastmnn.log 2>&1
        """


rule annotate_fastmnn:
    """annotate.R: resolution-vs-annotation comparison plots."""
    input:
        rds=f"{PROCESSED_DIR}/seurat_obj_clustered_1M.rds",
    output:
        pdf=f"{PROCESSED_DIR}/resolutions_comparison_final.pdf",
    log:
        "logs/annotate_fastmnn.log",
    conda:
        "envs/environment.yaml"
    shell:
        "mkdir -p logs && cd {PROCESSED_DIR} && Rscript annotate.R > ../../../../logs/annotate_fastmnn.log 2>&1"
