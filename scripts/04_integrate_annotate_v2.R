# ══════════════════════════════════════════════════════════════════════════════
# 04_integrate_annotate.R
# Filter, combine, integrate (SCTransform), annotate with Azimuth
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Fix here() root before anything else ────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
PROJECT_ROOT<- if (length(args) >= 1) args[1] else getwd()
if (PROJECT_ROOT == "test"){
PROJECT_ROOT <- "/g/data/pq84/sarahhp/projects/scauto/test_data"}

if (all(!file.exists(file.path(PROJECT_ROOT, ".here")),
        !grepl("\\.Rproj",list.files(PROJECT_ROOT)))) {
  stop("project_root (", PROJECT_ROOT, ") doesn't look right — ",
       "expected to find .Rproj or .here file in dir. ",
       "Did you cd into the project directory before running?")
}

setwd(PROJECT_ROOT)
try(detach("package:here", unload = TRUE), silent = TRUE)
library(here)

if (here() != PROJECT_ROOT) {
  stop(glue::glue(
    "here() still wrong: {here()}\n",
    "Start R from: {PROJECT_ROOT}"
  ))
}
cat("✔ Project root:", here(), "\n")

utility_file <- "/g/data/pq84/sarahhp/projects/scauto/utility_seuratv5.R" # new utility file
if (!file.exists(utility_file)) stop(glue("utility.R not found: {utility_file}"))

options(future.globals.maxSize = 32 * 1024^3)  # 8 GiB, adjust to data size

# ── 1. Libraries ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(scuttle)
  library(Seurat)
  library(Azimuth)
  library(msigdbr)
  library(dplyr)
  library(glue)
  library(patchwork)
  library(ggplot2)
  library(scales)
})

int_methods = c("RPCA","Harmony")
for (m in int_methods) {
  dir.create(here("output", "figures", m), recursive = TRUE, showWarnings = FALSE)
}
dir.create(here("output", "figures", "CellType"), recursive = TRUE, showWarnings = FALSE)

# ── 2. Load central donor metadata ─────────────────────────────────────────────
donor_meta_file <- here("data", "meta_data", "donor_meta.csv")

if (!file.exists(donor_meta_file)) {
  stop(glue(
    "donor_meta.csv not found: {donor_meta_file}\n",
    "Run 00_create_donor_metadata.R first."
  ))
}

donor_meta <- read.csv(donor_meta_file, stringsAsFactors = FALSE)

# All derived from CSV — nothing hard-coded
cohort_levels  <- c("Healthy",
                    unique(donor_meta$cohort[donor_meta$cohort != "Healthy"]))
healthy_donors <- donor_meta$donor_short[donor_meta$cohort == "Healthy"]
donor_levels   <- c(donor_meta$donor_short, "doublet", "unassigned")

cat("✔ Donor metadata loaded:\n")
print(donor_meta[, c("donor_short", "cohort", "disease", "pool")])
cat(glue("Cohort levels: {paste(cohort_levels, collapse = ' > ')} \n"))

# ── Helper: attach cohort / disease / sex to Seurat object ────────────────────
attach_donor_meta <- function(seu_obj, meta = donor_meta) {
  idx             <- match(as.character(seu_obj$genetic_donor), meta$donor_short)
  seu_obj$cohort  <- factor(meta$cohort[idx], levels = cohort_levels)
  seu_obj$disease <- meta$disease[idx]
  seu_obj$sex     <- meta$sex[idx]
  seu_obj
}

# ── 3. Gene filter function ────────────────────────────────────────────────────
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
symbol_colname =  intersect(c("db_gene_symbol", "human_gene_symbol"), names(c2_sets))[1]

filter_genes <- function(sce) {
  unannot = rownames(sce)[is.na(rowData(sce)$ID)]
  if (length(unannot) > 0) {
      warning(sprintf("Removing %d genes without Gene 'ID' annotations",
              length(unannot)))
      unannot_filt = rownames(sce) %in% unannot
      sce = sce[!unannot_filt,]
              }
  mito   <- rownames(sce)[rowData(sce)$CHR == "MT"]
  ribo   <- union(
    grep("^RP(S|L)", rownames(sce), value = TRUE),
    c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ][symbol_colname]
  )
  #sex    <- rownames(sce)[rowData(sce)$CHR %in% c("X", "Y") ]
  pseudo <- rownames(sce)[grepl("pseudogene",
                                rowData(sce)$ENSEMBL.GENEBIOTYPE,
                                ignore.case = TRUE) ]
  
  uninformative <- rownames(sce) %in% c(mito, ribo, #sex,
                                        pseudo)
  junk          <- sce$genetic_donor %in% c("doublet", "unassigned") |
    sce$cell_status   == "empty_droplet"
  
  cat(glue("  Uninformative genes removed: {sum(uninformative)} \n"))
  cat(glue("  Junk cells removed:          {sum(junk)} \n\n"))
  
  sce[!uninformative, !junk]
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN — build or load annotated Seurat object
# ══════════════════════════════════════════════════════════════════════════════

pool_names = unique(donor_meta$pool)
#Add "Pool_" prefix  per name if pool names are short like "A", "B"
full_pool_names <- ifelse(
  stringr::str_length(pool_names) < 2,
  paste0("Pool_", pool_names),
  pool_names
)

#Naming convention for output file
#Only save file as "Pools_{pool1}_{pool2}" if  there are 3 or less pools
if (length(pool_names)<4){
    joint_name = paste0("Pools_",paste(pool_names, collapse="_"))
} else {
    joint_name = paste0("All_pools_", format(Sys.Date("%Y-%m-%d")))
}

out_file <- here("data", "SCEs", glue("{joint_name}.annotated.SEU.rds"))

if (file.exists(out_file)) {
  cat("✔ Annotated object exists — loading\n")
  seu <- readRDS(out_file)
  
} else {
  for (pool in full_pool_names){
    
    
    # ── Load and filter each pool ────────────────────────────────────────────
    sce_list <- lapply(full_pool_names, function(pool) {
      cat(glue("\n{strrep('─', 50)}\n  Filtering Pool {pool}\n{strrep('─', 50)}\n"))
      
      sce_file <- here("data", "SCEs", glue("{pool}.preprocessed.SCE.rds"))
      if (!file.exists(sce_file)) stop(glue("Missing: {sce_file}"))
      
      sce_this    <- readRDS(sce_file)
      sceFlt_this <- filter_genes(sce_this)
      sceFlt_this$pool <- pool
      colnames(sceFlt_this) <- paste0(colnames(sceFlt_this), "_", pool)
      
      cat(glue("  Pool {pool} cells after filter: {ncol(sceFlt_this)}\n"))
      sceFlt_this
    })
    names(sce_list) <- full_pool_names
    
    # ── Shared genes across all pools ────────────────────────────────────────
    shared_genes <- Reduce(intersect, lapply(sce_list, rownames))
    cat(glue("\n  Shared genes across {length(sce_list)} pools: {length(shared_genes)}\n"))
    sce_list <- lapply(sce_list, function(x) x[shared_genes, ])
    
    # ── Remove low-abundance genes (pooled across all pools) ─────────────────
    combined_counts <- do.call(cbind, lapply(sce_list, counts))
    keep <- rowSums(combined_counts > 0) > 20
    cat(glue("  Genes kept (>20 cells expressing): {sum(keep)}\n"))
    sce_list <- lapply(sce_list, function(x) x[keep, ])
    
    # ── Combine into a single SCE ────────────────────────────────────────────
    combined_counts  <- do.call(cbind, lapply(sce_list, counts))
    combined_colData <- do.call(rbind, lapply(sce_list, colData))
    
    sce_combined <- SingleCellExperiment(
      list(counts = combined_counts),
      colData = combined_colData,
      rowData = rowData(sce_list[[1]])   # identical across pools post-filtering
    )
    
    cat(glue("\n  Combined: {ncol(sce_combined)} cells × {nrow(sce_combined)} genes\n"))
  }
  
  # ── Convert to Seurat ────────────────────────────────────────────────────────
  cnt           <- counts(sce_combined)
  rownames(cnt) <- make.unique(rowData(sce_combined)$Symbol)
  
  seu_combined <- CreateSeuratObject(
    counts    = cnt,
    meta.data = data.frame(colData(sce_combined))
  )
  
  # ── SCTransform and integration across donors ─────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Integrating across donors\n{strrep('─', 50)}\n"))
  
  source(utility_file)
  
  seu_integrated <- intDat(seu_combined, split = "genetic_donor", type = "RNA")
  #v5 object contains latent dimensions & UMAPs for each integration method
  #integrated is no longer its own assay; instead SCT assay contains normalised counts
  
  # ── Dimensionality reduction ─────────────────────────────────────────────────
  # calculating the pre-integration UMAP (SCT)
  DefaultAssay(seu_integrated) <- "SCT" #integrated is no longer its own assay
  seu_integrated <- seu_integrated %>%
    RunUMAP(dims = 1:30, verbose = FALSE, 
            reduction = "pca", reduction.name="sct.umap")
  
  cat(glue("  PCA and UMAP complete\n"))
  
  # ── Azimuth annotation ───────────────────────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Running Azimuth annotation\n{strrep('─', 50)}\n"))
  seu <- RunAzimuth(seu_integrated, reference = "pbmcref")
  
  # ── Attach cohort / disease / sex from central CSV ────────────────────────────
  seu <- attach_donor_meta(seu)
  
  # Verify metadata attached
  if (all(is.na(seu$cohort))) {
    stop("cohort is all NA — check donor_short in donor_meta.csv matches seu$genetic_donor")
  }
  
  cat("\n  Cohort summary:\n")
  print(table(seu$cohort))
  cat("\n  Cell type summary (l2):\n")
  print(table(seu$predicted.celltype.l2))
  
  # ── Save ─────────────────────────────────────────────────────────────────────
  saveRDS(seu, out_file)
  cat(glue("  ✔ Saved: {out_file}\n"))
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════
cat(glue("\n{strrep('─', 50)}\n  Generating plots\n{strrep('─', 50)}\n"))

# Cohort colours derived from levels — no hard-coding
cohort_colours <- setNames(
  c("#4DAF4A", "#E41A1C", "#FF7F00", "#984EA3", "#377EB8")[seq_along(cohort_levels)],
  cohort_levels
)

for (m in int_methods){
  DefaultDimReduc(seu) = paste(tolower(m), "umap", sep=".")
  # ── UMAP: integration QC ──────────────────────────────────────────────────────
  pdf(here("output", "figures", m, "UMAP_integration.pdf"),
      width = 18, height = 6)
  print(
    DimPlot(seu, group.by = "pool",          pt.size = 0.3, shuffle = T)  +
      DimPlot(seu, group.by = "genetic_donor", pt.size = 0.3, shuffle=T)  +
      DimPlot(seu, group.by = "cohort",        pt.size = 0.3,
              cols = cohort_colours, shuffle=T)
  )
  dev.off()
  cat(sprintf("  ✔ %s/UMAP_integration.pdf\n",m))
  
  # ── UMAP: cell types ──────────────────────────────────────────────────────────
  pdf(here("output", "figures",m, "UMAP_celltypes_l1.pdf"),
      width = 10, height = 8)
  print(
    DimPlot(seu, group.by = "predicted.celltype.l1",
            label = TRUE, pt.size = 0.3, repel = TRUE) + NoLegend()
  )
  dev.off()
  
  pdf(here("output", "figures", m, "UMAP_celltypes_l2.pdf"),
      width = 12, height = 10)
  print(
    DimPlot(seu, group.by = "predicted.celltype.l2",
            label = TRUE, pt.size = 0.3, repel = TRUE) + NoLegend()
  )
  dev.off()
  cat(sprintf("  ✔ %s/UMAP_celltypes.pdf\n", m))
  
  # ── UMAP: per donor ───────────────────────────────────────────────────────────
  pdf(here("output", "figures", m, "UMAP_per_donor.pdf"),
      width = 18, height = 10)
  print(
    DimPlot(seu, group.by = "predicted.celltype.l2",
            split.by = "genetic_donor",
            label = TRUE, pt.size = 0.3, repel = TRUE, ncol = 3) +
      NoLegend()
  )
  dev.off()
  cat("  ✔ UMAP_per_donor.pdf\n")
}



# ── Cell type × cohort composition barplot ─────────────────────────────────────
celltype_cohort <- as.data.frame(table(
  CellType = seu$predicted.celltype.l2,
  Cohort   = seu$cohort
)) %>%
  group_by(Cohort) %>%
  mutate(pct = 100 * Freq / sum(Freq)) %>%
  ungroup() %>%
  filter(Freq > 0)

p_comp <- ggplot(celltype_cohort,
                 aes(x = Cohort, y = pct, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_x_discrete(limits = cohort_levels) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = NULL, y = "% of cells", fill = "Cell type (l2)",
       title = "Cell type composition per cohort") +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        axis.text.x     = element_text(angle = 30, hjust = 1),
        legend.position = "right",
        legend.text     = element_text(size = 7))

ggsave(here("output", "figures", "CellType", "composition_by_cohort.pdf"),
       p_comp, width = 10, height = 7)
cat("  ✔ composition_by_cohort.pdf\n")

# ── Summary tables ────────────────────────────────────────────────────────────
cat("\n  Final cell counts:\n")
cat(glue("  Total cells: {ncol(seu)}\n"))
cat("\n  Cell type (l2) × cohort:\n")
print(table(seu$predicted.celltype.l2, seu$cohort))
cat("\n  Donor × cohort:\n")
print(table(Donor = seu$genetic_donor, Cohort = seu$cohort))

cat(glue("\n{strrep('═', 50)}\n  INTEGRATION & ANNOTATION COMPLETE\n{strrep('═', 50)}\n"))
cat(glue("  Output: {out_file}\n\n"))

