# ══════════════════════════════════════════════════════════════════════════════
# 05_DE_analysis_celltype.R
# Cell-type-aware Sanity + limma DE analysis
# Scripts are from /g/data/pq84/single_cell/10X/IBM/blood.R
#
# USAGE: Set CELLTYPE at the top, then source the script.
#   Rscript 05_DE_analysis_celltype.R [CELLTYPE] [--wd /path/to/project_root]
#   CELLTYPE defaults to "B"; --wd defaults to the current directory.
#   --wd test  → uses the test_data project
#
# Supported CELLTYPE values:
#   "all_cells"  — no cell-type subsetting
#   "CD4"        — CD4 TCM, TEM, Naive, CTL
#   "CD8"        — CD8 TCM, TEM, Naive, CTL
#   "Monocytes"  — CD14 Mono, CD16 Mono, cDC1, cDC2
#   "Naive"      — CD4 Naive + CD8 Naive (combined naive T)
#   "NK"         — NK, NK_CD56bright
#   "B"          — B naive, B intermediate, B memory, Plasmablast
# ══════════════════════════════════════════════════════════════════════════════

args <- commandArgs(trailingOnly = TRUE)
args <- args[args != "--args"]  # tolerate a stray --args

WD_FLAGS <- c("--wd", "--project")
flag_idx <- which(args %in% WD_FLAGS)

if (length(flag_idx) > 0) {
  val_idx <- flag_idx[1] + 1
  if (val_idx > length(args) || startsWith(args[val_idx], "--")) {
    stop("--wd needs a directory path after it")
  }
  PROJECT_ROOT <- args[val_idx]
} else {
  PROJECT_ROOT <- getwd()
}

if (PROJECT_ROOT == "test") {
  PROJECT_ROOT <- "/g/data/pq84/sarahhp/projects/scauto/test_data"
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = TRUE)

# Positional args = anything that isn't a flag or a flag's value
positional <- args[setdiff(seq_along(args), c(flag_idx, flag_idx + 1))]
CELLTYPE   <- if (length(positional) >= 1) positional[1] else "B"

VALID_CELLTYPES <- c("all_cells", "CD4", "CD8", "Monocytes", "Naive", "NK", "B")
if (!CELLTYPE %in% VALID_CELLTYPES) {
  stop("Unknown CELLTYPE '", CELLTYPE, "'. Choose one of: ",
       paste(VALID_CELLTYPES, collapse = ", "))
}

# ── Project root check ─────────────────────────────────────────────────────────
has_marker <- file.exists(file.path(PROJECT_ROOT, ".here")) ||
              any(grepl("\\.Rproj$", list.files(PROJECT_ROOT)))
if (!has_marker) {
  stop("project_root (", PROJECT_ROOT, ") doesn't look right — ",
       "expected a .Rproj or .here file in that directory.")
}

setwd(PROJECT_ROOT)
try(detach("package:here", unload = TRUE), silent = TRUE)
library(here)

if (normalizePath(here()) != PROJECT_ROOT) {
  stop(glue::glue("here() resolved to {here()}, expected {PROJECT_ROOT}"))
}
cat("✔ Project root:", here(), "\n")
cat("✔ Cell type:   ", CELLTYPE, "\n")

# ── 0. Parameters ──────────────────────────────────────────────────────────────
SANITY_EXE <- "/g/data/pq84/andrew/software/Sanity/bin/Sanity"

N_CORES    <- 4
FDR_THRESH <- 0.05

# ── 1. Libraries ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(statmod)
  library(doParallel)
  library(foreach)
  library(glue)
})

sep <- strrep("═", 60)
getwd()
# ══════════════════════════════════════════════════════════════════════════════
# LOAD CENTRAL DONOR METADATA
# ══════════════════════════════════════════════════════════════════════════════
donor_meta_file <- here("data", "meta_data", "donor_meta.csv")

if (!file.exists(donor_meta_file)) {
  stop(glue(
    "donor_meta.csv not found at: {donor_meta_file}\n",
    "Please run 00_create_donor_metadata.R first."
  ))
}

donor_meta <- read.csv(donor_meta_file, stringsAsFactors = FALSE)

cohort_levels  <- c("Healthy",
                    unique(donor_meta$cohort[donor_meta$cohort != "Healthy"]))
healthy_donors <- donor_meta$donor_short[donor_meta$cohort == "Healthy"]

cat("Donor metadata loaded:\n")
print(donor_meta[, c("donor_short", "cohort", "disease", "pool")])
cat(glue("\nCohort levels: {paste(cohort_levels, collapse=' > ')}\n"))

# Helper: attach cohort / disease / sex to any Seurat object
attach_donor_meta <- function(seu_obj, meta = donor_meta) {
  idx             <- match(as.character(seu_obj$genetic_donor), meta$donor_short)
  seu_obj$cohort  <- factor(meta$cohort[idx], levels = cohort_levels)
  seu_obj$disease <- meta$disease[idx]
  seu_obj$sex     <- meta$sex[idx]
  seu_obj
}

cat(glue("\n{sep}\n  Cell type : {CELLTYPE}\n{sep}\n\n"))

# ── 2. Cell-type definitions ───────────────────────────────────────────────────
CELLTYPE_MAP <- list(
  all_cells  = NULL,
  CD4        = c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"),
  CD8        = c("CD8 TCM", "CD8 TEM", "CD8 Naive", "CD8 CTL"),
  Monocytes  = c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2"),
  Naive      = c("CD4 Naive", "CD8 Naive"),
  NK         = c("NK", "NK_CD56bright"),
  B          = c("B naive", "B intermediate", "B memory", "Plasmablast")
)

if (!CELLTYPE %in% names(CELLTYPE_MAP)) {
  stop(glue(
    "Unknown CELLTYPE '{CELLTYPE}'.\n",
    "Choose from: {paste(names(CELLTYPE_MAP), collapse=', ')}"
  ))
}

target_subtypes <- CELLTYPE_MAP[[CELLTYPE]]

# ── 3. Load Seurat object ──────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
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

SEU_PATH <- here("data", "SCEs", glue("{joint_name}.annotated.SEU.rds"))
seu <- readRDS(SEU_PATH)
#VariableFeatures(seu) <- rownames(seu@assays$integrated)

# ── 4. Subset by cell type ─────────────────────────────────────────────────────
if (!is.null(target_subtypes)) {
  seu_sub <- subset(seu, predicted.celltype.l2 %in% target_subtypes)
  cat(glue("Subsetting to {CELLTYPE}: {ncol(seu_sub)} cells retained\n"))
} else {
  seu_sub <- seu
  cat(glue("Using all cells: {ncol(seu_sub)} cells\n"))
}

# ── 5. Attach donor metadata ───────────────────────────────────────────────────
# Replaces hard-coded case_when() — adds cohort, disease, sex
seu_sub <- attach_donor_meta(seu_sub)

cat("\nCohort breakdown:\n")
print(table(seu_sub$cohort))
cat("\nDonor breakdown:\n")
print(table(seu_sub$genetic_donor))

# Verify metadata attached correctly
if (all(is.na(seu_sub$cohort))) {
  stop("cohort is all NA — attach_donor_meta() failed. Check donor_meta.csv donor_short names match seu$genetic_donor.")
}

# ── 6. Read count QC summary ───────────────────────────────────────────────────
DefaultAssay(seu_sub) <- "RNA"
seu_sub <- JoinLayers(seu_sub)

count_matrix_full <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")

per_donor <- seu_sub@meta.data %>%
  mutate(
    cell_id      = rownames(.),
    total_counts = Matrix::colSums(count_matrix_full),
    n_genes      = Matrix::colSums(count_matrix_full > 0)
  ) %>%
  group_by(genetic_donor, cohort, disease) %>%
  summarise(
    n_cells           = n(),
    total_reads       = sum(total_counts),
    mean_reads_cell   = round(mean(total_counts),   1),
    median_reads_cell = round(median(total_counts), 1),
    mean_genes_cell   = round(mean(n_genes),        1),
    .groups = "drop"
  ) %>%
  mutate(cohort = factor(cohort, levels = cohort_levels)) %>%
  arrange(cohort)

cat("\nRead count summary:\n")
print(per_donor)

# ── 7. Output directories ──────────────────────────────────────────────────────
#Making `CELLTYPE` the root dir for DE analysis
out_dir    <- here("output", "DE_analysis", CELLTYPE)
sanity_dir <- here(out_dir, "sanity")
dir.create(sanity_dir, recursive = TRUE, showWarnings = FALSE)

#Adding separate dir for limma output
limma_dir <- here(out_dir, "limma")
dir.create(limma_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(per_donor,
          here("output", "DE_analysis", glue("read_counts_per_donor_{CELLTYPE}.csv")),
          row.names = FALSE)

# ── 8. Build Sanity count matrix ───────────────────────────────────────────────
#VariableFeatures(seu_sub) <- rownames(seu_sub@assays$integrated)
filtered_features          <- rownames(seu_sub@assays$SCT)

countmatrix <- as.matrix(
  GetAssayData(seu_sub, assay = "RNA", layer = "counts")[filtered_features, ]
)
cat(glue("\nCount matrix: {nrow(countmatrix)} genes × {ncol(countmatrix)} cells\n"))

forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity  <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

countmatrix_file <- here(sanity_dir, "countmatrix.txt")
if (!file.exists(countmatrix_file)) {
  write.table(forsanity, file = countmatrix_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("✔ Count matrix written\n")
} else {
  cat("✔ Count matrix already exists, skipping write\n")
}

# ── 9. Run Sanity ──────────────────────────────────────────────────────────────
likelihood_file <- here(sanity_dir, "likelihood.txt")
if (!file.exists(likelihood_file)) {
  if (!file.exists(SANITY_EXE)) stop("Sanity not found at: ", SANITY_EXE)
  cat("Running Sanity normalisation...\n")
  system(paste(SANITY_EXE, "-n 28 -e 1 -f", countmatrix_file, "-d", sanity_dir))
} else {
  cat("✔ Sanity output already exists, skipping\n")
}

# ── 10. Load Sanity LTQs ───────────────────────────────────────────────────────
ltq_file <- here(sanity_dir, "log_transcription_quotients.txt")

# Diagnostic info
cat(glue("File exists: {file.exists(ltq_file)}\n"))
cat(glue("File size  : {file.size(ltq_file)} bytes\n"))
cat("First line (120 chars):\n")
cat(substr(readLines(ltq_file, n = 1), 1, 120), "\n\n")

# Read — check.names=FALSE preserves "-" in barcode names
expr_raw <- read.table(
  ltq_file,
  sep              = "\t",
  header           = TRUE,
  row.names        = 1,
  stringsAsFactors = FALSE,
  check.names      = FALSE   # critical: stops "-" → "." substitution in names
)

cat(glue("Raw dim: {nrow(expr_raw)} rows × {ncol(expr_raw)} cols\n"))
cat("First rowname:", rownames(expr_raw)[1], "\n")
cat("First colname:", colnames(expr_raw)[1], "\n")

# Sanity output should be genes × cells
# Detect if transposed: barcodes contain "-" digits; gene names don't
looks_like_barcode <- function(x) grepl("[ACGT]{8,}|-[0-9]+|_[AB]$", x)

if (looks_like_barcode(rownames(expr_raw)[1]) &&
    !looks_like_barcode(colnames(expr_raw)[1])) {
  cat("⚠ Cells detected as rows — transposing to genes × cells\n")
  expr_raw <- t(expr_raw)
} else {
  cat("✔ Orientation correct: rows=genes, cols=cells\n")
}

expr <- data.matrix(expr_raw)
expr <- expr - min(expr)
cat(glue("Expression matrix: {nrow(expr)} genes × {ncol(expr)} cells\n"))
cat("Gene example:", rownames(expr)[1], "\n")
cat("Cell example:", colnames(expr)[1], "\n")

# ── 11. Alignment check — Seurat ↔ Sanity ─────────────────────────────────────
seurat_cells <- colnames(seu_sub)
sanity_cells <- colnames(expr)

cat(glue("\nSeurat: {length(seurat_cells)} cells — e.g. {seurat_cells[1]}\n"))
cat(glue("Sanity: {length(sanity_cells)} cells — e.g. {sanity_cells[1]}\n"))
cat(glue("Overlap: {length(intersect(seurat_cells, sanity_cells))}\n"))

if (identical(seurat_cells, sanity_cells)) {
  # Perfect match — nothing to do
  cat("✔ Cell order matches exactly\n")
  
} else if (length(intersect(seurat_cells, sanity_cells)) > 0) {
  # Partial overlap — align to common cells in Sanity order
  common  <- intersect(sanity_cells, seurat_cells)
  expr    <- expr[, common, drop = FALSE]
  # Use integer index — avoids Seurat subset() name-matching bug
  seu_sub <- seu_sub[, match(common, colnames(seu_sub))]
  cat(glue("✔ Realigned to {length(common)} common cells\n"))
  
} else {
  stop(glue(
    "No overlap between Seurat and Sanity cell names.\n",
    "Seurat example: {seurat_cells[1]}\n",
    "Sanity example: {sanity_cells[1]}\n",
    "Likely cause: countmatrix written before JoinLayers() changed cell names,\n",
    "or check.names=FALSE not applied. Delete countmatrix.txt and re-run."
  ))
}

# Hard stop if dimensions still don't match
stopifnot(
  "expr and seu_sub cell counts differ after alignment" =
    ncol(expr) == ncol(seu_sub)
)
cat(glue("✔ Final: {ncol(expr)} cells × {nrow(expr)} genes\n"))

# ── 12. Metadata vectors ───────────────────────────────────────────────────────
# Extract AFTER alignment — guaranteed same order as expr columns
cohort <- as.character(seu_sub$cohort)
donor  <- as.character(seu_sub$genetic_donor)

cohorts_present <- sort(unique(cohort))
cat("Cohorts present:", paste(cohorts_present, collapse = ", "), "\n")

# ── 13. Design matrix ──────────────────────────────────────────────────────────
design <- model.matrix(~0 + cohort)
colnames(design) <- gsub("cohort", "", colnames(design))
cat("Design columns:", paste(colnames(design), collapse = ", "), "\n")

# ── 14. duplicateCorrelation (parallel) ────────────────────────────────────────
corfit_file <- here(limma_dir, glue("corfit_{CELLTYPE}.RData"))

if (file.exists(corfit_file)) {
  cat("✔ Loading cached corfit\n")
  load(corfit_file)
} else {
  cat(glue("Calculating duplicateCorrelation on {N_CORES} cores...\n"))
  
  ngenes       <- nrow(expr)
  nbeta        <- ncol(design)
  MaxBlockSize <- max(table(donor))
  weights      <- getEAWP(expr)$weights
  Array        <- donor
  nafun        <- function(e) NA
  
  engine <- function(i) {
    y <- drop(expr[i, ])
    o <- is.finite(y)
    A <- factor(Array[o])
    nobs    <- sum(o)
    nblocks <- length(levels(A))
    if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < (nobs - 1L)) {
      y <- y[o]
      X <- design[o, , drop = FALSE]
      Z <- model.matrix(~0 + A)
      s <- if (!is.null(weights)) {
        w <- drop(weights[i, ])[o]
        tryCatch(suppressWarnings(
          mixedModel2Fit(y, X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp),
          error = nafun)
      } else {
        tryCatch(suppressWarnings(
          mixedModel2Fit(y, X, Z, only.varcomp = TRUE, maxit = 20)$varcomp),
          error = nafun)
      }
      if (!is.na(s[1])) s[2] / sum(s)
    }
  }
  
  registerDoParallel(N_CORES)
  rho <- foreach(i = seq_len(ngenes), .combine = "c") %dopar% { engine(i) }
  stopImplicitCluster()
  
  rhomax <- 0.99
  rhomin <- 1 / (1 - MaxBlockSize) + 0.01
  m <- min(rho, 0, na.rm = TRUE); if (m < rhomin) rho[rho < rhomin] <- rhomin
  m <- max(rho, 0, na.rm = TRUE); if (m > rhomax) rho[rho > rhomax] <- rhomax
  
  arho   <- atanh(rho)
  mrho   <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
  corfit <- list(consensus.correlation = mrho, cor = mrho,
                 atanh.correlations    = arho)
  
  cat(glue("Consensus correlation (rho): {round(mrho, 4)}\n"))
  save(corfit, file = corfit_file)
}

# ── 15. Build contrasts ────────────────────────────────────────────────────────
if (!"Healthy" %in% cohorts_present) {
  stop("Healthy cohort not found in this cell-type subset — check cell numbers per donor.")
}

disease_present <- intersect(
  donor_meta$cohort[donor_meta$cohort != "Healthy"],
  cohorts_present
)

contrast_strings <- setNames(
  paste0(disease_present, " - Healthy"),
  paste0(disease_present, "_vs_Healthy")
)

cat("Contrasts to fit:\n")
for (nm in names(contrast_strings)) cat("  ", nm, "\n")

cm <- makeContrasts(contrasts = contrast_strings, levels = design)
colnames(cm) <- names(contrast_strings)
library(ggrepel)
# ── 16. Fit limma model ────────────────────────────────────────────────────────
cat("Fitting limma model...\n")
fit <- lmFit(expr, design,
             block       = donor,
             correlation = corfit$consensus.correlation)
save(fit, file = here(limma_dir, glue("limmafit_{CELLTYPE}.RData")))

fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# ── 17. Extract results, GSEA rank files & volcano plots ──────────────────────
results_list <- list()

for (contrast_name in colnames(cm)) {
  cat(glue("\nExtracting: {contrast_name}\n"))
  
  tt    <- topTable(fit, coef = contrast_name, number = Inf, adjust.method = "BH")
  n_sig <- sum(tt$adj.P.Val < FDR_THRESH)
  cat(glue("  Significant (FDR < {FDR_THRESH}): {n_sig}\n"))
  
  tt$gene <- sub(".*_", "", rownames(tt))
  #Saving DE tables to DE_analysis/CELLTYPE
  write.csv(tt, here(out_dir, glue("{contrast_name}_topTable.csv")))
  
  # GSEA rank file
  gseatt <- data.frame(gene = tt$gene, t = tt$t) %>%
    filter(!duplicated(gene)) %>%
    arrange(desc(t))
  write.table(gseatt,
              file  = here(limma_dir, glue("{contrast_name}_gsea.rnk")),
              sep   = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Volcano plot
  disease_label <- gsub("_vs_Healthy", "", contrast_name)
  volcano_df <- data.frame(
    gene  = tt$gene,
    logFC = tt$logFC,
    logP  = -log10(tt$P.Value),
    sig   = tt$adj.P.Val < FDR_THRESH
  )
  top_genes <- head(volcano_df[order(tt$adj.P.Val), ], 20)
  
  p <- ggplot(volcano_df, aes(x = logFC, y = logP)) +
    geom_point(aes(color = sig), size = 1.5, alpha = 0.7) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    box.padding = 0.5, size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("grey70", "red"), guide = "none") +
    geom_vline(xintercept = 0,            linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
    annotate("text", x = -Inf, y = Inf,
             label = "← Higher in Healthy",
             hjust = -0.05, vjust = 1.5,
             color = "steelblue", fontface = "bold", size = 3.5) +
    annotate("text", x =  Inf, y = Inf,
             label = paste0("Higher in ", disease_label, " →"),
             hjust = 1.05, vjust = 1.5,
             color = "firebrick", fontface = "bold", size = 3.5) +
    labs(
      x        = "Log2 Fold Change",
      y        = "-log10(P-value)",
      title    = glue("{CELLTYPE}: {contrast_name}"),
      subtitle = glue("{n_sig} genes at FDR < {FDR_THRESH}")
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
  
  ggsave(here(limma_dir, glue("{contrast_name}_volcano.pdf")),
         plot = p, width = 10, height = 8)
  
  results_list[[contrast_name]] <- tt
}

# ── 18. Final summary ──────────────────────────────────────────────────────────
cat(glue("\n{sep}\n  DE ANALYSIS COMPLETE — {CELLTYPE}\n{sep}\n"))
cat(glue("  Cells analysed : {ncol(expr)}\n"))
cat(glue("  Genes tested   : {nrow(expr)}\n\n"))

for (contrast_name in names(results_list)) {
  n_sig <- sum(results_list[[contrast_name]]$adj.P.Val < FDR_THRESH)
  cat(sprintf("  %-35s %d significant genes\n", contrast_name, n_sig))
}

saveRDS(results_list, here(limma_dir, glue("all_DE_results_{CELLTYPE}.rds")))
cat(glue("\n✔ All results saved to: {limma_dir}\n\n"))

