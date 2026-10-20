# ══════════════════════════════════════════════════════════════════════════════
# 01_load_cellranger.R — Load count matrices (onek1k) → SCE
# ══════════════════════════════════════════════════════════════════════════════

setwd("/g/data/tb21/sp0427/projects/onek1k")

library(here)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(glue)
library(data.table)

here()
pool_names = c("GSM5899873")

# Verify paths exist before the loop
for (i in pool_names) {
  p <- here("raw_scRNAseq", i, paste0(i, "_OneK1K_Pool1_RawCounts.csv.gz"))
  cat(i, ":", ifelse(file.exists(p), "✔ EXISTS", "✘ NOT FOUND"), "\n")
  cat("  ", p, "\n")
  
  b <- here("raw_scRNAseq", i, paste0(i, "_OneK1K_scRNA_Sample1_Individual_Barcodes.csv.gz"))
  cat(i, ":", ifelse(file.exists(b), "✔ EXISTS", "✘ NOT FOUND"), "\n")
  cat("  ", b, "\n")
}

dir.create(here("data", "SCEs"),     recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "barcodes"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "QC"), recursive = TRUE, showWarnings = FALSE)

for (i in pool_names) {

  out_file <- here("data", "SCEs", glue("{i}.CellRanger.SCE.rds"))

  if (file.exists(out_file)) {
    cat("✔ Pool", i, "already loaded — skipping\n")
    next
  }

  cat("\n── Loading Pool", i, "──────────────────────────\n")


  # Path to counts matrix
  counts_path <- here("raw_scRNAseq", i, paste0(i, "_OneK1K_Pool1_RawCounts.csv.gz"))

  counts = fread(counts_path, header = TRUE)
  stopifnot(!anyDuplicated(rownames(counts)))
  gene_ids <- counts[[1]]
  counts_mat <- as.matrix(counts[, -1, with = FALSE])
  rownames(counts_mat) <- gene_ids
  counts_mat = Matrix::Matrix(counts_mat, sparse = TRUE) 
  
  
  barcode_path =  here("raw_scRNAseq", i, paste0(i, "_OneK1K_scRNA_Sample1_Individual_Barcodes.csv.gz"))
  meta <- fread(barcode_path, header = TRUE)
  meta <- meta[match(colnames(counts_mat), meta$Barcode), ]
  head(meta)
  stopifnot(all(colnames(counts_mat) == meta$Barcode))
  stopifnot(sum(is.na(meta$`Individual ID`)) ==0)
  
  sce <- SingleCellExperiment(assays = list(counts = counts_mat),
                              colData = DataFrame(genetic_donor = meta$`Individual ID`,
                                                  barcode       = meta$Barcode))
  
  
  stopifnot(!anyDuplicated(colnames(sce)))

  # Split alternative experiments if present
  if ("Type" %in% colnames(rowData(sce))) {
    sce <- splitAltExps(sce, rowData(sce)$Type, "Gene Expression")
  }

  cat("  Cells loaded:", ncol(sce), "\n")
  cat("  Genes loaded:", nrow(sce), "\n")

  saveRDS(sce, out_file, compress = "xz")
  cat("  ✔ Saved:", out_file, "\n")

  # Save barcodes
  writeLines(
    colnames(sce),
    here("data", "barcodes", glue("Pool_{i}_barcodes.txt"))
  )
}
