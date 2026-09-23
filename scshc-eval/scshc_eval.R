#!/usr/bin/env Rscript
# Evaluate cluster significance with scSHC over one or more cluster TSVs.
#
# scSHC::testClusters builds a pseudobulk hierarchy over the *existing* clusters
# and tests each split for significance. Non-significant splits are merged, so
# the read-off is: input K vs surviving K -> how many clusters are "real", and
# which of the input clusters collapsed together.
#
# Usage:
#   Rscript scshc_eval.R <counts.h5ad> <out_dir> <clusters.tsv> [clusters.tsv ...]
# Env:
#   SCSHC_CORES   parallel workers for the null simulation (default 4)
#
# The h5ad must carry raw integer counts in layers/counts (anndata CSR).
# Each TSV needs a cell_id column + one label column (2nd column, or named
# `cluster`/`truths`). Counts are subset+reordered to each TSV's cells.

suppressPackageStartupMessages({
  library(hdf5r); library(Matrix); library(scSHC); library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: scshc_eval.R <counts.h5ad> <out_dir> <clusters.tsv> [...]")
h5 <- args[[1]]; out_dir <- args[[2]]; tsvs <- args[-(1:2)]
CORES <- as.integer(Sys.getenv("SCSHC_CORES", "4"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load raw counts as genes x cells dgCMatrix ----
# anndata CSR-by-cell has the same p/i/x layout as a CSC genes x cells matrix.
message("loading counts from ", h5)
f <- H5File$new(h5, "r")
g <- f[["layers/counts"]]
sh <- as.integer(h5attributes(g)$shape)          # (n_cells, n_genes)
counts <- new("dgCMatrix",
  p = g[["indptr"]][], i = g[["indices"]][],
  x = as.numeric(g[["data"]][]), Dim = c(sh[2], sh[1]))
rownames(counts) <- f[["var/_index"]][]
colnames(counts) <- f[["obs/_index"]][]
f$close_all()
message(sprintf("counts: %d genes x %d cells", nrow(counts), ncol(counts)))

# short, unique id from the omnibenchmark path (fi-/nr-/fe-/pc-/nn-/cl- tokens)
make_id <- function(path) {
  toks <- regmatches(path, gregexpr("(fi|nr|fe|pc|nn|cl)-[a-z0-9]+", path))[[1]]
  if (length(toks)) paste(toks, collapse = "_") else tools::file_path_sans_ext(basename(path))
}

label_col <- function(dt) {
  for (nm in c("cluster", "truths")) if (nm %in% names(dt)) return(nm)
  setdiff(names(dt), "cell_id")[1]                # fall back to first non-id column
}

summary_path <- file.path(out_dir, "summary.tsv")
if (!file.exists(summary_path))
  cat("id\tn_cells\tinput_K\tsurviving_K\tretained_frac\telapsed_s\tstatus\n", file = summary_path)

for (tsv in tsvs) {
  id <- make_id(tsv)
  res <- tryCatch({
    cl <- fread(tsv, header = TRUE)
    lc <- label_col(cl)
    cl <- cl[cl$cell_id %in% colnames(counts)]
    data <- counts[, cl$cell_id, drop = FALSE]
    ids <- as.character(cl[[lc]])
    if (length(unique(ids)) < 2) stop("fewer than 2 clusters")

    t0 <- Sys.time()
    r <- testClusters(data, ids, cores = CORES)
    dt <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)

    new_labels <- r[[1]]
    in_k  <- length(unique(ids))
    out_k <- length(unique(new_labels))
    merges <- tapply(ids, new_labels, function(v) sort(unique(v)))
    tree_txt <- capture.output(print(r[[2]]))

    writeLines(jsonlite::toJSON(list(
      id = id, tsv = tsv, n_cells = length(ids),
      input_K = in_k, surviving_K = out_k,
      retained_frac = round(out_k / in_k, 3),
      elapsed_s = dt,
      merges = merges,            # surviving group -> original labels it absorbed
      tree = tree_txt             # printed significance tree (node names carry p-values)
    ), auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, paste0(id, ".json")))

    cat(sprintf("%s\t%d\t%d\t%d\t%.3f\t%.1f\tok\n",
                id, length(ids), in_k, out_k, out_k / in_k, dt),
        file = summary_path, append = TRUE)
    message(sprintf("[ok]  %s  K %d -> %d  (%ss)", id, in_k, out_k, dt))
    TRUE
  }, error = function(e) {
    cat(sprintf("%s\tNA\tNA\tNA\tNA\tNA\terror: %s\n", id, conditionMessage(e)),
        file = summary_path, append = TRUE)
    message(sprintf("[err] %s  %s", id, conditionMessage(e)))
    FALSE
  })
}
message("done -> ", summary_path)
