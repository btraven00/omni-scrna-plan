#!/usr/bin/env Rscript
# Validate a NORM `{dataset}_normalized.h5` or a FEAT
# `{dataset}_normalized_selected.h5` (same 10x-style CSC layout, indexed by
# gene: shape = c(n_genes, n_cells), indices hold 0-based gene ids).
#
# A gene with no expression at all has variance exactly 0, so HVG selection
# ranks it against ties, PCA loads it with no signal and scran-style size
# factors can divide by zero -- all of which surface far downstream as a
# confusing number rather than an error. Fail here instead.

suppressPackageStartupMessages(library(rhdf5))

validate_file <- function(path) {
  on.exit(h5closeAll())
  ok <- TRUE

  tryCatch({
    shape <- as.integer(h5read(path, "matrix/shape"))
    n_genes <- shape[1]
    genes <- as.character(h5read(path, "matrix/genes"))

    if (length(genes) != n_genes) {
      stop(sprintf("matrix/genes has %d entries, matrix/shape says %d",
                   length(genes), n_genes))
    }

    data <- h5read(path, "matrix/data")
    indices <- as.integer(h5read(path, "matrix/indices"))
    # explicit zeros are stored by some normalizations, so don't count them
    nnz_per_gene <- tabulate(indices[data != 0] + 1L, nbins = n_genes)
    empty <- which(nnz_per_gene == 0)

    if (length(empty) > 0) {
      stop(sprintf("%d/%d genes have zero expression in every cell: %s%s",
                   length(empty), n_genes,
                   paste(genes[head(empty, 5)], collapse = ", "),
                   if (length(empty) > 5) ", ..." else ""))
    }

  }, error = function(e) {
    message(sprintf("FAIL\t%s\t%s", path, e$message))
    ok <<- FALSE
  })

  if (ok) {
    message(sprintf("OK\t%s", path))
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (identical(args, "--selftest")) {
  # ponytail: one check, not a suite. Writes a 3-gene matrix where gene g2 is
  # present only as an explicit zero, and asserts we catch it.
  p <- tempfile(fileext = ".h5")
  h5createFile(p)
  h5createGroup(p, "matrix")
  h5write(as.integer(c(3, 2)), p, "matrix/shape")
  h5write(c("g1", "g2", "g3"), p, "matrix/genes")
  h5write(c("c1", "c2"), p, "matrix/barcodes")
  h5write(c(1.5, 0.0, 2.5, 3.5), p, "matrix/data")
  h5write(as.integer(c(0, 1, 2, 0)), p, "matrix/indices")
  h5write(as.integer(c(0, 3, 4)), p, "matrix/indptr")
  out <- capture.output(validate_file(p), type = "message")
  stopifnot(grepl("^FAIL", out), grepl("1/3 genes have zero expression", out),
            grepl("g2", out))
  unlink(p)
  message("selftest OK")
} else {
  invisible(lapply(args, validate_file))
}
