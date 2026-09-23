#!/usr/bin/env Rscript
# Is scSHC (the CLUSTERER, not testClusters) feasible at scale, and is it
# invariant to the order cells are stored in?
#
# Two questions in one pass, because both are answered by the same runs:
#   * elapsed_s vs n_cells decides whether a 50-permutation ECDF is affordable
#     at all, or whether scSHC has to be compared on a subsample.
#   * ARI(identity, permuted) is the same stability variable order_probe.py
#     measures for the graph clusterers -- ground-truth-free, so it needs no
#     label set and transfers off pbmc.
#
# scSHC consumes COUNTS, not the PCA embedding, so the permutation has to be
# applied to the count columns here rather than to a stored graph upstream.
#
#   SCSHC_CORES=16 Rscript scshc_order_probe.R data/data.h5ad out_order "2000 5000"
suppressMessages({ library(hdf5r); library(Matrix); library(scSHC) })

# Exact ARI from the contingency table. Inline rather than via mclust: the only
# environment that has scSHC is snakemake-managed and is not ours to add to.
adjustedRandIndex <- function(a, b) {
  t <- table(a, b); n <- sum(t)
  s <- sum(choose(t, 2))
  ra <- sum(choose(rowSums(t), 2)); cb <- sum(choose(colSums(t), 2))
  e <- ra * cb / choose(n, 2)
  (s - e) / ((ra + cb) / 2 - e)
}

args   <- commandArgs(trailingOnly = TRUE)
h5     <- if (length(args) >= 1) args[[1]] else "data/data.h5ad"
outdir <- if (length(args) >= 2) args[[2]] else "out_order"
sizes  <- as.integer(strsplit(if (length(args) >= 3) args[[3]] else "2000 5000", " +")[[1]])
perms  <- as.integer(Sys.getenv("SCSHC_PERMS", "2"))
CORES  <- as.integer(Sys.getenv("SCSHC_CORES", "8"))
ALPHA  <- as.numeric(Sys.getenv("SCSHC_ALPHA", "0.05"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("loading counts from ", h5)
f  <- H5File$new(h5, "r"); g <- f[["layers/counts"]]
sh <- as.integer(h5attributes(g)$shape)
counts <- new("dgCMatrix", p = g[["indptr"]][], i = g[["indices"]][],
              x = as.numeric(g[["data"]][]), Dim = c(sh[2], sh[1]))
rownames(counts) <- f[["var/_index"]][]; colnames(counts) <- f[["obs/_index"]][]
f$close_all()
message(sprintf("counts: %d genes x %d cells", nrow(counts), ncol(counts)))

sp <- file.path(outdir, "summary.tsv")
if (!file.exists(sp))
  cat("n_cells\tperm\tk\telapsed_s\tari_vs_identity\tstatus\n", file = sp)

set.seed(42)
for (n in sizes) {
  # one fixed subsample per size; the permutation varies WITHIN it, so cell
  # composition is not confounded with storage order.
  sub  <- counts[, sample(ncol(counts), min(n, ncol(counts)))]
  base <- NULL
  for (p in 0:perms) {
    ord <- if (p == 0) seq_len(ncol(sub)) else
           sample(ncol(sub))                       # relabelling only: same cells
    t0  <- proc.time()[["elapsed"]]
    res <- tryCatch(scSHC(sub[, ord], alpha = ALPHA, cores = CORES),
                    error = function(e) e)
    el  <- proc.time()[["elapsed"]] - t0
    if (inherits(res, "error")) {
      cat(sprintf("%d\t%d\tNA\t%.1f\tNA\t%s\n", n, p, el, class(res)[1]), file = sp, append = TRUE)
      message(sprintf("[fail] n=%d perm=%d after %.1fs: %s", n, p, el, conditionMessage(res)))
      next
    }
    lab <- as.integer(factor(res[[1]]))
    back <- integer(length(lab)); back[ord] <- lab   # undo the relabelling
    if (p == 0) base <- back
    ari <- if (p == 0) 1.0 else adjustedRandIndex(base, back)
    cat(sprintf("%d\t%d\t%d\t%.1f\t%.4f\tok\n", n, p, length(unique(lab)), el, ari),
        file = sp, append = TRUE)
    message(sprintf("[ok] n=%d perm=%d k=%d %.1fs ARI=%.4f", n, p, length(unique(lab)), el, ari))
  }
}
message("done -> ", sp)
