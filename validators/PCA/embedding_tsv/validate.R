#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

validate_file <- function(path) {
  ok <- TRUE

  tryCatch({
    df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)

    if (!identical(colnames(df)[1:2], c("cell_id", "PC1"))) {
      stop(sprintf("columns must start cell_id, PC1; got %s", paste(head(colnames(df), 2), collapse = ", ")))
    }
    if (nrow(df) == 0) stop("empty file")
    if (ncol(df) < 10) stop("too few columns")
    if (anyNA(df[1, -1])) stop("first data row contains non-numeric values in PC columns")

  }, error = function(e) {
    message(sprintf("FAIL\t%s\t%s", path, e$message))
    ok <<- FALSE
  })

  if (ok) {
    message(sprintf("OK\t%s", path))
  }
  ok
}

ok <- vapply(args, validate_file, logical(1))
if (!all(ok)) quit(status = 1)
