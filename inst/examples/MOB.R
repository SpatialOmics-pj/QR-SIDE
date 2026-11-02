# Run all 282 spatial spots in the bundled MOB dataset with the tutorial's
# original model dimensions. The example marker table is derived from the
# included single-cell reference; no biological conclusions are asserted.
run_mob_example <- function(data_file = NULL, max_iter = 50L, verbose = FALSE) {
  if (is.null(data_file)) {
    data_file <- tempfile(fileext = ".RData")
    on.exit(unlink(data_file), add = TRUE)
    utils::download.file(
      "https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/spdata/MOB/Rep12_MOB_count_matrix-1.RData",
      data_file, mode = "wb", quiet = !verbose)
  }
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  data <- new.env()
  load(data_file, envir = data)
  counts <- as.matrix(Matrix::t(data$MOB_raw))
  coordinates <- do.call(rbind, strsplit(rownames(counts), "x", fixed = TRUE))
  pos <- data.frame(x = as.numeric(coordinates[, 1]), y = as.numeric(coordinates[, 2]),
                    row.names = rownames(counts))
  markers <- utils::read.delim(system.file("examples", "MOB_markers.tsv", package = "QRSIDE"))
  qc <- QRSIDE::QRSIDE_QC(counts, pos, threshold = 20)
  fit <- QRSIDE::run_QRSIDE(qc[[1]], qc[[2]], markers, Num_topic = 7,
                            Num_HVG = 1000, dim_embed = 20, max_marker = 4,
                            markerflag = FALSE, max_iter = max_iter, verbose = verbose)
  stopifnot(nrow(fit$pi) == nrow(qc[[1]]), ncol(fit$pi) == 7L,
            identical(dim(fit$beta), c(7L, 5L)),
            all(is.finite(fit$pi)), all(is.finite(fit$beta)),
            max(abs(rowSums(fit$pi) - 1)) < 1e-6,
            max(abs(rowSums(fit$beta) - 1)) < 1e-6)
  message("MOB full-data example passed: ", nrow(fit$pi), " spots, 7 topics, 5 cell types.")
  fit
}
