library(QRSIDE)
set.seed(17)
n <- 60L
p <- 80L
counts <- matrix(rpois(n * p, 0.7), n, p,
                 dimnames = list(paste0("spot", seq_len(n)), paste0("gene", seq_len(p))))
groups <- rep(1:2, each = n / 2L)
counts[, 1:6] <- counts[, 1:6] + 1
counts[groups == 1, 1:3] <- counts[groups == 1, 1:3] + 15
counts[groups == 2, 4:6] <- counts[groups == 2, 4:6] + 15
counts[groups == 1, 7:25] <- counts[groups == 1, 7:25] + matrix(rpois(30 * 19, 3), 30)
counts[groups == 2, 26:45] <- counts[groups == 2, 26:45] + matrix(rpois(30 * 20, 3), 30)
pos <- data.frame(x = rep(seq_len(10), 6), y = rep(seq_len(6), each = 10),
                  row.names = rownames(counts))
markers <- list(type_A = paste0("gene", 1:3), type_B = paste0("gene", 4:6))
table <- data.frame(gene = unlist(markers, use.names = FALSE),
                    label = rep(names(markers), lengths(markers)))
check_fit <- function(fit, topics) {
  stopifnot(identical(dim(fit$pi), c(n, as.integer(topics))),
            identical(dim(fit$beta), c(as.integer(topics), 2L)),
            all(is.finite(fit$pi)), all(is.finite(fit$beta)),
            all(fit$pi >= 0), all(fit$beta >= -1e-10),
            max(abs(rowSums(fit$pi) - 1)) < 1e-6,
            max(abs(rowSums(fit$beta) - 1)) < 1e-6,
            identical(rownames(fit$pi), rownames(counts)),
            identical(rownames(fit$pos), rownames(counts)),
            setequal(fit$order, 1:2), setequal(fit$cell_types, names(markers)))
}
for (embedding in c(1L, 2L, 5L)) {
  fit <- run_QRSIDE(counts, pos, table, 2, 50, embedding, 3, FALSE, max_iter = 8)
  check_fit(fit, 2)
  if (embedding == 2L) table_fit <- fit
  cat("PASS table input, full pipeline, embedding", embedding, "\n")
}
sparse_fit <- run_QRSIDE(Matrix::Matrix(counts, sparse = TRUE), pos, markers,
                          2, 50, 2, 3, TRUE, max_iter = 8)
stopifnot(isTRUE(all.equal(table_fit, sparse_fit)))
cat("PASS sparse matrix and equivalent marker formats produce the same fit\n")
fit <- run_QRSIDE(counts, pos, markers, 2, 50, 2, 3, TRUE)
check_fit(fit, 2)
cat("PASS list input, full pipeline, default 50 iterations\n")
reordered <- run_QRSIDE(counts, pos[rev(seq_len(n)), , drop = FALSE], markers,
                         2, 50, 2, 3, TRUE)
stopifnot(isTRUE(all.equal(fit, reordered)))
cat("PASS shuffled coordinate rows preserve the same fitted result\n")
single <- run_QRSIDE(counts, pos, markers, 1, 50, 2, 3, TRUE, max_iter = 8)
check_fit(single, 1)
cat("PASS single spatial topic\n")
expect_error <- function(expr, pattern) {
  err <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), fixed = TRUE))
}
expect_error(run_QRSIDE(counts / 2, pos, markers, 2, 50, 2, 3, TRUE), "raw counts")
expect_error(run_QRSIDE(counts, pos[-1, ], markers, 2, 50, 2, 3, TRUE), "same unique spot")
expect_error(run_QRSIDE(counts, pos, markers, 2, 50, 1000, 3, TRUE), "dim_embed")
expect_error(run_QRSIDE(counts, pos, list(A = "missing", B = "gene4"),
                       2, 50, 2, 3, TRUE), "No expressed")
expect_error(run_QRSIDE(counts, pos, markers, 2, 5000, 2, 3, TRUE), "Num_HVG")
qc_counts <- matrix(c(5, 0, 0, 1), 2, dimnames = list(c("a", "b"), c("g1", "g2")))
qc_pos <- data.frame(x = c(2, 1), y = c(0, 0), row.names = c("b", "a"))
qc <- QRSIDE_QC(qc_counts, qc_pos, 3)
stopifnot(identical(dim(qc[[1]]), c(1L, 1L)), identical(rownames(qc[[2]]), "a"), qc[[2]]$x == 1)
cat("PASS input validation and QC dimension/alignment regressions\n")
coords <- matrix(c(0, 0, 1, 0, 2, 0, 10, 10), byrow = TRUE, ncol = 2)
neighbor <- as.matrix(get(".spatial_neighbors", asNamespace("QRSIDE"))(coords))
expected <- matrix(0, 4, 4)
expected[1, 2] <- expected[2, 1] <- expected[2, 3] <- expected[3, 2] <- 1
stopifnot(identical(unname(neighbor), expected), !"SpatialDecon" %in% loadedNamespaces(),
          !"SC.MEB" %in% loadedNamespaces())
cat("PASS spatial neighbors and removal of companion-package dependencies\n")
