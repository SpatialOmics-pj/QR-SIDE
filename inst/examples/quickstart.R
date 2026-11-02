# Reproducible synthetic example for checking installation; no biological inference.
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

qc <- QRSIDE_QC(counts, pos, threshold = 0)
fit <- run_QRSIDE(qc[[1]], qc[[2]], markers, Num_topic = 2, Num_HVG = 50,
                  dim_embed = 2, max_marker = 3, markerflag = TRUE)
stopifnot(all(is.finite(fit$pi)), all(is.finite(fit$beta)))
print(dim(fit$pi))
print(fit$beta)
message("QR-SIDE example completed successfully.")
