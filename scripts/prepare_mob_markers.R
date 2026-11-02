# Reproduce example marker annotations from the bundled reference.
# These annotations are an installation example, not validated biological results.
suppressPackageStartupMessages(library(SingleCellExperiment))
load("spdata/MOB/Rep12_MOB_count_matrix-1.RData")
load("scdata/MOB/MOB.dge.sceset.RData")
reference <- SummarizedExperiment::assay(sce, "counts")
labels <- as.character(SummarizedExperiment::colData(sce)$cellType)
types <- c("GC", "PGC", "M/TC", "OSNs", "EPL-IN")
common <- intersect(rownames(reference), rownames(MOB_raw))
means <- sapply(types, function(label) Matrix::rowMeans(reference[common, labels == label, drop = FALSE]))
markers <- setNames(vector("list",length(types)), types)
used <- character()
for(i in seq_along(types)) {
 other <- apply(means[, -i, drop = FALSE], 1, max)
 score <- log2((means[, i] + 0.1) / (other + 0.1))
 ok <- means[,i] > 0.1 & Matrix::rowSums(MOB_raw[common, , drop = FALSE]) > 100
 candidates <- common[order(score, decreasing = TRUE)]
 candidates <- candidates[candidates %in% common[ok] & !candidates %in% used]
 markers[[i]] <- head(candidates, 4)
 used <- c(used, markers[[i]])
}
write.table(data.frame(gene=unlist(markers, use.names=FALSE), label=rep(types,lengths(markers))),
            "inst/examples/MOB_markers.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("MOB count dimensions:",dim(MOB_raw),"\n")
print(markers)
