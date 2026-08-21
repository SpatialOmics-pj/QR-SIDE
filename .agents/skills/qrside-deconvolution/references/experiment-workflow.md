# QR-SIDE experiment workflow

Read this reference when preparing or running an experiment, validating files, or diagnosing an installation or model failure.

## Expected objects

### Spatial counts

`sp_counts` must have spots in rows and genes in columns. Preserve raw counts for the model input.

```r
stopifnot(is.matrix(sp_counts) || is.data.frame(sp_counts))
sp_counts <- as.matrix(sp_counts)
stopifnot(!is.null(rownames(sp_counts)), !is.null(colnames(sp_counts)))
stopifnot(!anyDuplicated(rownames(sp_counts)))
stopifnot(!anyDuplicated(colnames(sp_counts)))
stopifnot(all(is.finite(sp_counts)), all(sp_counts >= 0))
stopifnot(all(abs(sp_counts - round(sp_counts)) < sqrt(.Machine$double.eps)))
```

### Spatial coordinates

`pos` must have one row per spot and two finite coordinate columns. Align by identifiers rather than relying on current row order.

```r
stopifnot(is.matrix(pos) || is.data.frame(pos))
stopifnot(!is.null(rownames(pos)), ncol(pos) == 2L)
stopifnot(all(rownames(sp_counts) %in% rownames(pos)))
pos <- pos[rownames(sp_counts), , drop = FALSE]
stopifnot(identical(rownames(sp_counts), rownames(pos)))
stopifnot(all(is.finite(as.matrix(pos))))
```

### Marker table

Use a table when markers are stored as gene/cell-type pairs:

```r
stopifnot(is.data.frame(markerframe))
stopifnot(all(c("gene", "label") %in% names(markerframe)))
markerframe <- unique(markerframe[c("gene", "label")])
cell_types <- unique(markerframe$label)
marker_overlap <- setNames(
  lapply(cell_types, function(label) {
    intersect(markerframe$gene[markerframe$label == label], colnames(sp_counts))
  }),
  cell_types
)
stopifnot(length(marker_overlap) > 0L, all(lengths(marker_overlap) > 0L))
markerflag <- FALSE
```

### Marker list

Use a list when marker vectors are already grouped by cell type:

```r
stopifnot(is.list(markerframe), length(markerframe) > 0L)
markerframe <- lapply(markerframe, function(x) unique(intersect(x, colnames(sp_counts))))
stopifnot(all(lengths(markerframe) > 0L))
markerflag <- TRUE
```

Named list elements are strongly preferred so downstream reports retain human-readable cell-type labels.

## Reproducible run template

Replace paths and parameters with values supplied by the user or documented in the project. Do not treat the example parameter values as universal defaults.

```r
library(QRSIDE)

set.seed(1)

sp_counts <- readRDS("data/sp_counts.rds")
pos <- readRDS("data/pos.rds")
markerframe <- readRDS("data/markers.rds")

# Run the validation checks above before QC.
qc_threshold <- 20
qc <- QRSIDE_QC(sp_counts, pos, threshold = qc_threshold)
sp_counts_qc <- qc[[1]]
pos_qc <- qc[[2]]

params <- list(
  Num_topic = 5,
  Num_HVG = 2000,
  dim_embed = 15,
  max_marker = 50,
  markerflag = FALSE,
  qc_threshold = qc_threshold,
  seed = 1
)

dir.create("qrside-results", recursive = TRUE, showWarnings = FALSE)

fit <- run_QRSIDE(
  sp_counts = sp_counts_qc,
  pos = pos_qc,
  markerframe = markerframe,
  Num_topic = params$Num_topic,
  Num_HVG = params$Num_HVG,
  dim_embed = params$dim_embed,
  max_marker = params$max_marker,
  markerflag = params$markerflag
)

saveRDS(fit, "qrside-results/fit.rds")
saveRDS(params, "qrside-results/parameters.rds")
writeLines(capture.output(sessionInfo()), "qrside-results/session-info.txt")
```

When executing through an agent, redirect the full R console output to a timestamped log without hiding the live exit status.

## Result checks

```r
stopifnot(is.list(fit))
stopifnot(all(c("pi", "beta", "order", "pos") %in% names(fit)))
stopifnot(nrow(fit$pi) == nrow(sp_counts_qc))
stopifnot(nrow(fit$pos) == nrow(sp_counts_qc))
stopifnot(all(is.finite(fit$pi)), all(is.finite(fit$beta)))
```

Report:

- input and retained spot/gene counts;
- marker overlap by cell type;
- all parameter values and the random seed;
- dimensions of `pi`, `beta`, `order`, and `pos`;
- output paths and session information;
- warnings, convergence evidence available in the console log, and any limitations.

## Troubleshooting

- Missing `CountDeconvolution`: the wrong package named `SpatialDecon` is installed. Rerun the official QR-SIDE installer and recheck the namespace.
- Compilation failure: confirm R >= 4.2.3 and the platform toolchain described in the QR-SIDE README. Preserve the first compiler error, not only the final nonzero exit message.
- Empty marker error: compare marker identifiers with `colnames(sp_counts)` and report overlap per cell type. Do not silently change gene identifier systems.
- Spot mismatch: align `pos` to counts by row names and report unmatched identifiers; never align by row position alone.
- Long installation: dependency compilation and the repository archive can take several minutes. Do not terminate a progressing install solely because it is quiet.
- Model failure after QC: report retained dimensions, parameters, marker overlaps, the R error, and `sessionInfo()` before changing the experiment.
