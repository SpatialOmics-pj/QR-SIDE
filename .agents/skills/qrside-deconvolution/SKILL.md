---
name: qrside-deconvolution
description: Run QR-SIDE spatial transcriptomics cell-type deconvolution experiments in R, including one-click installation, input validation, QC, reproducible model execution, and output handoff. Use when a request mentions QR-SIDE or QRSIDE, spatial cell-type deconvolution without a quantitative single-cell reference, qualitative marker-gene references, marker-list deconvolution, or an agent-ready QR-SIDE workflow. Do not use for generic scRNA-seq analysis or methods that require a quantitative reference matrix.
---

# QR-SIDE Deconvolution

Use QR-SIDE to turn spot-by-gene spatial counts, spot coordinates, and qualitative cell-type markers into a reproducible deconvolution result.

## Workflow

1. Inspect the available data and the user's requested outcome. Identify the count matrix, spatial coordinates, marker annotations, output location, and any existing experiment configuration or notebook before choosing parameters.
2. Check whether R >= 4.2.3 and `QRSIDE` are available. If QR-SIDE is missing and installation is in scope, use the bundled installer instead of reconstructing its dependency steps:

   ```sh
   Rscript -e 'source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R"); install_qrside()'
   ```

   Respect the environment's network and filesystem permissions. Prefer a user-writable R library; pass `lib` to `install_qrside()` when the default library is not writable.
3. Verify the installation before processing data:

   ```sh
   Rscript -e 'library(QRSIDE); stopifnot(exists("run_QRSIDE", envir = asNamespace("QRSIDE"), inherits = FALSE)); stopifnot(exists("CountDeconvolution", envir = asNamespace("SpatialDecon"), inherits = FALSE))'
   ```

4. Validate inputs before starting an expensive fit:
   - `sp_counts` is a nonnegative, spot-by-gene count matrix with unique spot and gene identifiers.
   - `pos` has the same spot identifiers, exactly two coordinate columns, finite values, and one row per retained spot.
   - Marker identifiers use the same gene namespace as `colnames(sp_counts)` and every cell type has at least one overlapping marker.
   - `Num_HVG` does not exceed the usable gene count.
   - The selected marker format matches `markerflag`.
5. Use `QRSIDE_QC()` only after aligning counts and coordinates by spot identifier. Record the threshold and report how many spots and genes remain.
6. Reuse parameters from the user's protocol, configuration, or repository when present. Do not silently invent scientific parameters. If values are absent, inspect the data and ask for `Num_topic`, `Num_HVG`, `dim_embed`, `max_marker`, and the QC threshold unless the user explicitly authorizes a documented pilot choice.
7. Write a reproducible R script before running the fit. Set a seed, save the raw `run_QRSIDE()` result as RDS, save the exact parameters and input summaries, capture `sessionInfo()`, and retain the console log. Do not overwrite existing experiment outputs without authorization.
8. Verify that the result contains `pi`, `beta`, `order`, and `pos`, that dimensions are consistent with the retained spots and requested topics, and that values are finite. Present raw model outputs separately from biological interpretation.

## Marker formats

- With `markerflag = FALSE`, `markerframe` must be a data frame containing `gene` and `label` columns. `max_marker` limits markers per cell type.
- With `markerflag = TRUE`, `markerframe` must be a list with one marker vector per cell type. QR-SIDE infers the number of cell types from `length(markerframe)`. Preserve list order and names; the function still requires `max_marker`, but does not use it in this mode.

## Scientific boundaries

- `Num_topic` is the number of latent spatial topics/domains, not the number of cell types.
- Do not substitute the Bioconductor package named `SpatialDecon`; QR-SIDE requires its project-specific companion package with `CountDeconvolution()`.
- Do not claim cell-type labels or biological conclusions that are not supported by the marker definitions and model outputs.
- Do not transform `pi` and `beta` into a new spot-by-cell-type estimate unless the requested analysis defines and documents that transformation.

For input schemas, preflight checks, an experiment template, output checks, and troubleshooting, read [references/experiment-workflow.md](references/experiment-workflow.md).
