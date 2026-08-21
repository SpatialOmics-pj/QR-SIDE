# QR-SIDE

**Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics**

QR-SIDE performs spatial cell-type deconvolution without requiring a quantitative single-cell reference. It uses qualitative marker-gene information, a mixture Poisson model, latent spot-separable topics, and a Potts model to incorporate spatial continuity.

## Installation

QR-SIDE contains compiled C++ code and currently relies on a companion package used by the original implementation. The recommended installation method is therefore the bundled installer below, which installs the required CRAN/Bioconductor dependencies, the correct companion `SpatialDecon` package, `SC.MEB`, and QR-SIDE itself.

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside()
```

Then load the package with:

```r
library(QRSIDE)
```

### Why use the installer?

QR-SIDE requires a project-specific companion package named `SpatialDecon`. This package is **not the same package as the Bioconductor package with the same name**. Installing the wrong `SpatialDecon` can lead to errors such as:

```text
could not find function "CountDeconvolution"
```

The installer explicitly downloads the QR-SIDE companion package and verifies that `CountDeconvolution()` is available before installing QR-SIDE.

### Development version

To install a non-default branch, first source the installer and then specify the branch:

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside(ref = "fix/installation")
```

### System requirements

- R >= 4.2.3
- A working C/C++ toolchain for source-package compilation
- On macOS, install the Xcode Command Line Tools if compilation tools are missing
- On Windows, install the Rtools version corresponding to your R version if needed

## Main inputs

`run_QRSIDE()` takes the following main inputs:

1. `sp_counts`: spot-by-gene spatial transcriptomics count matrix.
2. `pos`: spot-by-2 spatial coordinate matrix.
3. `markerframe`: marker-gene information. Its format depends on `markerflag`.
4. `Num_topic`: number of spatial domains/topics.
5. `Num_HVG`: number of highly variable genes used for representation learning.
6. `dim_embed`: latent embedding dimension.
7. `max_marker`: maximum number of marker genes used per cell type when `markerflag = FALSE`.
8. `markerflag`: controls the marker input format.

When `markerflag = FALSE`, `markerframe` is expected to contain `gene` and `label` columns. When `markerflag = TRUE`, `markerframe` should be a list in which each element contains the marker genes for one cell type; the number of cell types is inferred automatically from `length(markerframe)`.

## Quick start

See `tutorial/MOB.ipynb` for an example workflow.

## Troubleshooting

If installation fails during compilation, first check that a compiler toolchain is available. If a Bioconductor dependency fails to install, run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::valid()
```

If you previously installed another package named `SpatialDecon` and encounter a missing `CountDeconvolution()` error, rerun the QR-SIDE installer above; it reinstalls and validates the companion package required by QR-SIDE.

## License

GPL-3
