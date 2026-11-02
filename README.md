# QR-SIDE

[![QR-SIDE installation check](https://github.com/SpatialOmics-pj/QR-SIDE/actions/workflows/install-check.yaml/badge.svg)](https://github.com/SpatialOmics-pj/QR-SIDE/actions/workflows/install-check.yaml)

**Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics**

QR-SIDE estimates cell-type composition in spatial transcriptomics data without requiring a quantitative single-cell reference. It combines qualitative marker-gene information with non-marker genes, a mixture Poisson model, latent spot-separable topics, and a Potts model that encourages spatially continuous domains.

## Highlights

- Uses marker-gene annotations instead of a quantitative single-cell reference matrix.
- Supports marker annotations as either a two-column table or a list grouped by cell type.
- Uses both marker and non-marker genes during model fitting.
- Jointly estimates spatial topics and cell-type composition.
- Includes a one-command installer for all R, Bioconductor, and model dependencies.

## Coding-agent skill

QR-SIDE includes an Agent Skills-compatible workflow for discovering, installing, and running reproducible spatial deconvolution experiments. Codex automatically discovers the skill when working inside this repository, and compatible coding agents can import the same skill directory.

To install the skill for use in another project, ask Codex:

```text
$skill-installer Install the qrside-deconvolution skill from https://github.com/SpatialOmics-pj/QR-SIDE/tree/main/.agents/skills/qrside-deconvolution
```

Then invoke it explicitly, or describe a matching QR-SIDE deconvolution task and allow the agent to select it automatically:

```text
$qrside-deconvolution Validate my spatial counts, coordinates, and marker genes, install QR-SIDE if needed, and run a reproducible deconvolution experiment.
```

See the [QR-SIDE deconvolution skill](.agents/skills/qrside-deconvolution/SKILL.md) for its supported workflow and scientific safeguards.

## Installation

Use R >= 4.2.3 in a fresh R/RStudio session. Install the R build tools for your system first: Xcode Command Line Tools on macOS (QR-SIDE uses R's bundled Fortran runtime when available), Rtools on Windows, or R development tools on Linux.

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside()
library(QRSIDE)
```

The installer sets up CRAN/Bioconductor dependencies and downloads the small, versioned R source package. It does not download the large example-data repository. Windows/macOS dependencies use binaries when available. To choose a writable library:

```r
install_qrside(lib = "/path/to/R/library")
```

QR-SIDE 1.0.1 includes its core deconvolution and spatial-neighbor implementations. It no longer needs the separate `SpatialDecon` or `SC.MEB` packages; the unrelated Bioconductor package named `SpatialDecon` cannot interfere with the installation.

### Download the package directly

Download [QRSIDE_1.0.1.tar.gz](https://github.com/SpatialOmics-pj/QR-SIDE/raw/refs/heads/main/QRSIDE_1.0.1.tar.gz). If dependencies are already installed:

```r
install.packages("/path/to/QRSIDE_1.0.1.tar.gz", repos = NULL, type = "source")
```

To set up dependencies and then install the downloaded file:

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside(package_file = "/path/to/QRSIDE_1.0.1.tar.gz")
```

### Verify a complete run

This synthetic example runs the complete model and checks finite outputs. It is an installation check, not a biological analysis.

```r
library(QRSIDE)
source(system.file("examples", "quickstart.R", package = "QRSIDE"))
```

The real-data MOB example downloads the spatial count file (about 3 MB), uses example markers derived from the bundled reference, and runs all spatial spots with the original tutorial's 7 topics, 1,000 HVGs, and 20 embedding dimensions:

```r
source(system.file("examples", "MOB.R", package = "QRSIDE"))
fit <- run_mob_example()
saveRDS(fit, "MOB_fit.rds")
```

Example markers support reproducibility testing and require independent review before biological interpretation. Their derivation is in `scripts/prepare_mob_markers.R`.

### Linux system libraries

For Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libglpk-dev libgsl-dev
```

### Install a development branch

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside(ref = "branch-name")
```

## Input data

The main function is:

```r
run_QRSIDE(
  sp_counts,
  pos,
  markerframe,
  Num_topic,
  Num_HVG,
  dim_embed,
  max_marker,
  markerflag
)
```

| Argument | Description |
| --- | --- |
| `sp_counts` | Finite nonnegative integer raw counts with unique spot row names and gene column names. |
| `pos` | Spatial coordinates with the same spots in rows and two coordinate columns. Row names must match `sp_counts`; coordinates are aligned by those identifiers. |
| `markerframe` | Marker-gene annotations in one of the two formats shown below. |
| `Num_topic` | Number of latent spatial topics/domains. |
| `Num_HVG` | Number of highly variable genes used for representation learning. |
| `dim_embed` | Latent embedding dimension. |
| `max_marker` | Maximum number of markers used per cell type when `markerflag = FALSE`. |
| `markerflag` | Selects the marker input format: table (`FALSE`) or list (`TRUE`). |

### Marker table

Set `markerflag = FALSE` and provide a data frame containing `gene` and `label` columns:

```r
markerframe <- data.frame(
  gene = c("marker_gene_1", "marker_gene_2", "marker_gene_3"),
  label = c("cell_type_A", "cell_type_A", "cell_type_B")
)
```

### Marker list

Set `markerflag = TRUE` and provide one marker vector per cell type. QR-SIDE infers the number of cell types from the length of the list. At least two cell types are required.

```r
markerframe <- list(
  cell_type_A = c("marker_gene_1", "marker_gene_2"),
  cell_type_B = c("marker_gene_3", "marker_gene_4")
)
```

Marker genes must use the same identifiers as `colnames(sp_counts)`. Each cell type needs at least one expressed marker not already assigned to an earlier cell type. Marker genes remain in the model even when not selected as HVGs. Choose `dim_embed` no larger than the number of variable non-marker genes or the number of spots minus one. Each retained spot needs positive marker and non-marker counts.

## Example workflow

```r
library(QRSIDE)

# sp_counts: spot-by-gene count matrix
# pos: spot-by-two coordinate matrix or data frame

qc <- QRSIDE_QC(sp_counts, pos, threshold = 20)
sp_counts_qc <- qc[[1]]
pos_qc <- qc[[2]]

fit <- run_QRSIDE(
  sp_counts = sp_counts_qc,
  pos = pos_qc,
  markerframe = markerframe,
  Num_topic = 5,
  Num_HVG = 2000,
  dim_embed = 15,
  max_marker = 50,
  markerflag = FALSE
)
```

`run_QRSIDE()` returns a list containing:

| Element | Description |
| --- | --- |
| `pi` | Estimated topic weights for the spatial spots. |
| `beta` | Estimated cell-type composition associated with the topics. |
| `order` | Cell-type column matching used by the fitted model. |
| `pos` | Spatial coordinates used in the analysis, in the input spot order. |
| `cell_types` | Cell-type names in the fitted `beta` column order. |

Optional arguments `max_iter = 50L`, `seed = 999L`, and `verbose = FALSE` control iteration limits, reproducible initialization, and diagnostics. The existing eight positional arguments remain supported.

For a runnable MOB example, use `run_mob_example()` above. [tutorial/MOB.ipynb](tutorial/MOB.ipynb) is the extended Seurat/CARD workflow and requires those optional packages.

## Troubleshooting

### A package fails while compiling

Confirm that the build tools listed under [Installation](#installation) are installed, restart R, and rerun `install_qrside()`.

### A Bioconductor dependency is unavailable

Check whether the R and Bioconductor installation is consistent:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::valid()
```

### `CountDeconvolution()` cannot be found

This indicates an older QR-SIDE installation. Restart R and rerun the installer to update to 1.0.1 or later. The current package contains this implementation internally.

### Installation still fails

Open a [GitHub issue](https://github.com/SpatialOmics-pj/QR-SIDE/issues) and include the R version, operating system, complete installation output, and the result of `sessionInfo()`.

## License

QR-SIDE is distributed under the [GPL-3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
