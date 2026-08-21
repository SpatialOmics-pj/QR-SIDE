# QR-SIDE

[![QR-SIDE installation check](https://github.com/SpatialOmics-pj/QR-SIDE/actions/workflows/install-check.yaml/badge.svg)](https://github.com/SpatialOmics-pj/QR-SIDE/actions/workflows/install-check.yaml)

**Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics**

QR-SIDE estimates cell-type composition in spatial transcriptomics data without requiring a quantitative single-cell reference. It combines qualitative marker-gene information with non-marker genes, a mixture Poisson model, latent spot-separable topics, and a Potts model that encourages spatially continuous domains.

## Highlights

- Uses marker-gene annotations instead of a quantitative single-cell reference matrix.
- Supports marker annotations as either a two-column table or a list grouped by cell type.
- Uses both marker and non-marker genes during model fitting.
- Jointly estimates spatial topics and cell-type composition.
- Includes a one-command installer for all R, Bioconductor, and companion-package dependencies.

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

QR-SIDE requires R 4.2.3 or later and compiles C++ source code. The recommended installer sets up the required CRAN and Bioconductor packages, installs compatible versions of `STdeconvolve` and `SC.MEB`, installs the QR-SIDE companion `SpatialDecon` package, and then validates the QR-SIDE installation.

Run the following commands in R or RStudio:

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside()
```

Load the installed package:

```r
library(QRSIDE)
```

To install into a specific R library and allow Bioconductor dependency updates, use:

```r
install_qrside(lib = "/path/to/R/library", upgrade = TRUE)
```

### Why the bundled installer is required

QR-SIDE uses a project-specific companion package named `SpatialDecon`. It is **not** the Bioconductor package with the same name. Installing the other package causes errors such as:

```text
could not find function "CountDeconvolution"
```

The bundled installer downloads the correct companion package and verifies that both `SpatialDecon::CountDeconvolution()` and `QRSIDE::run_QRSIDE()` are available.

### Install another branch

Once the installer is available on `main`, a development branch can be installed by passing its name:

```r
source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
install_qrside(ref = "branch-name")
```

## System requirements

- R >= 4.2.3
- A C++14-capable build toolchain
- macOS: Xcode Command Line Tools and the compilers required by the installed R version
- Windows: the Rtools release corresponding to the installed R version
- Linux: standard R development tools and system libraries required by compiled CRAN/Bioconductor packages

For Ubuntu, the system libraries used by the clean-install test are:

```bash
sudo apt-get update
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libglpk-dev
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
| `sp_counts` | Spatial count matrix with spots in rows and genes in columns. |
| `pos` | Spatial coordinates with the same spots in rows and two coordinate columns. Row names must match `sp_counts`. |
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

Set `markerflag = TRUE` and provide one marker vector per cell type. QR-SIDE infers the number of cell types from the length of the list.

```r
markerframe <- list(
  cell_type_A = c("marker_gene_1", "marker_gene_2"),
  cell_type_B = c("marker_gene_3", "marker_gene_4")
)
```

Marker genes must use the same identifiers as the column names of `sp_counts`.

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
| `pos` | Spatial coordinates used in the analysis. |

For a complete analysis example, see [tutorial/MOB.ipynb](tutorial/MOB.ipynb).

## Troubleshooting

### A package fails while compiling

Confirm that the build tools listed under [System requirements](#system-requirements) are installed, restart R, and rerun `install_qrside()`.

### A Bioconductor dependency is unavailable

Check whether the R and Bioconductor installation is consistent:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::valid()
```

### `CountDeconvolution()` cannot be found

Another package named `SpatialDecon` is installed. Rerun the QR-SIDE installer; it replaces and validates the companion package required by QR-SIDE.

### Installation still fails

Open a [GitHub issue](https://github.com/SpatialOmics-pj/QR-SIDE/issues) and include the R version, operating system, complete installation output, and the result of `sessionInfo()`.

## License

QR-SIDE is distributed under the [GPL-3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
