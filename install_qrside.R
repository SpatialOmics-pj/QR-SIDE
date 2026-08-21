# One-click installer for QR-SIDE
#
# Run from R/RStudio with:
# source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
# install_qrside()

install_qrside <- function(
    repo = "SpatialOmics-pj/QR-SIDE",
    ref = "main",
    lib = .libPaths()[1],
    upgrade = FALSE,
    verbose = TRUE
) {
  if (getRversion() < "4.2.3") {
    stop("QR-SIDE requires R >= 4.2.3. Current R version: ", getRversion())
  }

  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib <- normalizePath(lib, mustWork = TRUE)
  .libPaths(unique(c(lib, .libPaths())))

  repos <- getOption("repos")
  cran_repo <- unname(repos["CRAN"])
  if (!length(cran_repo) || is.na(cran_repo) || identical(cran_repo, "@CRAN@")) {
    repos["CRAN"] <- "https://cloud.r-project.org"
    options(repos = repos)
  }

  msg <- function(...) {
    if (isTRUE(verbose)) message(...)
  }

  cran_pkgs <- c(
    "BH", "Rcpp", "RcppArmadillo", "clue", "combinat", "mclust", "purrr"
  )

  bioc_pkgs <- c(
    "BiocSingular", "SingleCellExperiment", "STdeconvolve", "scater", "scran"
  )

  missing_cran <- cran_pkgs[
    !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_cran)) {
    msg("Installing CRAN dependencies: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, lib = lib)
  }

  # QR-SIDE 1.0.0 requires MASS >= 7.3-60. Older R installations may ship
  # an earlier recommended version and therefore need an explicit update.
  if (!requireNamespace("MASS", quietly = TRUE) ||
      utils::packageVersion("MASS") < numeric_version("7.3-60")) {
    msg("Installing MASS >= 7.3-60...")
    try(install.packages("MASS", lib = lib), silent = TRUE)

    if (!requireNamespace("MASS", quietly = TRUE) ||
        utils::packageVersion("MASS") < numeric_version("7.3-60")) {
      mass_archive <- paste0(
        "https://cran.r-project.org/src/contrib/Archive/MASS/",
        "MASS_7.3-60.tar.gz"
      )
      install.packages(
        mass_archive,
        repos = NULL,
        type = "source",
        lib = lib
      )
    }
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    msg("Installing BiocManager...")
    install.packages("BiocManager", lib = lib)
  }

  missing_bioc <- bioc_pkgs[
    !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_bioc)) {
    msg("Installing Bioconductor dependencies: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(
      missing_bioc,
      lib = lib,
      ask = FALSE,
      update = isTRUE(upgrade)
    )
  }

  # QR-SIDE currently uses getneighborhood_fast() from SC.MEB.
  if (!requireNamespace("SC.MEB", quietly = TRUE)) {
    # Use GitHub's archive endpoint rather than install_github(). The latter
    # queries the GitHub API first and can fail in shared/CI environments when
    # the unauthenticated API rate limit is exhausted.
    scmeb_url <- paste0(
      "https://github.com/Shufeyangyi2015310117/SC.MEB/",
      "archive/refs/heads/master.tar.gz"
    )
    msg("Installing SC.MEB from its source archive...")
    install.packages(
      scmeb_url,
      repos = NULL,
      type = "source",
      lib = lib
    )
  }

  # IMPORTANT: QR-SIDE uses its own companion package named SpatialDecon.
  # This is different from the Bioconductor package with the same name.
  # The companion tarball is kept on main and is independent of QRSIDE ref.
  spatialdecon_url <- paste0(
    "https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/",
    "main/SpatialDecon_1.0.tar.gz"
  )

  spatialdecon_tar <- tempfile(fileext = ".tar.gz")
  on.exit(unlink(spatialdecon_tar), add = TRUE)

  msg("Installing the QR-SIDE companion SpatialDecon package...")
  utils::download.file(
    spatialdecon_url,
    spatialdecon_tar,
    mode = "wb",
    quiet = !isTRUE(verbose)
  )
  install.packages(
    spatialdecon_tar,
    repos = NULL,
    type = "source",
    lib = lib
  )

  if (!requireNamespace("SpatialDecon", quietly = TRUE) ||
      !exists(
        "CountDeconvolution",
        envir = asNamespace("SpatialDecon"),
        inherits = FALSE
      )) {
    stop(
      paste(
        "The required QR-SIDE companion SpatialDecon package was not installed",
        "correctly. In particular, CountDeconvolution() is missing.",
        "Do not substitute the Bioconductor package of the same name."
      )
    )
  }

  qrside_url <- paste0(
    "https://codeload.github.com/",
    repo,
    "/tar.gz/",
    utils::URLencode(ref, reserved = TRUE)
  )
  msg("Installing QRSIDE from its source archive...")
  install.packages(
    qrside_url,
    repos = NULL,
    type = "source",
    lib = lib
  )

  if (!requireNamespace("QRSIDE", quietly = TRUE)) {
    stop("QRSIDE installation failed.")
  }

  if (!exists("run_QRSIDE", envir = asNamespace("QRSIDE"), inherits = FALSE)) {
    stop("QRSIDE installed, but run_QRSIDE() was not found.")
  }

  msg("QR-SIDE installation completed successfully.")
  msg("Load it with: library(QRSIDE)")

  invisible(TRUE)
}
