# Run in a fresh R session:
# source("https://raw.githubusercontent.com/SpatialOmics-pj/QR-SIDE/main/install_qrside.R")
# install_qrside()
install_qrside <- function(repo = "SpatialOmics-pj/QR-SIDE", ref = "main",
                           lib = .libPaths()[1], upgrade = FALSE,
                           verbose = TRUE, package_file = NULL) {
  if (getRversion() < "4.2.3") stop("QR-SIDE requires R >= 4.2.3.")
  if ("QRSIDE" %in% loadedNamespaces())
    stop("Restart R/RStudio before installing or updating QRSIDE, then rerun the installer.")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (file.access(lib, 2L) != 0L) {
    user_lib <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1L]]
    if (!length(user_lib) || !nzchar(user_lib[1L]))
      stop("Choose a writable R library using install_qrside(lib = '/path/to/library').")
    lib <- path.expand(user_lib[1L])
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib <- normalizePath(lib, mustWork = TRUE)
  if (file.access(lib, 2L) != 0L) stop("The selected R library is not writable: ", lib)
  .libPaths(unique(c(lib, .libPaths())))
  old_timeout <- getOption("timeout", 60)
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  repos <- getOption("repos")
  if (!length(repos["CRAN"]) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
    repos["CRAN"] <- "https://cloud.r-project.org"
    options(repos = repos)
  }
  msg <- function(...) if (isTRUE(verbose)) message(...)
  available <- function(p) requireNamespace(p, quietly = TRUE)
  cran <- c("BH", "Rcpp", "RcppArmadillo", "clue", "ggplot2",
            "mclust", "quadprog", "reshape2", "scatterpie", "slam",
            "topicmodels", "viridis")
  missing <- cran[!vapply(cran, available, logical(1))]
  if (length(missing)) {
    msg("Installing CRAN dependencies: ", paste(missing, collapse = ", "))
    # Use binary dependencies on Windows/macOS when supplied for this R version.
    # QRSIDE itself is installed from source and still needs the R build tools.
    if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
      binary <- tryCatch(rownames(available.packages(type = "binary")), error = function(e) character())
      ready <- intersect(missing, binary)
      if (length(ready)) install.packages(ready, lib = lib, type = "binary")
      missing <- missing[!vapply(missing, available, logical(1))]
    }
    if (length(missing)) install.packages(missing, lib = lib)
  }
  if (!available("BiocManager")) install.packages("BiocManager", lib = lib)
  bioc <- c("BiocParallel", "fgsea", "SingleCellExperiment", "SummarizedExperiment", "scuttle", "scran")
  missing <- bioc[!vapply(bioc, available, logical(1))]
  if (length(missing))
    BiocManager::install(missing, lib = lib, ask = FALSE, update = isTRUE(upgrade), force = TRUE)
  if (!available("STdeconvolve") || utils::packageVersion("STdeconvolve") < "1.3.1") {
    msg("Installing STdeconvolve...")
    install.packages("https://codeload.github.com/JEFworks-Lab/STdeconvolve/tar.gz/f37bf14f39df995b4bd9dab402a763d18c50abeb",
                     repos = NULL, type = "source", lib = lib)
  }
  required <- c(cran, bioc, "STdeconvolve")
  missing <- required[!vapply(required, available, logical(1))]
  if (length(missing)) stop("Dependencies could not be loaded: ", paste(missing, collapse = ", "),
                            ". Review the installation output and system requirements in README.")

  if (is.null(package_file)) {
    # Download the small R source package, not the large example-data repository.
    base <- paste0("https://raw.githubusercontent.com/", repo, "/", ref, "/")
    description_file <- tempfile(fileext = ".dcf")
    package_file <- tempfile(fileext = ".tar.gz")
    on.exit(unlink(c(description_file, package_file)), add = TRUE)
    utils::download.file(paste0(base, "DESCRIPTION"), description_file, quiet = !verbose)
    version <- read.dcf(description_file, fields = "Version")[[1L]]
    utils::download.file(paste0(base, "QRSIDE_", version, ".tar.gz"),
                         package_file, mode = "wb", quiet = !verbose)
  }
  if (!file.exists(package_file)) stop("Source package not found: ", package_file)
  msg("Installing QRSIDE...")
  install.packages(package_file, repos = NULL, type = "source", lib = lib)
  if (!requireNamespace("QRSIDE", quietly = TRUE) || utils::packageVersion("QRSIDE") < "1.0.1")
    stop("QRSIDE installation did not complete. Review the compilation output above.")
  ns <- asNamespace("QRSIDE")
  stopifnot(is.function(QRSIDE::run_QRSIDE),
            is.function(get(".count_deconvolution", ns)),
            is.function(get(".joint_st_count_decon", ns)))
  msg("QR-SIDE installed successfully. Load it with library(QRSIDE).")
  invisible(TRUE)
}
