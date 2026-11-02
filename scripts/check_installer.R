# A failed package install must never be reported as a successful old installation.
source("install_qrside.R")
broken <- tempfile(fileext = ".tar.gz")
writeLines("deliberately invalid source package", broken)
err <- tryCatch({install_qrside(package_file = broken, verbose = FALSE); NULL}, error = identity)
unlink(broken)
stopifnot(inherits(err, "error"), grepl("installation failed", conditionMessage(err), fixed = TRUE))
cat("PASS failed installation is detected even when QRSIDE is already installed.\n")
