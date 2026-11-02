version <- read.dcf("DESCRIPTION", fields = "Version")[[1L]]
archive <- paste0("QRSIDE_", version, ".tar.gz")
stopifnot(file.exists(archive))
tmp <- tempfile("qrside-archive-")
dir.create(tmp)
utils::untar(archive, exdir = tmp)
files <- list.files(c("R", "src", "man", "tests", "inst/examples"),
                    recursive = TRUE, full.names = TRUE)
files <- files[!grepl("\\.(o|so|dll|dylib)$", files)]
files <- c("NAMESPACE", files)
for (file in files) {
  packaged <- file.path(tmp, "QRSIDE", file)
  if (!file.exists(packaged) || tools::md5sum(file) != tools::md5sum(packaged))
    stop("Download archive is stale: ", file)
}
fields <- c("Version", "Imports", "LinkingTo", "Depends")
normalize <- function(x) gsub("[[:space:]]+", " ", as.character(x))
stopifnot(identical(normalize(read.dcf("DESCRIPTION", fields = fields)),
                    normalize(read.dcf(file.path(tmp, "QRSIDE", "DESCRIPTION"), fields = fields))))
unlink(tmp, recursive = TRUE)
cat("Download archive matches current source.\n")
