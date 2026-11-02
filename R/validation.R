.positive_integer <- function(x, name, allow_zero = FALSE) {
  minimum <- if (allow_zero) 0 else 1
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max)
    stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
  invisible(TRUE)
}

.validate_spatial_input <- function(counts, pos) {
  if (!(is.matrix(counts) || is.data.frame(counts) || inherits(counts, "Matrix")))
    stop("Counts must be a matrix, data frame, or Matrix object.", call. = FALSE)
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || !all(dim(counts) > 0L) ||
      any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)))
    stop("Counts must contain finite nonnegative integers (raw counts).", call. = FALSE)
  for (margin in c("row", "column")) {
    ids <- if (margin == "row") rownames(counts) else colnames(counts)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("Counts need unique, nonempty ", margin, " names.", call. = FALSE)
  }
  pos <- as.matrix(pos)
  if (!is.numeric(pos) || ncol(pos) != 2L || any(!is.finite(pos)))
    stop("pos must have exactly two finite numeric coordinate columns.", call. = FALSE)
  ids <- rownames(pos)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !setequal(ids, rownames(counts)))
    stop("pos must have the same unique spot row names as counts.", call. = FALSE)
  storage.mode(counts) <- "double"
  list(counts = counts, pos = as.data.frame(pos[rownames(counts), , drop = FALSE]))
}

.prepare_markers <- function(markers, is_list, genes, max_marker) {
  if (is_list) {
    if (!is.list(markers) || is.data.frame(markers) || !length(markers))
      stop("markerflag = TRUE requires a list of marker vectors.", call. = FALSE)
    if (is.null(names(markers))) names(markers) <- paste0("cell_type", seq_along(markers))
  } else {
    if (!is.data.frame(markers) || !all(c("gene", "label") %in% names(markers)))
      stop("Marker tables require gene and label columns.", call. = FALSE)
    labels <- as.character(markers$label)
    if (anyNA(labels) || any(!nzchar(labels)))
      stop("Marker cell-type labels must be nonempty.", call. = FALSE)
    types <- unique(labels)
    markers <- setNames(lapply(types, function(x) as.character(markers$gene[labels == x])), types)
  }
  if (anyNA(names(markers)) || any(!nzchar(names(markers))) || anyDuplicated(names(markers)))
    stop("Marker list names must be unique, nonempty cell-type labels.", call. = FALSE)
  used <- character()
  for (i in seq_along(markers)) {
    values <- as.character(markers[[i]])
    if (anyNA(values) || any(!nzchar(values)))
      stop("Marker genes must be nonempty strings.", call. = FALSE)
    values <- setdiff(intersect(unique(values), genes), used)
    if (!is_list) values <- head(values, max_marker)
    if (!length(values))
      stop("No expressed, distinct marker genes remain for cell type ", names(markers)[i], ".", call. = FALSE)
    markers[[i]] <- values
    used <- c(used, values)
  }
  markers
}
