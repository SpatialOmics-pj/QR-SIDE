# Correlation-based factor scores used by the original QR-SIDE initializer.
factpc <- function(X, m = 2) {
  X <- as.matrix(X)
  S <- stats::cor(X)
  eig <- eigen(S, symmetric = TRUE)
  A <- sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), `*`)
  signs <- ifelse(colSums(A) < 0, -1, 1)
  A <- sweep(A, 2L, signs, `*`)[, seq_len(m), drop = FALSE]
  scores <- scale(X) %*% solve(S + diag(1e-6, ncol(X))) %*% A
  list(loadings = A, scores = scores)
}

#' @title 基因表达数据质量控制过滤
#' @description
#'   对空间转录组数据的基因表达矩阵进行质量控制过滤，基于基因和位点的总表达量阈值筛选数据，
#'   并保留与过滤后数据匹配的空间坐标信息。
#'
#' @param geneexpr 矩阵或数据框，形状为 "spot by gene"，即行代表空间位点，列代表基因。
#' @param pos 矩阵或数据框，形状为 "spot by 2"，即每一行对应一个空间位点的二维坐标（例如 x 和 y）。
#' @param threshold 数值，过滤阈值。行（位点）或列（基因）的总表达量必须大于此值才会被保留。
#'
#' @return 返回一个包含两个元素的列表：
#'   \itemize{
#'     \item \code{[[1]]} 矩阵，过滤后的基因表达矩阵（行和列均满足总表达量 > threshold）。
#'     \item \code{[[2]]} 数据框，与过滤后表达矩阵行名匹配的空间坐标信息，形状为 "spot by 2"。
#'   }
#' @export

QRSIDE_QC <- function(geneexpr, pos, threshold) {
  input <- .validate_spatial_input(geneexpr, pos)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0)
    stop("threshold must be a finite nonnegative number.", call. = FALSE)
  keep_spots <- rowSums(input$counts) > threshold
  keep_genes <- colSums(input$counts) > threshold
  filtered <- input$counts[keep_spots, keep_genes, drop = FALSE]
  if (!all(dim(filtered) > 0L))
    stop("QC removed every spot or every gene; lower threshold or check the counts.", call. = FALSE)
  list(filtered, input$pos[rownames(filtered), , drop = FALSE])
}


Createsce <- function(X, n) {
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = t(X)))
}
