#' Run QR-SIDE spatial cell-type deconvolution
#' @param sp_counts Nonnegative integer counts, with named spots in rows and genes in columns.
#' @param pos Two numeric coordinate columns with the same spot identifiers.
#' @param markerframe A gene/label data frame, or a named list of marker vectors.
#' @param Num_topic Number of spatial topics.
#' @param Num_HVG Number of highly variable genes selected before retaining markers.
#' @param dim_embed Latent embedding dimension.
#' @param max_marker Maximum markers per cell type for table input.
#' @param markerflag Use a marker list (TRUE) or table (FALSE).
#' @param max_iter Maximum iterations in each model-fitting stage (at least 2).
#' @param seed Random seed for reproducible initialization.
#' @param verbose Print model diagnostics.
#' @return A list containing pi, beta, order, pos, and cell_types.
#' @export
run_QRSIDE <- function(sp_counts, pos, markerframe, Num_topic, Num_HVG,
                       dim_embed, max_marker, markerflag, max_iter = 50L,
                       seed = 999L, verbose = FALSE) {
  input <- .validate_spatial_input(sp_counts, pos)
  X <- input$counts
  pos <- input$pos
  if (nrow(X) < 3L || any(rowSums(X) == 0))
    stop("At least three spots with positive total counts are required.", call. = FALSE)
  X <- X[, colSums(X) > 0, drop = FALSE]
  .positive_integer(Num_topic, "Num_topic")
  .positive_integer(Num_HVG, "Num_HVG")
  .positive_integer(dim_embed, "dim_embed")
  .positive_integer(max_marker, "max_marker")
  .positive_integer(max_iter, "max_iter")
  if (max_iter < 2L) stop("max_iter must be at least 2.", call. = FALSE)
  if (Num_topic > nrow(X)) stop("Num_topic exceeds the number of spots.", call. = FALSE)
  if (Num_HVG > ncol(X)) stop("Num_HVG exceeds the number of genes.", call. = FALSE)
  if (!is.logical(markerflag) || length(markerflag) != 1L || is.na(markerflag))
    stop("markerflag must be TRUE or FALSE.", call. = FALSE)
  .positive_integer(seed, "seed", allow_zero = TRUE)
  set.seed(seed)

  markers <- .prepare_markers(markerframe, markerflag, colnames(X), max_marker)
  cell_types <- names(markers)
  K <- length(markers)
  if (K < 2L) stop("At least two cell types with usable markers are required.", call. = FALSE)
  if (K > min(dim(X))) stop("Too many cell types for the count matrix.", call. = FALSE)
  marker_genes <- unlist(markers, use.names = FALSE)

  # Marker genes must remain available even when not selected as HVGs.
  sce <- Createsce(X, nrow(X))
  sce <- scuttle::logNormCounts(sce)
  dec <- scran::modelGeneVar(sce)
  hvgs <- scran::getTopHVGs(dec, n = Num_HVG)
  genes <- unique(c(hvgs, marker_genes))
  logX <- t(as.matrix(SummarizedExperiment::assay(sce, "logcounts")[genes, , drop = FALSE]))
  nonmarkers <- setdiff(genes, marker_genes)
  variable <- vapply(nonmarkers, function(g) stats::sd(logX[, g]) > 0, logical(1))
  nonmarkers <- nonmarkers[variable]
  if (dim_embed > min(length(nonmarkers), nrow(X) - 1L))
    stop("dim_embed must not exceed the number of variable non-marker genes or spots minus one.", call. = FALSE)
  X_u <- X[, nonmarkers, drop = FALSE]
  X_m <- X[, marker_genes, drop = FALSE]
  if (any(rowSums(X_u) == 0) || any(rowSums(X_m) == 0))
    stop("Every spot needs positive marker and non-marker counts; filter low-quality spots before fitting.", call. = FALSE)
  n <- nrow(X)
  p_m <- ncol(X_m)
  p <- ncol(X_u)
  q <- dim_embed
  rho <- matrix(0, p_m, K)
  rho[cbind(seq_len(p_m), rep(seq_len(K), lengths(markers)))] <- 1
  delta <- rho * 3

  # Retain the caller's spot set and order throughout initialization.
  clean <- STdeconvolve::cleanCounts(t(X), min.lib.size = 0,
                                    min.reads = 0, min.detected = 0)
  corpus <- tryCatch(
    STdeconvolve::restrictCorpus(clean, removeAbove = 1, removeBelow = 0.05,
                                verbose = verbose),
    error = function(e) NULL
  )
  lda_genes <- if (is.null(corpus)) character() else rownames(corpus)
  lda_genes <- unique(c(lda_genes, marker_genes))
  lda_counts <- X[, lda_genes, drop = FALSE]
  if (ncol(lda_counts) < K || any(rowSums(lda_counts) == 0))
    stop("Insufficient expressed genes for LDA initialization after corpus filtering.", call. = FALSE)
  # Only one K is fitted, so model selection and parallel workers are unnecessary.
  lda <- topicmodels::LDA(slam::as.simple_triplet_matrix(lda_counts), k = K,
                         control = list(seed = as.integer(seed), verbose = 0,
                                        keep = 1, estimate.alpha = TRUE))
  result <- STdeconvolve::getBetaTheta(lda, perc.filt = 0.05,
                                     betaScale = 1000, verbose = verbose)
  deconProp <- result$theta[rownames(X), , drop = FALSE]
  markerexpr <- t(result$beta[, marker_genes, drop = FALSE])
  rho_lda <- matrix(0, p_m, K)
  rho_lda[cbind(seq_len(p_m), max.col(markerexpr, ties.method = "first"))] <- 1

  if (Num_topic == 1L) {
    labels <- rep(1L, n)
  } else {
    gmm <- mclust::Mclust(deconProp, G = Num_topic, verbose = verbose)
    if (is.null(gmm) || length(unique(gmm$classification)) != Num_topic)
      stop("Spatial-topic initialization failed; reduce Num_topic or check the marker/count variation.", call. = FALSE)
    labels <- gmm$classification
  }
  beta <- t(vapply(seq_len(Num_topic), function(i)
    colMeans(deconProp[labels == i, , drop = FALSE]), numeric(K)))
  beta <- beta / rowSums(beta)
  factors <- factpc(logX[, nonmarkers, drop = FALSE], q)
  Z_int <- factors$scores
  W_int <- factors$loadings
  pi_initial <- matrix(1e-6, n, Num_topic)
  pi_initial[cbind(seq_len(n), labels)] <- 1
  pi_initial <- pi_initial / rowSums(pi_initial)
  pi0 <- rep(1 / Num_topic, Num_topic)
  mu <- .initial_mu(Z_int, labels, beta, K)

  initial <- .count_deconvolution(
    X_u, Z_int, W_int, K = K, T = Num_topic, beta = beta, mu = t(mu),
    pi = pi_initial, pi0 = matrix(pi0), Lambda = rep(1, p), alpha = rep(0, n),
    epsilon2 = matrix(1, n, q), max_iter = max_iter, eps = 1e-5, verbose = verbose,
    fixed_loglambda_mu = FALSE, fixed_loglambda_s2 = FALSE,
    fixed_z_mu = FALSE, fixed_z_s2 = FALSE, fixed_m = TRUE, fixed_b = TRUE,
    fixed_pi = TRUE, fixed_W = FALSE, fixed_alpha = FALSE, fixed_Lambda = FALSE,
    fixed_pi0 = TRUE, fixed_epsilon2 = TRUE, fixed_mu = TRUE,
    fixed_sigma2 = TRUE, fixed_beta2 = TRUE, alpha_spot_specific = TRUE, cutoff = 1e6
  )

  distances <- outer(seq_len(K), seq_len(K), Vectorize(function(i, j)
    sqrt(sum((rho[, i] - rho_lda[, j])^2))))
  permutation <- order(as.integer(clue::solve_LSAP(distances)))
  rho <- rho[, permutation, drop = FALSE]
  delta <- delta[, permutation, drop = FALSE]
  mu <- .initial_mu(initial$z_mu, labels, beta, K)
  adjacency <- .spatial_neighbors(as.matrix(pos), cutoff = 1.3)
  y <- max.col(initial$pi, ties.method = "first")
  marker_fit <- DeconvMarker(X_m, rho, Num_topic, beta, delta, pi_initial,
                            matrix(pi0), rep(1, p_m), matrix(0, p_m, K),
                            max_iter = max_iter, eps = 1e-5, verbose = verbose)
  m <- .initial_mu(marker_fit$loglambda_mu, labels, beta, K)
  fitted <- .joint_st_count_decon(
    X_m, rho, X_u, y, adjacency, initial$z_mu, initial$W,
    K = K, T = Num_topic, beta = beta, mu = t(mu), pi = pi_initial, pi0 = pi0,
    Lambda = initial$Lambda, alpha = initial$alpha, epsilon2 = matrix(1, n, q),
    beta0_y = rep(0, Num_topic), beta_y_grid = seq(0, 4, 0.2),
    Lambda_m = rep(1, p_m), delta = delta, m = m,
    max_iter = max_iter, maxIter_ICM = 50L, eps = 1e-4, verbose = verbose,
    fixed_loglambda_mu = FALSE, fixed_loglambda_s2 = FALSE,
    fixed_z_mu = FALSE, fixed_z_s2 = FALSE, fixed_m = FALSE, fixed_b = FALSE,
    fixed_pi = FALSE, fixed_W = FALSE, fixed_alpha = FALSE, fixed_Lambda = FALSE,
    fixed_pi0 = FALSE, fixed_epsilon2 = FALSE, fixed_mu = FALSE,
    fixed_sigma2 = FALSE, fixed_beta2 = FALSE, fixed_beta_y = FALSE,
    alpha_spot_specific = TRUE, cutoff = 300
  )
  if (any(!is.finite(fitted$pi)) || any(!is.finite(fitted$beta)))
    stop("Model fitting produced non-finite values; check counts and model dimensions.", call. = FALSE)
  rownames(fitted$pi) <- rownames(X)
  colnames(fitted$pi) <- rownames(fitted$beta) <- paste0("topic", seq_len(Num_topic))
  colnames(fitted$beta) <- cell_types[permutation]
  list(pi = fitted$pi, beta = fitted$beta, order = permutation, pos = pos,
       cell_types = cell_types[permutation])
}

.initial_mu <- function(Z, labels, beta, K) {
  means <- do.call(rbind, lapply(seq_len(nrow(beta)), function(i)
    colMeans(Z[labels == i, , drop = FALSE])))
  t(solve(crossprod(beta) + diag(1e-5, K), crossprod(beta, means)))
}
