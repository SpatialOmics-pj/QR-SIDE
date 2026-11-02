#include <RcppArmadillo.h>

// Same strict Euclidean radius rule used by the original SC.MEB dependency.
// [[Rcpp::export(name = ".spatial_neighbors")]]
arma::sp_mat spatial_neighbors(const arma::mat& positions, double cutoff = 1.3) {
  arma::sp_mat adjacency(positions.n_rows, positions.n_rows);
  const double radius2 = cutoff * cutoff;
  for (arma::uword j = 0; j < positions.n_rows; ++j) {
    for (arma::uword i = j + 1; i < positions.n_rows; ++i) {
      const double dx = positions(i, 0) - positions(j, 0);
      if (std::abs(dx) >= cutoff) continue;
      const double dy = positions(i, 1) - positions(j, 1);
      if (dx * dx + dy * dy < radius2) {
        adjacency(i, j) = 1;
        adjacency(j, i) = 1;
      }
    }
  }
  return adjacency;
}
