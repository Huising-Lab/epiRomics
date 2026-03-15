// greedy_align.cpp — Rcpp implementation of greedy distance-based peak alignment
// Used by .align_tf_peaks() in epiRomics_hic.R
//
// Computes the distance delta matrix between enhancer and contact peak midpoints,
// then performs greedy matching (closest pairs first, each peak matched at most once).
//
// 10-50x faster than the R implementation by:
//   1. Avoiding base::outer() allocation for the full matrix
//   2. Computing matrix minimum in-place without Inf masking
//   3. No R interpreter overhead per iteration

#include <Rcpp.h>
using namespace Rcpp;

//' Greedy distance-based peak alignment (C++ implementation)
//'
//' Given enhancer distances and contact distances from their respective
//' folding boundaries, finds greedy matching of closest pairs within tolerance.
//'
//' @param enh_dists numeric vector of enhancer peak distances from boundary
//' @param con_dists numeric vector of contact peak distances from boundary
//' @param tolerance maximum delta for a pair to be considered aligned
//' @param keylock_tol tighter tolerance for key-lock designation
//' @return List with: n_aligned, aligned_enh (1-based), aligned_con (1-based),
//'         best_delta, is_key_lock
// [[Rcpp::export]]
List greedy_align_cpp(NumericVector enh_dists, NumericVector con_dists,
                      double tolerance, double keylock_tol) {
  int n_enh = enh_dists.size();
  int n_con = con_dists.size();
  int max_pairs = std::min(n_enh, n_con);

  // Pre-allocate output vectors
  IntegerVector aligned_enh(max_pairs);
  IntegerVector aligned_con(max_pairs);
  int n_aligned = 0;
  double min_delta = R_PosInf;

  // Matched flags
  std::vector<bool> matched_enh(n_enh, false);
  std::vector<bool> matched_con(n_con, false);

  for (int iter = 0; iter < max_pairs; iter++) {
    // Find unmatched pair with smallest delta
    double best_val = R_PosInf;
    int best_i = -1, best_j = -1;

    for (int i = 0; i < n_enh; i++) {
      if (matched_enh[i]) continue;
      for (int j = 0; j < n_con; j++) {
        if (matched_con[j]) continue;
        double delta = std::abs(enh_dists[i] - con_dists[j]);
        if (delta < best_val) {
          best_val = delta;
          best_i = i;
          best_j = j;
        }
      }
    }

    if (best_val > tolerance || best_i < 0) break;

    matched_enh[best_i] = true;
    matched_con[best_j] = true;
    aligned_enh[n_aligned] = best_i + 1;  // 1-based for R
    aligned_con[n_aligned] = best_j + 1;
    n_aligned++;

    if (best_val < min_delta) min_delta = best_val;
  }

  // Trim to actual size
  IntegerVector out_enh(n_aligned);
  IntegerVector out_con(n_aligned);
  for (int k = 0; k < n_aligned; k++) {
    out_enh[k] = aligned_enh[k];
    out_con[k] = aligned_con[k];
  }

  return List::create(
    Named("n_aligned") = n_aligned,
    Named("aligned_enh") = out_enh,
    Named("aligned_con") = out_con,
    Named("best_delta") = (n_aligned > 0) ? min_delta : R_PosInf,
    Named("is_key_lock") = (n_aligned > 0) && (min_delta <= keylock_tol)
  );
}
