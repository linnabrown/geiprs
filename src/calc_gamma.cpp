#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_gamma_Rcpp(NumericMatrix x, IntegerVector ix, IntegerVector iy, int bn) {
  // Re-implemented the sparsegl:::calc_gamma function using Rcpp. 
  // This modified version does not work for "ncols > 2" in the original code. 
  // Helper lambda to calculate largest squared singular value for 2-column case
  auto maxeig2 = [](NumericMatrix& mat, int col1, int col2) {
    double mat11 = 0.0, mat12 = 0.0, mat22 = 0.0;

    // Cross-product of the matrix (mat' * mat)
    for (int i = 0; i < mat.nrow(); ++i) {
      mat11 += mat(i, col1) * mat(i, col1);
      mat12 += mat(i, col1) * mat(i, col2);
      mat22 += mat(i, col2) * mat(i, col2);
    }

    // Trace and determinant
    double tr = mat11 + mat22;
    double det = mat11 * mat22 - mat12 * mat12;

    // Return the largest eigenvalue
    return (tr + sqrt(tr * tr - 4 * det)) / 2;
  };

  NumericVector gamma(bn);

  for (int g = 0; g < bn; ++g) {
    int start = ix[g] - 1; // Adjust 1-based index from R to 0-based for C++
    int end = iy[g] - 1;
    int ncols = end - start + 1;

    if (ncols == 2) {
      // Compute largest singular value directly without creating a submatrix
      gamma[g] = maxeig2(x, start, end);
    } else {
      // Sum of squares of the selected columns
      double sumsq = 0.0;
      for (int i = 0; i < x.nrow(); ++i) {
        for (int j = start; j <= end; ++j) {
          sumsq += x(i, j) * x(i, j);
        }
      }
      gamma[g] = sumsq;
    }
  }

  return gamma / x.nrow();
}