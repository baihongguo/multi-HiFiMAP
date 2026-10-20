// multi-HiFiMAP C++ backend (version 1.1.0)
// Loaded by run_multi_HiFiMAP.R through Rcpp::sourceCpp(); no manual compilation is needed.
//
// HiFiMAPCalculator keeps the sparse local IBD matrix X of the current SNP in memory and, for each new SNP,
// applies the sparse difference dX and returns the IBD-side moment terms of the RV-coefficient test
// (T, T2, S2, S3, U, R, B, T3) computed from the non-zero elements of X only; the double-centered matrix
// HXH is never formed.

// 64-bit indices so that very dense IBD matrices (more than 2^31 non-zeros) are supported
#define ARMA_64BIT_WORD

#include <RcppArmadillo.h>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]

class HiFiMAPCalculator {
private:
    arma::sp_mat X;          // local IBD matrix of the current SNP
    arma::mat R_mat;         // n x k scaled residuals
    arma::mat K_inv;         // k x k inverse phenotype correlation matrix

    arma::vec c_means;       // column means of X
    double Tstat_num;        // numerator of the test statistic, tr(r' X r K^-1)
    double n_dbl;
    int n;

public:
    HiFiMAPCalculator(arma::sp_mat X_init, arma::mat R_init, arma::mat K_init)
        : X(X_init), R_mat(R_init), K_inv(K_init) {

        n = X.n_rows;
        n_dbl = (double)n;

        arma::mat col_sums_mat = arma::mat(arma::sum(X, 0));
        c_means = col_sums_mat.t() / n_dbl;

        arma::mat temp_prod = R_mat.t() * X;
        arma::mat p_by_p_matrix = temp_prod * R_mat;
        Tstat_num = arma::accu(p_by_p_matrix % K_inv.t());
    }

    // X * v with the stored X, so that R does not need a second copy of X for this product
    arma::vec multiply_X_vec(const arma::vec& v) {
        return X * v;
    }

    // Apply dX (the change in the IBD matrix from the previous SNP) and return the moment terms.
    // T3_num is the change in tr((HXH)^3), estimated in R with Hutchinson's estimator.
    Rcpp::List update_and_calculate(const arma::sp_mat& dX, double T3_num) {
        auto start_total = std::chrono::high_resolution_clock::now();

        // 1. Update column means and the test statistic numerator
        arma::mat dX_col_sums = arma::mat(arma::sum(dX, 0));
        arma::vec d_means = dX_col_sums.t() / n_dbl;
        c_means += d_means;
        double g_mean = arma::mean(c_means);

        arma::mat d_temp_prod = R_mat.t() * dX;
        arma::mat d_p_by_p = d_temp_prod * R_mat;
        Tstat_num += arma::accu(d_p_by_p % K_inv.t());

        // 2. Update X
        X = X + dX;

        // 3. Sum of squares of HXH from one pass over the non-zeros (no copy of X);
        //    the same pass collects the row sums of squares used for R below
        arma::vec rowsq(n, arma::fill::zeros);
        double sum_sq_X = 0.0;
        for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
            double x_sq = (*it) * (*it);
            rowsq(it.row()) += x_sq;
            sum_sq_X += x_sq;
        }
        double sum_Xc_sq = sum_sq_X - 2.0 * n_dbl * arma::accu(arma::pow(c_means, 2)) + n_dbl * n_dbl * g_mean * g_mean;

        double Tstat = Tstat_num / sum_Xc_sq;
        double sum_Xc_sq_2 = sum_Xc_sq * sum_Xc_sq;
        double sum_Xc_sq_3 = sum_Xc_sq_2 * sum_Xc_sq;

        // 4. T, T2, S2, S3
        double tr_X = arma::accu(X.diag());
        double tr_Xc = tr_X - n_dbl * g_mean;
        double T = tr_Xc / sum_Xc_sq;
        double T2 = 1.0 / sum_Xc_sq;

        arma::vec X_diag_dense = arma::vec(X.diag());
        arma::vec diag_Xc = X_diag_dense - 2.0 * c_means + g_mean;

        double S2 = arma::accu(arma::pow(diag_Xc, 2)) / sum_Xc_sq_2;
        double S3 = arma::accu(arma::pow(diag_Xc, 3)) / sum_Xc_sq_3;

        // 5. U = sum_ij (HXH)_ij^3: zero entries of X contribute (g - c_i - c_j)^3, which sums to
        //    n^2 g^3 - 2 n sum_i c_i^3; non-zero entries are corrected one by one
        double sum_c2 = arma::accu(arma::pow(c_means, 2));
        double sum_c3 = arma::accu(arma::pow(c_means, 3));
        double sum_Xc_cubed_zeros = n_dbl * n_dbl * g_mean * g_mean * g_mean - 2.0 * n_dbl * sum_c3;
        double sum_Xc_cubed_nonzeros = 0;
        for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
            double term_common = -c_means(it.row()) - c_means(it.col()) + g_mean;
            double xc_ij = *it + term_common;
            sum_Xc_cubed_nonzeros += (xc_ij * xc_ij * xc_ij) - (term_common * term_common * term_common);
        }
        double U = (sum_Xc_cubed_zeros + sum_Xc_cubed_nonzeros) / sum_Xc_sq_3;

        // 6. B = d' (HXH) d with d = diag(HXH)
        double sum_d = arma::accu(diag_Xc);
        double B_num = arma::as_scalar(diag_Xc.t() * (X * diag_Xc))
                     - 2.0 * sum_d * arma::dot(c_means, diag_Xc)
                     + g_mean * sum_d * sum_d;
        double B = B_num / sum_Xc_sq_3;

        // 7. R = d' diag((HXH)^2), with diag((HXH)^2)_i = sum_j (X_ij - c_i - c_j + g)^2 for symmetric X
        arma::vec Xc_vec = X * c_means;
        arma::vec diag_Xc_sq = rowsq - 2.0 * Xc_vec - n_dbl * arma::square(c_means)
                             + 2.0 * n_dbl * g_mean * c_means + sum_c2 - n_dbl * g_mean * g_mean;
        double R = arma::dot(diag_Xc, diag_Xc_sq) / sum_Xc_sq_3;

        double T3 = T3_num / sum_Xc_sq_3;

        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_total - start_total;

        return Rcpp::List::create(
            Rcpp::Named("Tstat") = Tstat,
            Rcpp::Named("T") = T,
            Rcpp::Named("T2") = T2,
            Rcpp::Named("S2") = S2,
            Rcpp::Named("S3") = S3,
            Rcpp::Named("U") = U,
            Rcpp::Named("R") = R,
            Rcpp::Named("B") = B,
            Rcpp::Named("T3") = T3,
            Rcpp::Named("time_total_s") = elapsed.count()
        );
    }
};

RCPP_MODULE(HiFiMAP_Module) {
    Rcpp::class_<HiFiMAPCalculator>("HiFiMAPCalculator")
        .constructor<arma::sp_mat, arma::mat, arma::mat>()
        .method("update_and_calculate", &HiFiMAPCalculator::update_and_calculate)
        .method("multiply_X_vec", &HiFiMAPCalculator::multiply_X_vec)
    ;
}
