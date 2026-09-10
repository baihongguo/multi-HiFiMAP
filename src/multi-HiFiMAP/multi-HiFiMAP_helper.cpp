// ADD THIS LINE AT THE VERY TOP:
#define ARMA_64BIT_WORD

// SAVE AS: HiFiMAP_Stateful.cpp
#include <RcppArmadillo.h>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]

class HiFiMAPCalculator {
private:
    arma::sp_mat X;          // Persists in memory
    arma::mat R_mat;         // Persists in memory
    arma::mat K_inv;         // Persists in memory
    
    // Cached statistics
    arma::vec c_means;       // Column means
    double Tstat_num;        // Numerator of Tstat
    double n_dbl;            // Dimension as double
    int n;                   // Dimension as int

public:
    // Constructor
    HiFiMAPCalculator(arma::sp_mat X_init, arma::mat R_init, arma::mat K_init) 
        : X(X_init), R_mat(R_init), K_inv(K_init) {
        
        n = X.n_rows;
        n_dbl = (double)n;

        // 1. Initial Means
        arma::mat col_sums_mat = arma::mat(arma::sum(X, 0));
        c_means = col_sums_mat.t() / n_dbl;

        // 2. Initial Tstat_num
        arma::mat temp_prod = R_mat.t() * X; 
        arma::mat p_by_p_matrix = temp_prod * R_mat;
        Tstat_num = arma::accu(p_by_p_matrix % K_inv.t());
    }

    // NEW HELPER: Perform X * v using the internal stored X
    // This allows R to get products without holding the full matrix X
    arma::vec multiply_X_vec(const arma::vec& v) {
        return X * v;
    }

    // The update function
    Rcpp::List update_and_calculate(const arma::sp_mat& dX, double T3_num) {
        auto start_total = std::chrono::high_resolution_clock::now();

        // 1. Update Means
        arma::mat dX_col_sums = arma::mat(arma::sum(dX, 0));
        arma::vec d_means = dX_col_sums.t() / n_dbl;
        c_means += d_means;
        double g_mean = arma::mean(c_means);
        
        // 2. Update Tstat_num
        arma::mat d_temp_prod = R_mat.t() * dX;
        arma::mat d_p_by_p = d_temp_prod * R_mat;
        double dTstat_num = arma::accu(d_p_by_p % K_inv.t());
        Tstat_num += dTstat_num;

        // 3. Update X Internally (The only place X grows)
        X = X + dX; 

        // 4. Calculate Denominator (sum_Xc_sq)
        // Note: X % X is efficient for sparse matrices
        double sum_sq_X = arma::accu(X % X);
        double sum_Xc_sq = sum_sq_X - 2.0 * n_dbl * arma::accu(arma::pow(c_means, 2)) + n_dbl * n_dbl * g_mean * g_mean;
        
        double Tstat = Tstat_num / sum_Xc_sq;
        
        double sum_Xc_sq_2 = sum_Xc_sq * sum_Xc_sq;
        double sum_Xc_sq_3 = sum_Xc_sq_2 * sum_Xc_sq;

        // 5. Calculate T, T2, S2, S3
        double tr_X = arma::accu(X.diag());
        double tr_Xc = tr_X - n_dbl * g_mean;
        double T = tr_Xc / sum_Xc_sq;
        double T2 = 1.0 / sum_Xc_sq; 

        // Explicit cast to dense vector for diag operations
        arma::vec X_diag_dense = arma::vec(X.diag());
        arma::vec diag_Xc = X_diag_dense - 2.0 * c_means + g_mean;

        double S2 = arma::accu(arma::pow(diag_Xc, 2)) / sum_Xc_sq_2;
        double S3 = arma::accu(arma::pow(diag_Xc, 3)) / sum_Xc_sq_3;

        // 6. Calculate U
        double sum_c1 = arma::accu(c_means);
        double sum_c2 = arma::accu(arma::pow(c_means, 2));
        double sum_c3 = arma::accu(arma::pow(c_means, 3));
        
        double sum_Xc_cubed_zeros = -n_dbl*sum_c3 - 3*sum_c1*sum_c2 + 3*g_mean*n_dbl*sum_c2 - 3*g_mean*sum_c1*sum_c1 + 3*g_mean*g_mean*n_dbl*sum_c1 - n_dbl*n_dbl*g_mean*g_mean*g_mean;
        double sum_Xc_cubed_nonzeros = 0;

        arma::vec r_means = c_means; 
        for(arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
            int r = it.row();
            int c = it.col();
            double x_ij = *it;
            double term_common = -r_means(r) - c_means(c) + g_mean;
            double xc_ij = x_ij + term_common;
            double xc_ij_zero = term_common;
            sum_Xc_cubed_nonzeros += (xc_ij*xc_ij*xc_ij) - (xc_ij_zero*xc_ij_zero*xc_ij_zero);
        }
        double U = (sum_Xc_cubed_zeros + sum_Xc_cubed_nonzeros) / sum_Xc_sq_3;

        // 7. Calculate B
        double B_num = arma::as_scalar(diag_Xc.t() * (X * diag_Xc))
                     - 2 * arma::dot(diag_Xc, (diag_Xc % c_means))
                     + g_mean * arma::accu(arma::pow(diag_Xc, 2));
        double B = B_num / sum_Xc_sq_3;

        // 8. Calculate R
        arma::sp_mat X_sq = X % X; 
        arma::vec diag_X_sq = arma::vec(X_sq.diag()); // Explicit cast
        arma::vec diag_Xc_sq = diag_X_sq 
                             - 2 * (X_diag_dense % c_means)
                             + arma::pow(c_means, 2)
                             - 2 * g_mean * (X_diag_dense - c_means)
                             + g_mean * g_mean;
        double R_num = arma::dot(diag_Xc, diag_Xc_sq);
        double R = R_num / sum_Xc_sq_3;

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