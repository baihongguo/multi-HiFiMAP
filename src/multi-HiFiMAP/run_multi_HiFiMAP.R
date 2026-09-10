Sys.setenv(MKL_NUM_THREADS = 1)
library(Rcpp)
library(PearsonDS)
library(data.table)
library(Matrix)
library(RcppArmadillo)
library(MASS)

# Load the optimized C++ stateful backend
sourceCpp("src/HiFiMAP_Stateful.cpp")

# -------------------------------------------------------------------------
# 1. SETUP AND CHUNK ARGUMENTS
# -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 9) {
  stop("Error: Arguments required: <chr_num> <chunk_start> <chunk_end> <chunk_id> <ibd_dir> <res_file> <cor_file> <id_file> <out_file>")
}

chr_num     <- as.numeric(args[1])
chunk_start <- as.numeric(args[2])
chunk_end   <- as.numeric(args[3])
chunk_id    <- as.numeric(args[4])
ibd_dir     <- args[5]
res_path    <- args[6]
cor_path    <- args[7]
id_path     <- args[8]
output_file <- args[9]

time_file <- gsub("\\.txt$", ".time", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

N_hutchinson <- 10

# -------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -------------------------------------------------------------------------
read_diff_as_dK <- function(diff_file, n_hap) {
  if (!file.exists(diff_file)) return(NULL)
  dat <- scan(diff_file, what = list(integer(), integer(), integer()), quiet = TRUE)
  ops <- dat[[1]]; rows <- dat[[2]] + 1; cols <- dat[[3]] + 1
  if (length(ops) == 0) return(NULL)
  vals <- ifelse(ops == 1, 1, -1)
  dt <- data.table(i = rows, j = cols, x = vals)
  dt <- dt[, .(x = sum(x)), by = .(i, j)]
  dt <- dt[x != 0]
  if (nrow(dt) == 0) return(NULL)
  dK <- sparseMatrix(i = dt$i, j = dt$j, x = dt$x, dims = c(n_hap, n_hap))
  return(forceSymmetric(dK, uplo = "U"))
}

# --- EXACT HUTCHINSON ESTIMATOR ---
estimate_dynamic_traces_hutchinson <- function(X_sparse, XZ, dX, dXZ) {
  dX_XZ <- dX %*% XZ
  tr_X_dX_X <- mean(diag(crossprod(XZ, dX_XZ)))
  X_dXZ <- X_sparse %*% dXZ
  tr_X2_dX <- mean(diag(crossprod(XZ, X_dXZ)))
  
  Z_dX2 <- crossprod(dXZ, dX)
  tr_dX2_X <- mean(diag(Z_dX2 %*% XZ))
  tr_dX_X_dX <- mean(diag(crossprod(dXZ, X_dXZ)))
  
  tr_dX3 <- mean(diag(Z_dX2 %*% dXZ))
  
  return(list(
    tr_X2dX = 2 * tr_X2_dX + tr_X_dX_X,
    tr_XdX2 = 2 * tr_dX2_X + tr_dX_X_dX,
    tr_dX3  = tr_dX3
  ))
}

hutchinson_trace_X3_initial <- function(X, Z) {
  dot_products <- diag(crossprod(Z, X %*% (X %*% (X %*% Z))))
  return(mean(dot_products))
}

calculate_p_value <- function(stats, n, stats_precalc) {
  Tstat <- stats$Tstat; T_val <- stats$T; T2 <- stats$T2
  S2 <- stats$S2; S3 <- stats$S3; U <- stats$U; R_val <- stats$R; B <- stats$B; T3 <- stats$T3
  Ts <- stats_precalc$Ts; T2s <- stats_precalc$T2s; S2s <- stats_precalc$S2s
  S3s <- stats_precalc$S3s; Us <- stats_precalc$Us; Rs <- stats_precalc$Rs
  Bs <- stats_precalc$Bs; T3s <- stats_precalc$T3s
  
  mean.rv <- T_val * Ts / (n - 1)
  term_A <- (n - 1) * T2 - T_val^2
  term_B <- (n - 1) * T2s - Ts^2
  if(term_A < 0) term_A <- 0
  if(term_B < 0) term_B <- 0
  
  temp1 <- 2 * term_A * term_B / (n - 1)^2 / (n + 1) / (n - 2)
  temp21 <- n * (n + 1) * S2 - (n - 1) * (T_val^2 + 2 * T2)
  temp22 <- n * (n + 1) * S2s - (n - 1) * (Ts^2 + 2 * T2s)
  temp23 <- (n + 1) * n * (n - 1) * (n - 2) * (n - 3)
  temp2 <- temp21 * temp22 / temp23
  variance.rv <- temp1 + temp2
  
  if (variance.rv <= 1e-20) return(1.0)
  
  t1 <- n^2 * (n + 1) * (n^2 + 15 * n - 4) * S3 * S3s
  t2 <- 4 * (n^4 - 8 * n^3 + 19 * n^2 - 4 * n - 16) * U * Us
  t3 <- 24 * (n^2 - n - 4) * (U * Bs + B * Us)
  t4 <- 6 * (n^4 - 8 * n^3 + 21 * n^2 - 6 * n - 24) * B * Bs
  t5 <- 12 * (n^4 - n^3 - 8 * n^2 + 36 * n - 48) * R_val * Rs
  t6 <- 12 * (n^3 - 2 * n^2 + 9 * n - 12) * (T_val * S2 * Rs + R_val * Ts * S2s)
  t7 <- 3 * (n^4 - 4 * n^3 - 2 * n^2 + 9 * n - 12) * T_val * Ts * S2 * S2s
  t81 <- (n^3 - 3 * n^2 - 2 * n + 8) * (R_val * Us + U * Rs); t82 <- (n^3 - 2 * n^2 - 3 * n + 12) * (R_val * Bs + B * Rs)
  t8 <- 24 * (t81 + t82)
  t9 <- 12 * (n^2 - n + 4) * (T_val * S2 * Us + U * Ts * S2s)
  t10 <- 6 * (2 * n^3 - 7 * n^2 - 3 * n + 12) * (T_val * S2 * Bs + B * Ts * S2s)
  t11 <- -2 * n * (n - 1) * (n^2 - n + 4) * ((2 * U + 3 * B) * S3s + (2 * Us + 3 * Bs) * S3)
  t12 <- -3 * n * (n - 1)^2 * (n + 4) * ((T_val * S2 + 4 * R_val) * S3s + (Ts * S2s + 4 * Rs) * S3)
  t13 <- 2 * n * (n - 1) * (n - 2) * ((T_val^3 + 6 * T_val * T2 + 8 * T3) * S3s + (Ts^3 + 6 * Ts * T2s + 8 * T3s) * S3)
  t14 <- T_val^3 * ((n^3 - 9 * n^2 + 23 * n - 14) * Ts^3 + 6 * (n - 4) * Ts * T2s + 8 * T3s)
  t15 <- 6 * T_val * T2 * ((n - 4) * Ts^3 + (n^3 - 9 * n^2 + 24 * n - 14) * Ts * T2s + 4 * (n - 3) * T3s)
  t16 <- 8 * T3 * (Ts^3 + 3 * (n - 3) * Ts * T2s + (n^3 - 9 * n^2 + 26 * n - 22) * T3s)
  t17 <- -16 * (T_val^3 * Us + U * Ts^3) - 6 * (T_val * T2 * Us + U * Ts * T2s) * (2 * n^2 - 10 * n + 16)
  t18 <- -8 * (T3 * Us + U * T3s) * (3 * n^2 - 15 * n + 16) - (T_val^3 * Bs + B * Ts^3) * (6 * n^2 - 30 * n + 24)
  t19 <- -6 * (T_val * T2 * Bs + B * Ts * T2s) * (4 * n^2 - 20 * n + 24) - 8 * (T3 * Bs + B * T3s) * (3 * n^2 - 15 * n + 24)
  t201 <- 24 * (T_val^3 * Rs + R_val * Ts^3) + 6 * (T_val * T2 * Rs + R_val * Ts * T2s) * (2 * n^2 - 10 * n + 24)
  t202 <- 8 * (T3 * Rs + R_val * T3s) * (3 * n^2 - 15 * n + 24) + (3 * n^2 - 15 * n + 6) * (T_val^3 * Ts * S2s + T_val * S2 * Ts^3)
  t203 <- 6 * (T_val * T2 * Ts * S2s + Ts * T2s * T_val * S2) * (n^2 - 5 * n + 6) + 48 * (T3 * Ts * S2s + T3s * T_val * S2)
  t20 <- -(n - 2) * (t201 + t202 + t203)
  
  temp31 <- t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 + t19 + t20
  temp32 <- n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5)
  mom3 <- temp31 / temp32
  
  skewness.rv <- (mom3 - 3 * mean.rv * variance.rv - mean.rv^3) / variance.rv^1.5
  m3 <- as.numeric(skewness.rv)
  
  Z_stat <- (Tstat - mean.rv) / sqrt(variance.rv)
  if (abs(m3) < 1e-6) return(1 - pnorm(Z_stat))
  
  shape <- 4 / m3^2; scale <- m3 / 2; location <- -2 / m3      
  PIIIpars <- list(shape = shape, location = location, scale = scale)
  
  pv <- tryCatch({
    ppearsonIII(Z_stat, params = PIIIpars, lower.tail = FALSE)
  }, error = function(e) {
    return(NA)
  })
  return(pv)
}

# -------------------------------------------------------------------------
# 3. DATA LOADING & PRECALCULATION
# -------------------------------------------------------------------------
cat("Loading Data for Chr", chr_num, "Chunk", chunk_id, "...\n")
if(chunk_id == 0) cat("chr\tpos\tn.ibd.segs\tp.value\n", file = output_file)

if(file.exists(cor_path)) {
  K <- as.matrix(read.table(cor_path))
  K <- solve(K)
} else {
  stop(paste("Correlation matrix file not found:", cor_path))
}

set.seed(12345)
if(!file.exists(res_path)) stop(paste("Residual matrix file not found:", res_path))
r = as.matrix(fread(res_path))

# Pre-calculate invariant statistics
rtr <- crossprod(r); M <- K %*% rtr; Ts <- sum(K * t(rtr)); T2s <- sum(M * t(M))
rK <- r %*% K; diagW <- rowSums(rK * r); S2s <- sum(diagW^2); S3s <- sum(diagW^3)
v <- crossprod(r, diagW); Bs <- crossprod(v, K %*% v); M2 <- M %*% K
diagW2 <- rowSums((r %*% M2) * r); Rs <- sum(diagW * diagW2)
T3s <- sum(diag(M %*% M %*% M)); W <- r%*% K %*% t(r); Us <- sum(W^3)
stats_precalc <- list(Ts=Ts, T2s=T2s, S2s=S2s, S3s=S3s, Us=Us, Rs=Rs, Bs=Bs, T3s=T3s)
rm(W); gc()

sites <- fread(file.path(ibd_dir, "sites.txt"))
samples_vcf <- fread(file.path(ibd_dir, "samples.txt"), header=F)$V1
N_vcf <- length(samples_vcf)

# Handle Subject ID Matching flexibly
if(id_path != "NONE" && file.exists(id_path)) {
  id_include <- fread(id_path, header=F)$V1
} else {
  cat("[WARNING] No subject ID file provided. Assuming rows in the residual matrix correspond exactly to samples.txt.\n")
  id_include <- samples_vcf
}

n <- length(id_include)
if(nrow(r) != n) stop(paste("Mismatch: Residual rows (", nrow(r), ") != Subject IDs (", n, ")"))

match_idx <- match(id_include, samples_vcf)
rows <- rep(1:n, each=2)
cols <- c(rbind(2*match_idx - 1, 2*match_idx))
P_proj <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, N_vcf * 2))

# -------------------------------------------------------------------------
# 4. SMART CHECKPOINT INITIALIZATION & FAST-FORWARD
# -------------------------------------------------------------------------
target_state <- chunk_start - 1
cat("Locating optimal IBD checkpoint for Target State:", target_state, "...\n")

all_mtx_files <- list.files(ibd_dir, pattern = "^ibd_mat_[0-9]+\\.mtx$")
if (length(all_mtx_files) == 0) stop("No checkpoint files found in", ibd_dir)

checkpoint_indices <- as.numeric(gsub("ibd_mat_|\\.mtx", "", all_mtx_files))
valid_checkpoints <- checkpoint_indices[checkpoint_indices <= target_state]

if (length(valid_checkpoints) == 0) valid_checkpoints <- c(0)
best_checkpoint <- max(valid_checkpoints)

cat("  -> Loading Checkpoint:", best_checkpoint, "\n")

mat_file <- file.path(ibd_dir, paste0("ibd_mat_", best_checkpoint, ".mtx"))
K_hap_start <- readMM(mat_file)
X_sparse <- P_proj %*% K_hap_start %*% t(P_proj)
diag(X_sparse) <- diag(X_sparse) / 2
rm(K_hap_start); gc()

# Fast-Forward Gap (Optimized: Native sparse matrix addition only)
if (best_checkpoint < target_state) {
  cat("  -> Bridging gap: Applying diffs from", best_checkpoint + 1, "to", target_state, "\n")
  
  for (j in (best_checkpoint + 1):target_state) {
    diff_file <- file.path(ibd_dir, paste0("delta_", j, ".diff"))
    if (file.exists(diff_file)) {
      dK_hap <- read_diff_as_dK(diff_file, N_vcf * 2)
      if (!is.null(dK_hap)) {
        dX <- P_proj %*% dK_hap %*% t(P_proj)
        diag(dX) <- diag(dX) / 2
        if (nnzero(dX) > 0) {
          X_sparse <- X_sparse + dX
        }
      }
    }
  }
}

cat("  -> Initializing Exact State at Target...\n")
n_ibd_curr <- round(sum(X_sparse) / 2)
Z <- matrix(sample(c(-1.0, 1.0), n * N_hutchinson, replace = TRUE), nrow = n, ncol = N_hutchinson)
XZ <- X_sparse %*% Z

col_sums <- colSums(X_sparse)
S1 <- sum(X_sparse)
S2_sum <- sum(col_sums^2)
temp_vec <- X_sparse %*% col_sums
S3_sum <- sum(col_sums * temp_vec)

raw_trace_X3 <- hutchinson_trace_X3_initial(X_sparse, Z)

term2 <- (3.0 / n) * S3_sum
term3 <- (3.0 / (as.numeric(n)^2)) * S1 * S2_sum
term4 <- (1.0 / (as.numeric(n)^3)) * S1^3
prev_T3_num <- raw_trace_X3 - term2 + term3 - term4

# -------------------------------------------------------------------------
# 5. CHUNK SCAN & DYNAMIC UPDATES
# -------------------------------------------------------------------------
t1 <- proc.time()
calc_obj <- new(HiFiMAPCalculator, X_sparse, r, K)

# Push the exact starting T3_num into the C++ object so baseline is correct
dX_zero <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)

if (chunk_id == 0) {
  cat("--- Processing Position 0 ---\n")
  stats_pos0 <- calc_obj$update_and_calculate(dX_zero, prev_T3_num)
  pval <- calculate_p_value(stats_pos0, n, stats_precalc)
  cat(chr_num, sites$BP[1], n_ibd_curr, pval, "\n", file = output_file, append = TRUE, sep = "\t")
} else {
  invisible(calc_obj$update_and_calculate(dX_zero, prev_T3_num))
}
rm(dX_zero); gc()

cat("Starting Main Chunk Loop (", chunk_start, "to", chunk_end, ")...\n")

for (i in chunk_start:chunk_end) {
  diff_file <- file.path(ibd_dir, paste0("delta_", i, ".diff"))
  
  if (file.exists(diff_file)) {
    dK_hap <- read_diff_as_dK(diff_file, N_vcf * 2)
    
    if (!is.null(dK_hap)) {
      dX <- P_proj %*% dK_hap %*% t(P_proj)
      diag(dX) <- diag(dX) / 2
      
      if (nnzero(dX) > 0) {
        
        v_1tdX <- colSums(dX)
        vec_X_vdX1 <- calc_obj$multiply_X_vec(v_1tdX)
        
        # 1. Hutchinson Dynamic Update
        dXZ <- dX %*% Z
        delta_est <- estimate_dynamic_traces_hutchinson(X_sparse, XZ, dX, dXZ)
        raw_trace_X3 <- raw_trace_X3 + delta_est$tr_X2dX + delta_est$tr_XdX2 + delta_est$tr_dX3
        
        # 2. Exact S3 Term Expansion Updates
        t_S3_1 <- as.numeric(crossprod(col_sums, dX %*% col_sums))
        t_S3_2 <- sum(vec_X_vdX1 * col_sums)
        t_S3_4 <- as.numeric(crossprod(v_1tdX, dX %*% col_sums))
        t_S3_5 <- as.numeric(crossprod(col_sums, dX %*% v_1tdX))
        t_S3_6 <- sum(vec_X_vdX1 * v_1tdX)
        t_S3_7 <- as.numeric(crossprod(v_1tdX, dX %*% v_1tdX))
        
        S3_sum <- S3_sum + t_S3_1 + 2*t_S3_2 + t_S3_4 + t_S3_5 + t_S3_6 + t_S3_7
        S1 <- S1 + sum(dX)
        S2_sum <- S2_sum + 2 * sum(col_sums * v_1tdX) + sum(v_1tdX^2)
        
        # 3. Calculate Global State
        term2 <- (3.0 / n) * S3_sum
        term3 <- (3.0 / (as.numeric(n)^2)) * S1 * S2_sum
        term4 <- (1.0 / (as.numeric(n)^3)) * S1^3
        T3_num <- raw_trace_X3 - term2 + term3 - term4
        
        dT3_num <- T3_num - prev_T3_num
        stats <- calc_obj$update_and_calculate(dX, dT3_num)
        
        prev_T3_num <- T3_num
        
        pval <- calculate_p_value(stats, n, stats_precalc)
        
        # 4. Update the actual matrices for the next loop
        n_ibd_curr <- n_ibd_curr + round(sum(dX) / 2)
        X_sparse <- X_sparse + dX
        XZ <- XZ + dXZ
        col_sums <- col_sums + v_1tdX
        
        cat(chr_num, sites$BP[i+1], n_ibd_curr, pval, "\n", file = output_file, append = TRUE, sep = "\t")
      }
    }
  }
}

t2 <- proc.time()
write.table(t(c(t2 - t1)), file = time_file, quote = FALSE, row.names = FALSE, col.names = TRUE)