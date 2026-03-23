## ============================================================================
##  source("01_generate_data.R")
##
##  Generates the synthetic MET (multi-environment trial) dataset used by
##  03_greta_model_showcase.R and 04_test_translator.R.
##
##  Run ONCE to create met_data.rds in the working directory:
##    setwd("path/to/asreml_to_greta_translator/")
##    source("01_generate_data.R")
##
##  Dataset: plant breeding multi-environment trial (MET)
##    N = 120 observations: 10 genotypes x 4 environments x 3 replicates
##    Includes spatial layout, temporal, genomic, and multi-trait columns.
##    All columns retained to support all 20 model categories.
##
##  Output:
##    met_data.rds  -- list(dat, M_geno, G_mat, A_mat,
##                         geno_levels, env_levels, n_geno, n_env, n_rep, n_marker)
## ============================================================================

set.seed(2024)

# ── Dimensions ────────────────────────────────────────────────────────────────
n_geno   <- 10   # genotypes (varieties)
n_env    <- 4    # environments (sites)
n_rep    <- 3    # replicates per env
n_row    <- 6    # field rows per env-rep
n_col    <- 5    # field cols per env-rep
n_time   <- 4    # measurement times
n_block  <- 3    # incomplete blocks per rep
n_marker <- 50   # SNP markers for genomic examples

obs_per_env <- n_geno * n_rep   # 30 obs per env
N <- n_env * obs_per_env        # 120 base obs

# ── Factor levels ─────────────────────────────────────────────────────────────
geno_levels <- paste0("G", sprintf("%02d", 1:n_geno))
env_levels  <- paste0("E", 1:n_env)
rep_levels  <- paste0("R", 1:n_rep)

# ── Base skeleton ─────────────────────────────────────────────────────────────
df_list <- vector("list", n_env)

for (e in seq_len(n_env)) {
  env_df <- expand.grid(
    geno = geno_levels,
    rep  = rep_levels,
    stringsAsFactors = FALSE
  )
  env_df$env   <- env_levels[e]
  env_df$trial <- e
  df_list[[e]] <- env_df
}

dat <- do.call(rbind, df_list)
dat <- dat[order(dat$env, dat$rep, dat$geno), ]
rownames(dat) <- NULL
N <- nrow(dat)  # 120

# ── Spatial layout: assign row/col within each env × rep ─────────────────────
n_plots <- n_geno  # 10 plots per rep-env block (2 rows × 5 cols)
n_r <- 2
n_c <- 5
layout <- expand.grid(row = 1:n_r, col = 1:n_c)

dat$row  <- rep(layout$row, times = n_env * n_rep)
dat$col  <- rep(layout$col, times = n_env * n_rep)

# Global field row/col (for spatial models over the whole trial)
dat$field_row <- dat$row + (as.integer(factor(dat$rep))  - 1) * n_r
dat$field_col <- dat$col + (as.integer(factor(dat$env))  - 1) * n_c

# ── Block assignment (incomplete blocks within rep) ───────────────────────────
# 10 plots / 3 blocks → blocks of size ~3-4
dat$block <- paste0("B", ceiling((match(dat$geno, geno_levels)) / 4))
dat$block_id <- as.integer(factor(paste(dat$env, dat$rep, dat$block)))

# ── True variance components ──────────────────────────────────────────────────
sigma2_g   <- 1.5       # genotypic variance
sigma2_env <- 2.0       # environment variance
sigma2_rep <- 0.3       # rep-in-env
sigma2_blk <- 0.4       # block-in-rep
sigma2_gxe <- 0.8       # GxE interaction
sigma2_e   <- 0.25      # residual

# Per-env residual SDs (heterogeneous)
env_resid_sd <- c(0.3, 0.6, 0.4, 0.5)

# ── Random effects ────────────────────────────────────────────────────────────
geno_eff  <- rnorm(n_geno,  0, sqrt(sigma2_g))
env_eff   <- rnorm(n_env,   0, sqrt(sigma2_env))
rep_eff   <- rnorm(n_env * n_rep, 0, sqrt(sigma2_rep))
block_eff <- rnorm(max(dat$block_id), 0, sqrt(sigma2_blk))

# GxE interaction matrix
GxE_mat <- matrix(rnorm(n_geno * n_env, 0, sqrt(sigma2_gxe)),
                  nrow = n_geno, ncol = n_env)
colnames(GxE_mat) <- env_levels
rownames(GxE_mat) <- geno_levels

# ── Spatial AR1×AR1 surface ───────────────────────────────────────────────────
rho_r <- 0.7
rho_c <- 0.5
make_ar1_cov <- function(n, rho) {
  outer(1:n, 1:n, function(i, j) rho^abs(i - j))
}
R_row <- make_ar1_cov(max(dat$field_row), rho_r)
R_col <- make_ar1_cov(max(dat$field_col), rho_c)
R_spatial <- kronecker(R_row, R_col) * 0.5
L_sp <- t(chol(R_spatial + diag(nrow(R_spatial)) * 1e-6))
sp_noise_full <- as.vector(L_sp %*% rnorm(nrow(R_spatial)))
dat$spatial_effect <- sp_noise_full[
  (dat$field_row - 1) * max(dat$field_col) + dat$field_col
]

# ── Assemble Yield (continuous, Gaussian) ────────────────────────────────────
grand_mean  <- 8.0
geno_idx    <- match(dat$geno, geno_levels)
env_idx     <- match(dat$env, env_levels)
rep_key     <- (env_idx - 1) * n_rep + match(dat$rep, rep_levels)
resid_sd    <- env_resid_sd[env_idx]

dat$mu_yield <-
  grand_mean +
  geno_eff[geno_idx] +
  env_eff[env_idx] +
  rep_eff[rep_key] +
  block_eff[dat$block_id] +
  GxE_mat[cbind(geno_idx, env_idx)] +
  dat$spatial_effect

dat$yield <- dat$mu_yield + rnorm(N, 0, resid_sd)
dat$yield <- round(dat$yield, 3)

# ── Second trait (for bivariate / multivariate models) ───────────────────────
# Correlated with yield via a shared genetic component
rho_traits <- 0.7
geno_eff2  <- rho_traits * geno_eff +
              sqrt(1 - rho_traits^2) * rnorm(n_geno, 0, sqrt(sigma2_g))
dat$trait2 <- grand_mean * 0.5 +
              geno_eff2[geno_idx] +
              rnorm(N, 0, 0.5)
dat$trait2 <- round(dat$trait2, 3)

# ── Covariate x (continuous, for regression / random regression) ──────────────
dat$x_cov  <- round(rnorm(N, mean = 5, sd = 1.5), 2)

# ── Time variable (for longitudinal / random regression) ─────────────────────
# Each observation "measured" at one of 4 time points
dat$time   <- rep(1:n_time, length.out = N)
dat$time_c <- dat$time - mean(dat$time)  # centred

# Repeat-measures version: expand to long longitudinal dataset embedded as cols
long_rows <- rep(seq_len(N), each = n_time)
dat_long_y <- sapply(1:n_time, function(t) {
  dat$mu_yield + t * 0.3 * geno_eff[geno_idx] + rnorm(N, 0, 0.4)
})
# Store first 4 columns as y_t1 ... y_t4
for (t in 1:n_time) {
  dat[[paste0("y_t", t)]] <- round(dat_long_y[, t], 3)
}

# ── Count response (Poisson) ──────────────────────────────────────────────────
lambda_count  <- exp(0.5 + 0.3 * geno_eff[geno_idx] + 0.2 * env_eff[env_idx])
dat$count_y   <- rpois(N, lambda = lambda_count)

# ── Binary response (Bernoulli) ───────────────────────────────────────────────
logit_p      <- -1 + 0.5 * geno_eff[geno_idx] + 0.3 * env_eff[env_idx]
dat$binary_y <- rbinom(N, 1, prob = plogis(logit_p))

# ── Proportion / rate response (Beta regression) ─────────────────────────────
mu_beta      <- plogis(logit_p)
phi_beta     <- 8
dat$prop_y   <- round(rbeta(N,
                        shape1 = mu_beta * phi_beta,
                        shape2 = (1 - mu_beta) * phi_beta), 4)
# Ensure strictly in (0,1)
dat$prop_y   <- pmin(pmax(dat$prop_y, 1e-4), 1 - 1e-4)

# ── Overdispersed count (Negative Binomial) ───────────────────────────────────
dat$nb_y <- rnbinom(N, size = 2, mu = lambda_count)

# ── Ordinal response (1–5 disease score) ─────────────────────────────────────
thresholds   <- c(-1.5, -0.5, 0.5, 1.5)
latent_ord   <- logit_p + rnorm(N, 0, 0.5)
dat$ordinal_y <- cut(latent_ord,
                     breaks = c(-Inf, thresholds, Inf),
                     labels = 1:5, ordered_result = TRUE)
dat$ordinal_y <- as.integer(dat$ordinal_y)

# ── Continuous positive (Gamma) ───────────────────────────────────────────────
mu_gam       <- exp(1.5 + 0.3 * geno_eff[geno_idx])
shape_gam    <- 3
dat$gamma_y  <- round(rgamma(N, shape = shape_gam,
                              rate  = shape_gam / mu_gam), 3)

# ── Survival / time-to-event ──────────────────────────────────────────────────
haz_rate     <- exp(-2 + 0.3 * geno_eff[geno_idx])
dat$surv_time   <- round(rexp(N, rate = haz_rate), 2)
dat$surv_status <- rbinom(N, 1, 0.8)   # 1 = event, 0 = censored

# ── Weights (for weighted model) ─────────────────────────────────────────────
dat$wt       <- round(runif(N, 0.5, 2.0), 3)

# ── Genomic marker matrix M (n_obs × n_marker) ───────────────────────────────
# Store as separate object; annotated to dat via geno_id
M_geno <- matrix(rbinom(n_geno * n_marker, 2, 0.3),
                 nrow = n_geno, ncol = n_marker,
                 dimnames = list(geno_levels,
                                 paste0("SNP", sprintf("%03d", 1:n_marker))))

# Genomic relationship matrix G = MM'/p (VanRaden 2008)
p_freq <- colMeans(M_geno) / 2
M_c    <- sweep(M_geno, 2, 2 * p_freq, "-")
G_mat  <- (M_c %*% t(M_c)) / (2 * sum(p_freq * (1 - p_freq)))
G_mat  <- G_mat + diag(n_geno) * 0.01   # regularise
dimnames(G_mat) <- list(geno_levels, geno_levels)

# Pedigree-based A (simulated tridiagonal for illustration)
A_mat <- diag(n_geno)
for (i in 2:n_geno) {
  A_mat[i, i-1] <- A_mat[i-1, i] <- 0.25
}
dimnames(A_mat) <- list(geno_levels, geno_levels)

# ── Integer indices (required for greta indexing) ─────────────────────────────
dat$geno_id  <- as.integer(factor(dat$geno,  levels = geno_levels))
dat$env_id   <- as.integer(factor(dat$env,   levels = env_levels))
dat$rep_id   <- as.integer(factor(dat$rep,   levels = rep_levels))
dat$rep_env_id <- as.integer(factor(paste(dat$env, dat$rep)))
# block nested in rep×env
dat$block_nested_id <- as.integer(factor(paste(dat$env, dat$rep, dat$block)))

# ── Factor versions ───────────────────────────────────────────────────────────
dat$geno_f   <- factor(dat$geno,  levels = geno_levels)
dat$env_f    <- factor(dat$env,   levels = env_levels)
dat$rep_f    <- factor(dat$rep,   levels = rep_levels)
dat$block_f  <- factor(paste(dat$env, dat$rep, dat$block))

# ── Marker matrix row index (for GWAS / rrBLUP) ──────────────────────────────
dat$M_row   <- dat$geno_id   # dat$M_row[i] → row of M_geno

# ── Summary ───────────────────────────────────────────────────────────────────
cat("Dataset dimensions:", nrow(dat), "rows ×", ncol(dat), "columns\n")
cat("Genotypes:", n_geno, "| Environments:", n_env,
    "| Reps:", n_rep, "| N:", N, "\n")
cat("Columns:\n")
print(names(dat))
cat("\nFirst 6 rows:\n")
print(head(dat))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(list(dat = dat, M_geno = M_geno, G_mat = G_mat, A_mat = A_mat,
             geno_levels = geno_levels, env_levels = env_levels,
             n_geno = n_geno, n_env = n_env, n_rep = n_rep,
             n_marker = n_marker),
        file = "met_data.rds")
cat("\nSaved: met_data.rds\n")
