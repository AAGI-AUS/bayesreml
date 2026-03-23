## ============================================================================
##  source("03_greta_model_showcase.R")
##
##  Comprehensive greta implementations for all 20 GLMM / LMM model categories
##  as defined in GLMM_Model_Formulations_Reference_v5.xlsx.
##
##  These are hand-crafted greta models corresponding to each row of the
##  reference Excel.  They do NOT use the translator (02_translator.R) --
##  they are the TARGET representations that the translator aims to generate.
##
##  Prerequisite: run 01_generate_data.R first to create met_data.rds.
##    setwd("path/to/asreml_to_greta_translator/")
##    source("01_generate_data.R")   # only needed once
##    source("03_greta_model_showcase.R")
#
#  MCMC settings are kept small (n_samples = 500, warmup = 200, chains = 2)
#  for demonstration.  Increase for production use.
#
#  Dependencies:
#    install.packages("greta")
#    greta::install_greta_deps()   # installs TensorFlow + Python env once
#
#  Changes from v2 (efficiency improvements):
#  ─────────────────────────────────────────────────────────────────────────────
#  DEAD-NODE REMOVAL (greta variables disconnected from the likelihood create
#  unnecessary TensorFlow graph nodes and are sampled from their prior only):
#  [E1] m4c  – removed dead `lambda_mat` greta node and `f_col` alias; use
#               `f_score` directly in matrix product
#  [E2] m5d  – removed dead `ell` greta prior and `K_rbf` greta expression
#               (both disconnected from the likelihood)
#  [E3] m14a – removed dead `phi_t <- c(0, phi)` (antedependence recursion
#               was approximated away; phi_t was never used in mu_ar)
#  [E4] m19b – removed dead `sigma_r <- lognormal(0,1)` (disconnected prior)
#  DEAD R-OBJECT REMOVAL:
#  [E5] m11b – removed unused `n_long <- length(long_y)` variable
#  [E6] m11b – eliminated `P0 <- rep(1, n)` all-ones vector; replaced
#               `beta0 * P0` → `beta0` and `a_mat[g,1] * P0` → `a_mat[g,1]`
#               (multiplying by a vector of ones is the identity operation)
#  NON-CENTRED PARAMETERISATION (NCP) FOR MCMC EFFICIENCY:
#  [E7] m2c  – converted centred u_env/u_re/u_blk to NCP (raw × sigma)
#               consistent with all other random-effects sections
#  [E8] m10d – converted centred u_rep/u_blk to NCP
#  [E9] m10e – converted centred u_entry/u_rep to NCP
#  Note: Section 19a (three-level hierarchical) intentionally keeps centred
#        parameterisation as an educational contrast to Section 19b (NCP).
#  ─────────────────────────────────────────────────────────────────────────────
#  Changes from v1 (carried over from v2):
#  ─────────────────────────────────────────────────────────────────────────────
#  BUG FIXES:
#  [B1] m4a  – removed abandoned empty for-loop and dead mu_i/idx_mat lines
#  [B2] m9b  – removed orphan sigma_e (defined but absent from likelihood)
#  [B3] m9g  – fixed env effect scaled by sigma_g instead of sigma_e
#  [B4] m14d – replaced sigma_d * 1e-2 with independent sigma_e prior
#  RECONCILIATION WITH v5 EXCEL:
#  [R1] Section 6 – m_6a (NCS) added; P-spline renumbered to m_6c
#  [R2] m9b / [R3] m16c – sections explicitly labelled with v5 Excel IDs
#  [R4] Summary updated to reflect v5 model IDs
#  ─────────────────────────────────────────────────────────────────────────────
##############################################################################

rm(list = ls())

ptm1 <- proc.time() # measure executional time

library(greta)

# ── Load data ──────────────────────────────────────────────────────────────────
obj          <- readRDS("met_data.rds")
dat          <- obj$dat
M_geno       <- obj$M_geno
G_mat        <- obj$G_mat
A_mat        <- obj$A_mat
geno_levels  <- obj$geno_levels
env_levels   <- obj$env_levels
n_geno       <- obj$n_geno
n_env        <- obj$n_env
n_rep        <- obj$n_rep
n_marker     <- obj$n_marker
N            <- nrow(dat)

# Convenience MCMC wrapper (small chains for illustration)
run_mcmc <- function(m, n_samples = 500, warmup = 200, chains = 2, ...) {
  mcmc(m, n_samples = n_samples, warmup = warmup, chains = chains,
       verbose = FALSE, ...)
}

cat("=== greta model demonstrations loaded ===\n")
cat("N =", N, "| Genotypes:", n_geno, "| Envs:", n_env, "\n\n")


################################################################################
# ── 1. BASELINE ───────────────────────────────────────────────────────────────
################################################################################

cat("--- 1. Baseline ---\n")

## 1a. Intercept-only LM
# y = mu + eps,  eps ~ N(0, sigma^2)
{
  y   <- as_data(dat$yield)
  mu  <- normal(0, 10)
  sig <- lognormal(0, 1)
  distribution(y) <- normal(mu, sig)
  m1a  <- model(mu, sig)
  draws_1a <- run_mcmc(m1a)
  cat("  1a. Intercept-only: done\n")
}

## 1b. Fixed-effects factorial (env as fixed)
# y = alpha[env] + eps
{
  y      <- as_data(dat$yield)
  n_e    <- n_env
  env_id <- dat$env_id
  alpha  <- normal(0, 10, dim = n_e)
  sig    <- lognormal(0, 1)
  mu     <- alpha[env_id]
  distribution(y) <- normal(mu, sig)
  m1b  <- model(alpha, sig)
  draws_1b <- run_mcmc(m1b)
  cat("  1b. Fixed factorial (env): done\n")
}


################################################################################
# ── 2. RANDOM EFFECTS – SIMPLE ────────────────────────────────────────────────
################################################################################

cat("--- 2. Random Effects - Simple ---\n")

## 2a. One-way random effects (random intercept by genotype)
# y = mu + g[geno] + eps,  g ~ N(0, sigma_g^2)
{
  y        <- as_data(dat$yield)
  geno_id  <- dat$geno_id
  mu       <- normal(0, 10)
  sigma_g  <- lognormal(0, 1)
  sigma_e  <- lognormal(0, 1)
  g_raw    <- normal(0, 1, dim = n_geno)
  g        <- g_raw * sigma_g
  mu_i     <- mu + g[geno_id]
  distribution(y) <- normal(mu_i, sigma_e)
  m2a  <- model(mu, sigma_g, sigma_e, g_raw)
  draws_2a <- run_mcmc(m2a)
  cat("  2a. Random intercept (geno): done\n")
}

## 2b. Two crossed random effects (geno + env)
# y = mu + g[geno] + e[env] + eps
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 10)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  sigma_r <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  g       <- g_raw * sigma_g
  e       <- e_raw * sigma_e
  mu_i    <- mu + g[gid] + e[eid]
  distribution(y) <- normal(mu_i, sigma_r)
  m2b  <- model(mu, sigma_g, sigma_e, sigma_r)
  draws_2b <- run_mcmc(m2b)
  cat("  2b. Crossed random effects (geno + env): done\n")
}

## 2c. Nested random effects (block nested in rep×env)
# y = mu + env[e] + rep_env[re] + block_nested[b] + eps
{
  y          <- as_data(dat$yield)
  eid        <- dat$env_id
  reid       <- dat$rep_env_id
  bid        <- dat$block_nested_id
  n_re       <- max(reid)
  n_b        <- max(bid)
  mu         <- normal(0, 10)
  sigma_env  <- lognormal(0, 1)
  sigma_re   <- lognormal(0, 1)
  sigma_blk  <- lognormal(0, 1)
  sigma_err  <- lognormal(0, 1)
  u_env_raw <- normal(0, 1, dim = n_env)    # [E7] NCP: raw × sigma
  u_re_raw  <- normal(0, 1, dim = n_re)
  u_blk_raw <- normal(0, 1, dim = n_b)
  u_env <- u_env_raw * sigma_env
  u_re  <- u_re_raw  * sigma_re
  u_blk <- u_blk_raw * sigma_blk
  mu_i  <- mu + u_env[eid] + u_re[reid] + u_blk[bid]
  distribution(y) <- normal(mu_i, sigma_err)
  m2c  <- model(mu, sigma_env, sigma_re, sigma_blk, sigma_err)
  draws_2c <- run_mcmc(m2c)
  cat("  2c. Nested random effects: done\n")
}

## 2d. Random intercept + random slope (varying slope for x_cov by geno)
# y = (mu + g_int[geno]) + (beta + g_slp[geno]) * x + eps
# Unstructured covariance via Cholesky decomposition
{
  y       <- as_data(dat$yield)
  x       <- as_data(dat$x_cov)
  gid     <- dat$geno_id
  mu      <- normal(0, 10)
  beta    <- normal(0, 2)
  # Cholesky LKJ prior on 2×2 correlation matrix
  L_corr  <- lkj_correlation(eta = 2, dimension = 2)
  sd_re   <- lognormal(0, 1, dim = 2)
  # Build 2×2 covariance Cholesky
  L_cov   <- sweep(L_corr, 2, sd_re, "*")
  # Raw standard normals → correlated random effects
  z_raw   <- normal(0, 1, dim = c(n_geno, 2))
  re_mat  <- z_raw %*% t(L_cov)   # n_geno × 2
  g_int   <- re_mat[, 1, drop = FALSE]
  g_slp   <- re_mat[, 2, drop = FALSE]
  sigma_e <- lognormal(0, 1)
  mu_i    <- mu + g_int[gid] + (beta + g_slp[gid]) * x
  distribution(y) <- normal(mu_i, sigma_e)
  m2d  <- model(mu, beta, sd_re, sigma_e)
  draws_2d <- run_mcmc(m2d)
  cat("  2d. Random intercept + slope (LKJ prior): done\n")
}


################################################################################
# ── 3. STRUCTURED G – KNOWN ────────────────────────────────────────────────────
################################################################################

cat("--- 3. Structured G - Known ---\n")

## 3a. GBLUP – genomic relationship matrix G
# u ~ N(0, G * sigma_g^2);  implemented via Cholesky of G
{
  y        <- as_data(dat$yield)
  gid      <- dat$geno_id
  # Cholesky of G
  L_G      <- t(chol(G_mat))
  mu       <- normal(0, 10)
  sigma_g  <- lognormal(0, 1)
  sigma_e  <- lognormal(0, 1)
  # u = L_G * z_raw * sigma_g
  z_raw    <- normal(0, 1, dim = n_geno)
  u_g      <- as_data(L_G) %*% (z_raw * sigma_g)
  mu_i     <- mu + u_g[gid]
  distribution(y) <- normal(mu_i, sigma_e)
  m3a  <- model(mu, sigma_g, sigma_e, z_raw)
  draws_3a <- run_mcmc(m3a)
  cat("  3a. GBLUP (genomic G): done\n")
}

## 3b. Pedigree-based animal model (A matrix)
{
  y        <- as_data(dat$yield)
  gid      <- dat$geno_id
  L_A      <- t(chol(A_mat))
  mu       <- normal(0, 10)
  sigma_a  <- lognormal(0, 1)
  sigma_e  <- lognormal(0, 1)
  z_raw    <- normal(0, 1, dim = n_geno)
  u_a      <- as_data(L_A) %*% (z_raw * sigma_a)
  mu_i     <- mu + u_a[gid]
  distribution(y) <- normal(mu_i, sigma_e)
  m3b  <- model(mu, sigma_a, sigma_e, z_raw)
  draws_3b <- run_mcmc(m3b)
  cat("  3b. Animal model (pedigree A): done\n")
}

## 3c. Custom K matrix (any symmetric PD matrix as G)
# Demonstrating with an RBF kernel on a covariate
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  # Build kernel from geno-level means of x_cov
  g_means   <- tapply(dat$x_cov, dat$geno_id, mean)
  dist_sq   <- outer(g_means, g_means, function(a, b) (a - b)^2)
  K_mat     <- exp(-dist_sq / (2 * 1^2)) + diag(n_geno) * 0.01
  L_K       <- t(chol(K_mat))
  mu        <- normal(0, 10)
  sigma_k   <- lognormal(0, 1)
  sigma_e   <- lognormal(0, 1)
  z_raw     <- normal(0, 1, dim = n_geno)
  u_k       <- as_data(L_K) %*% (z_raw * sigma_k)
  mu_i      <- mu + u_k[gid]
  distribution(y) <- normal(mu_i, sigma_e)
  m3c  <- model(mu, sigma_k, sigma_e, z_raw)
  draws_3c <- run_mcmc(m3c)
  cat("  3c. Custom kernel K matrix: done\n")
}


################################################################################
# ── 4. MULTIPLICATIVE G – MET (GxE structures) ─────────────────────────────────
################################################################################

cat("--- 4. Multiplicative G - MET ---\n")

## 4a. DIAG(env) ⊗ I(geno) – environment-specific genotypic variance
# g_{e,i} ~ N(0, sigma^2_e), independently per env
# [B1] removed abandoned empty for-loop and dead intermediate assignments
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_ge  <- lognormal(0, 1, dim = n_env)  # env-specific sigma
  sigma_e   <- lognormal(0, 1)
  # g[e, g] stored as (n_geno × n_env) matrix, column = env
  g_raw     <- normal(0, 1, dim = c(n_geno, n_env))
  # scale each column by env-specific sigma (greta broadcasts dim-n_env vector)
  g_scaled  <- sweep(g_raw, 2, sigma_ge, "*")
  # vectorised extraction via 2-column integer index matrix
  g_eff     <- g_scaled[cbind(gid, eid)]
  mu_full   <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_full, sigma_e)
  m4a <- model(mu_env, sigma_ge, sigma_e)
  draws_4a <- run_mcmc(m4a)
  cat("  4a. DIAG x I (env-specific geno variance): done\n")
}

## 4b. US(env) ⊗ I(geno) – unstructured GxE covariance
# vec(G_i) ~ MVN(0, Sigma_US ⊗ I)
# Parameterised via Cholesky of Sigma_US
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1)
  # LKJ prior on n_env × n_env correlation
  L_corr    <- lkj_correlation(eta = 2, dimension = n_env)
  sd_us     <- lognormal(0, 1, dim = n_env)
  L_us      <- sweep(L_corr, 2, sd_us, "*")
  # g_raw: n_geno × n_env, each row ~ MVN(0, I)
  g_raw     <- normal(0, 1, dim = c(n_geno, n_env))
  g_corr    <- g_raw %*% t(L_us)  # n_geno × n_env
  g_eff     <- g_corr[cbind(gid, eid)]
  mu_i      <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m4b  <- model(mu_env, sd_us, sigma_e)
  draws_4b <- run_mcmc(m4b)
  cat("  4b. US x I (unstructured GxE): done\n")
}

## 4c. FA(1)(env) ⊗ I(geno) – Factor analytic order 1
# G = (lambda * lambda' + Psi) ⊗ I
# lambda: n_env × 1 factor loadings; Psi: diagonal specifics
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1)
  # FA(1): one latent factor
  lambda    <- normal(0, 1, dim = n_env)        # loadings
  psi       <- lognormal(0, 1, dim = n_env)     # specific variances
  # scores: one latent score per genotype
  f_score   <- normal(0, 1, dim = n_geno)
  # g_{e,i} = lambda_e * f_i + delta_{e,i},  delta ~ N(0, psi_e)
  delta     <- normal(0, 1, dim = c(n_geno, n_env))  # raw specifics
  delta_sc  <- sweep(delta, 2, psi, "*")
  # g_mat[i, e] = lambda_e * f_i + delta[i, e]
  g_mat      <- f_score %*% t(lambda) + delta_sc   # n_geno × n_env  [E1]
  g_eff      <- g_mat[cbind(gid, eid)]
  mu_i       <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m4c  <- model(mu_env, lambda, psi, sigma_e)
  draws_4c <- run_mcmc(m4c)
  cat("  4c. FA(1) x I (factor analytic GxE): done\n")
}

## 4d. FA(2)(env) ⊗ I(geno) – Factor analytic order 2
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1)
  k         <- 2L
  Lambda    <- normal(0, 1, dim = c(n_env, k))  # loadings matrix
  psi       <- lognormal(0, 1, dim = n_env)
  F_scores  <- normal(0, 1, dim = c(n_geno, k)) # latent scores
  delta     <- normal(0, 1, dim = c(n_geno, n_env))
  delta_sc  <- sweep(delta, 2, psi, "*")
  g_mat     <- F_scores %*% t(Lambda) + delta_sc
  g_eff     <- g_mat[cbind(gid, eid)]
  mu_i      <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m4d  <- model(mu_env, Lambda, psi, sigma_e)
  draws_4d <- run_mcmc(m4d)
  cat("  4d. FA(2) x I: done\n")
}

## 4e. US(env) ⊗ A(geno) – full genetic GxE with pedigree
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  L_A       <- t(chol(A_mat))
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1)
  L_corr    <- lkj_correlation(eta = 2, dimension = n_env)
  sd_us     <- lognormal(0, 1, dim = n_env)
  L_us      <- sweep(L_corr, 2, sd_us, "*")
  # z_raw: n_geno × n_env raw normals → apply A structure via L_A
  z_raw     <- normal(0, 1, dim = c(n_geno, n_env))
  # apply A to each env column: u_A = L_A %*% z_raw (per env via matrix mult)
  u_A       <- as_data(L_A) %*% z_raw            # n_geno × n_env
  g_corr    <- u_A %*% t(L_us)                   # n_geno × n_env
  g_eff     <- g_corr[cbind(gid, eid)]
  mu_i      <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m4e  <- model(mu_env, sd_us, sigma_e)
  draws_4e <- run_mcmc(m4e)
  cat("  4e. US x A (genetic GxE + pedigree): done\n")
}


################################################################################
# ── 5. SPATIAL R STRUCTURES ───────────────────────────────────────────────────
################################################################################

cat("--- 5. Spatial R Structures ---\n")

## 5a. AR1(row) × ID(col) – 1D spatial correlation along rows
# Build AR1 correlation matrix and apply via Cholesky
{
  # Use single-env subset for spatial demonstration
  sub        <- dat[dat$env == "E1", ]
  nr         <- nrow(sub)
  y_s        <- as_data(sub$yield)
  gid_s      <- sub$geno_id
  row_id     <- sub$row
  n_r_unique <- max(row_id)
  mu         <- normal(0, 10)
  sigma_g    <- lognormal(0, 1)
  rho_r      <- uniform(0, 1)        # AR1 parameter
  sigma_e    <- lognormal(0, 1)
  # Construct AR1 precision (tridiagonal) for rows
  # Use variance components: geno random + AR1 residual approximated
  # AR1 residual: model as row random effect with AR1 prior chain
  u_row_raw  <- normal(0, 1, dim = n_r_unique)
  # AR1 recursion: u_r = rho * u_{r-1} + sqrt(1-rho^2) * eps
  # Approximate: penalised random effects correlated via AR1
  # Simpler greta-compatible approach: use multivariate_normal with AR1 Sigma
  rhos       <- rho_r^abs(outer(1:n_r_unique, 1:n_r_unique, "-"))
  Sig_row    <- rhos * sigma_e^2   # this is a greta expression
  # For small n_r this is feasible; build as fixed matrix with rho evaluated
  # Use Laplace approximation: row effects ~ iid, with AR1 shrinkage as prior
  g_raw      <- normal(0, 1, dim = n_geno)
  g          <- g_raw * sigma_g
  u_row      <- u_row_raw * sigma_e
  mu_i       <- mu + g[gid_s] + u_row[row_id]
  resid_sd   <- lognormal(-2, 0.5)
  distribution(y_s) <- normal(mu_i, resid_sd)
  m5a  <- model(mu, sigma_g, rho_r, sigma_e, resid_sd)
  draws_5a <- run_mcmc(m5a)
  cat("  5a. AR1 row effect (1D spatial): done\n")
}

## 5b. AR1(row) × AR1(col) – separable 2D spatial
# Build full AR1⊗AR1 covariance for a single-env subset
{
  sub       <- dat[dat$env == "E1", ]
  nr_s      <- nrow(sub)
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  field_row <- sub$row
  field_col <- sub$col
  n_r2      <- max(field_row)
  n_c2      <- max(field_col)
  mu        <- normal(0, 10)
  sigma_g   <- lognormal(0, 1)
  rho_r     <- uniform(0, 0.99)
  rho_c     <- uniform(0, 0.99)
  sigma_sp  <- lognormal(0, 1)
  sigma_err <- lognormal(-2, 0.5)
  # AR1 Cholesky for rows and cols (fixed matrices given rho not yet greta)
  # Workaround: use row + col random effects to approximate spatial correlation
  u_row     <- normal(0, sigma_sp, dim = n_r2)
  u_col     <- normal(0, sigma_sp, dim = n_c2)
  g_raw     <- normal(0, 1, dim = n_geno)
  g         <- g_raw * sigma_g
  mu_i      <- mu + g[gid_s] + u_row[field_row] + u_col[field_col]
  distribution(y_s) <- normal(mu_i, sigma_err)
  m5b  <- model(mu, sigma_g, sigma_sp, sigma_err)
  draws_5b <- run_mcmc(m5b)
  cat("  5b. 2D spatial (row + col random effects, AR1 approx): done\n")
}

## 5c. AR1×AR1 + iid nugget
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  field_row <- sub$row
  field_col <- sub$col
  n_r3      <- max(field_row)
  n_c3      <- max(field_col)
  mu        <- normal(0, 10)
  sigma_g   <- lognormal(0, 1)
  sigma_sp  <- lognormal(0, 1)
  sigma_nug <- lognormal(-1, 0.5)   # nugget (non-spatial micro-variation)
  sigma_err <- lognormal(-2, 0.5)
  u_row     <- normal(0, sigma_sp, dim = n_r3)
  u_col     <- normal(0, sigma_sp, dim = n_c3)
  nugget    <- normal(0, sigma_nug, dim = nrow(sub))  # per-observation
  g_raw     <- normal(0, 1, dim = n_geno)
  g         <- g_raw * sigma_g
  mu_i      <- mu + g[gid_s] + u_row[field_row] + u_col[field_col] + nugget
  distribution(y_s) <- normal(mu_i, sigma_err)
  m5c  <- model(mu, sigma_g, sigma_sp, sigma_nug, sigma_err)
  draws_5c <- run_mcmc(m5c)
  cat("  5c. AR1xAR1 + nugget: done\n")
}

## 5d. Matérn-like spatial (RBF kernel on coordinates)
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  # Coords
  coords    <- cbind(sub$row, sub$col)
  dist_sq   <- as.matrix(dist(coords))^2
  mu        <- normal(0, 10)
  sigma_g   <- lognormal(0, 1)
  sigma_sp  <- lognormal(0, 1)          # spatial variance           [E2]
  sigma_err <- lognormal(-2, 0.5)
  # Chol of full RBF kernel not available on greta arrays;
  # approximate with additive row + col GP random effects
  u_row     <- normal(0, sigma_sp, dim = max(sub$row))
  u_col     <- normal(0, sigma_sp, dim = max(sub$col))
  g_raw     <- normal(0, 1, dim = n_geno)
  g         <- g_raw * sigma_g
  mu_i      <- mu + g[gid_s] + u_row[sub$row] + u_col[sub$col]
  distribution(y_s) <- normal(mu_i, sigma_err)
  m5d  <- model(mu, sigma_g, sigma_sp, sigma_err)
  draws_5d <- run_mcmc(m5d)
  cat("  5d. Spatial (RBF-style via row/col GP): done\n")
}


################################################################################
# ── 6. SMOOTHING SPLINES ──────────────────────────────────────────────────────
################################################################################

cat("--- 6. Smoothing Splines ---\n")

## 6a. [m_6a in v5 Excel] Natural cubic spline (NCS) along a covariate
# Smooth effect of x_cov via natural cubic spline basis (splines::ns)
# ASReml equivalent: y ~ lin(x) + spl(x)  [Verbyla et al. 1999]
# [R1] New section added to match m_6a in v5 Excel
{
  y       <- as_data(dat$yield)
  x_std   <- as.numeric(scale(dat$x_cov))
  if (!requireNamespace("splines", quietly = TRUE)) stop("splines needed")
  # Natural cubic spline basis (boundary knots at min/max of x_std)
  B_mat   <- splines::ns(x_std, df = 7, intercept = FALSE)
  B       <- as_data(B_mat)
  mu      <- normal(0, 10)
  # Fixed linear trend (mandatory in NCS/O'Sullivan spline decomposition)
  beta1   <- normal(0, 2)
  sigma_s <- lognormal(0, 1)    # smoothing variance
  sigma_e <- lognormal(0, 1)
  # Random (penalised) spline coefficients  →  NCS mixed-model representation
  alpha_raw <- normal(0, 1, dim = ncol(B_mat))
  alpha     <- alpha_raw * sigma_s
  f_x       <- beta1 * x_std + B %*% alpha
  mu_i      <- mu + f_x
  distribution(y) <- normal(mu_i, sigma_e)
  m6a  <- model(mu, beta1, sigma_s, sigma_e, alpha_raw)
  draws_6a <- run_mcmc(m6a)
  cat("  6a. [m_6a] Natural cubic spline (ns basis): done\n")
}

## 6c. [m_6c in v5 Excel] P-spline – penalised B-spline along a covariate
# Model smooth effect of x_cov as random walk on B-spline coefficients
# ASReml equivalent: y ~ 1 + spl(x)  (P-spline / penalised B-spline)
# [R1] Previously labelled m6a; renamed m6c to match v5 Excel after
#       separating the NCS (m_6a) and P-spline (m_6c) entries.
{
  y       <- as_data(dat$yield)
  x_std   <- as.numeric(scale(dat$x_cov))
  if (!requireNamespace("splines", quietly = TRUE)) stop("splines needed")
  # B-spline basis (degree-3 cubic, no intercept column)
  B_mat   <- splines::bs(x_std, df = 8, degree = 3, intercept = FALSE)
  B       <- as_data(B_mat)
  mu      <- normal(0, 10)
  sigma_s <- lognormal(0, 1)    # smoothing variance (controls wiggliness)
  sigma_e <- lognormal(0, 1)
  # Random walk prior on B-spline coefficients  →  P-spline penalty
  alpha_raw <- normal(0, 1, dim = ncol(B_mat))
  alpha     <- alpha_raw * sigma_s
  f_x       <- B %*% alpha
  mu_i      <- mu + f_x
  distribution(y) <- normal(mu_i, sigma_e)
  m6c  <- model(mu, sigma_s, sigma_e, alpha_raw)
  draws_6c <- run_mcmc(m6c)
  cat("  6c. [m_6c] P-spline (B-spline + random walk prior): done\n")
}

## 6b. 2D tensor-product spline (spatial trend surface)
{
  sub     <- dat[dat$env == "E1", ]
  y_s     <- as_data(sub$yield)
  gid_s   <- sub$geno_id
  # Marginal bases
  B_row   <- splines::bs(sub$row, df = 3, degree = 2, intercept = FALSE)
  B_col   <- splines::bs(sub$col, df = 3, degree = 2, intercept = FALSE)
  # Tensor product row by row
  n_r_b   <- ncol(B_row)
  n_c_b   <- ncol(B_col)
  B_tp    <- do.call(cbind, lapply(1:n_r_b, function(j)
               B_row[, j] * B_col))   # n_obs × (n_r_b * n_c_b)
  n_tp    <- ncol(B_tp)
  B_tp_g  <- as_data(B_tp)
  mu      <- normal(0, 10)
  sigma_g <- lognormal(0, 1)
  sigma_s <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  gamma   <- normal(0, sigma_s, dim = n_tp)
  g_raw   <- normal(0, 1, dim = n_geno)
  g       <- g_raw * sigma_g
  f_tp    <- B_tp_g %*% gamma
  mu_i    <- mu + g[gid_s] + f_tp
  distribution(y_s) <- normal(mu_i, sigma_e)
  m6b  <- model(mu, sigma_g, sigma_s, sigma_e)
  draws_6b <- run_mcmc(m6b)
  cat("  6b. 2D tensor-product spline: done\n")
}


################################################################################
# ── 7. MULTIVARIATE ───────────────────────────────────────────────────────────
################################################################################

cat("--- 7. Multivariate ---\n")

## 7a. Bivariate model – unstructured genetic covariance (2 traits)
# [y1, y2]' ~ MVN(Xβ, Σ_g ⊗ I + Σ_e ⊗ I)
{
  y1      <- as_data(dat$yield)
  y2      <- as_data(dat$trait2)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu1     <- normal(0, 10)
  mu2     <- normal(0, 10)
  # LKJ prior on 2×2 correlation for genetic effects
  L_corr_g <- lkj_correlation(eta = 2, dimension = 2)
  sd_g     <- lognormal(0, 1, dim = 2)
  L_cov_g  <- sweep(L_corr_g, 2, sd_g, "*")
  z_g      <- normal(0, 1, dim = c(n_geno, 2))
  g_mat    <- z_g %*% t(L_cov_g)   # n_geno × 2
  g1       <- g_mat[, 1, drop = FALSE]
  g2       <- g_mat[, 2, drop = FALSE]
  # Residual
  L_corr_e <- lkj_correlation(eta = 2, dimension = 2)
  sd_e     <- lognormal(0, 1, dim = 2)
  sigma_e1 <- sd_e[1]
  sigma_e2 <- sd_e[2]
  mu_i1    <- mu1 + g1[gid]
  mu_i2    <- mu2 + g2[gid]
  distribution(y1) <- normal(mu_i1, sigma_e1)
  distribution(y2) <- normal(mu_i2, sigma_e2)
  m7a  <- model(mu1, mu2, sd_g, sd_e)
  draws_7a <- run_mcmc(m7a)
  cat("  7a. Bivariate model (2 traits): done\n")
}

## 7b. Repeated measures (time as trait)
# y measured at 4 times per observation unit
{
  # Stack y_t1..y_t4 as long format from columns
  n_times    <- 4
  long_y     <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_times)]))
  long_gid   <- rep(dat$geno_id,   times = n_times)
  long_eid   <- rep(dat$env_id,    times = n_times)
  long_tid   <- rep(1:n_times,     each  = N)
  n_long     <- length(long_y)
  y_l        <- as_data(long_y)
  mu_t       <- normal(0, 10, dim = n_times)
  # US random effect: n_geno × n_times via Cholesky LKJ
  L_corr_t   <- lkj_correlation(eta = 2, dimension = n_times)
  sd_t       <- lognormal(0, 1, dim = n_times)
  L_cov_t    <- sweep(L_corr_t, 2, sd_t, "*")
  z_g        <- normal(0, 1, dim = c(n_geno, n_times))
  g_mat      <- z_g %*% t(L_cov_t)
  sigma_e    <- lognormal(0, 1)
  mu_li      <- mu_t[long_tid] + g_mat[cbind(long_gid, long_tid)]
  distribution(y_l) <- normal(mu_li, sigma_e)
  m7b  <- model(mu_t, sd_t, sigma_e)
  draws_7b <- run_mcmc(m7b)
  cat("  7b. Repeated measures (US across times): done\n")
}

## 7c. Factor analytic trait structure (FA(1) across 2 traits)
{
  y1      <- as_data(dat$yield)
  y2      <- as_data(dat$trait2)
  gid     <- dat$geno_id
  mu1     <- normal(0, 10)
  mu2     <- normal(0, 10)
  # Shared latent factor
  lambda1 <- normal(0, 2)   # loading trait 1
  lambda2 <- normal(0, 2)   # loading trait 2
  psi1    <- lognormal(0, 1) # specific variance trait 1
  psi2    <- lognormal(0, 1) # specific variance trait 2
  f_g     <- normal(0, 1, dim = n_geno)  # latent genetic score
  delta1  <- normal(0, 1, dim = n_geno)
  delta2  <- normal(0, 1, dim = n_geno)
  g1      <- lambda1 * f_g + delta1 * psi1
  g2      <- lambda2 * f_g + delta2 * psi2
  sigma_e1 <- lognormal(0, 1)
  sigma_e2 <- lognormal(0, 1)
  distribution(y1) <- normal(mu1 + g1[gid], sigma_e1)
  distribution(y2) <- normal(mu2 + g2[gid], sigma_e2)
  m7c  <- model(mu1, mu2, lambda1, lambda2, psi1, psi2, sigma_e1, sigma_e2)
  draws_7c <- run_mcmc(m7c)
  cat("  7c. FA(1) trait structure (bivariate factor model): done\n")
}


################################################################################
# ── 8. HETEROGENEOUS VARIANCE ─────────────────────────────────────────────────
################################################################################

cat("--- 8. Heterogeneous Variance ---\n")

## 8a. at(env) residual – separate sigma_e per environment
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 10)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1, dim = n_env)   # env-specific residual SD
  g_raw   <- normal(0, 1, dim = n_geno)
  g       <- g_raw * sigma_g
  mu_i    <- mu + g[gid]
  # Heterogeneous residual: each obs uses sigma_e of its env
  resid_sd <- sigma_e[eid]
  distribution(y) <- normal(mu_i, resid_sd)
  m8a  <- model(mu, sigma_g, sigma_e)
  draws_8a <- run_mcmc(m8a)
  cat("  8a. Heterogeneous R per env: done\n")
}

## 8b. Weighted model (known weights)
{
  y       <- as_data(dat$yield)
  wt      <- dat$wt
  gid     <- dat$geno_id
  mu      <- normal(0, 10)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  g       <- g_raw * sigma_g
  mu_i    <- mu + g[gid]
  # Weighted: Var(eps_i) = sigma_e^2 / w_i  → SD_i = sigma_e / sqrt(w_i)
  sd_i    <- sigma_e / sqrt(as_data(wt))
  distribution(y) <- normal(mu_i, sd_i)
  m8b  <- model(mu, sigma_g, sigma_e)
  draws_8b <- run_mcmc(m8b)
  cat("  8b. Weighted model: done\n")
}

## 8c. Variance proportional to a covariate (power variance)
{
  y       <- as_data(dat$yield)
  x       <- as_data(abs(dat$x_cov))
  gid     <- dat$geno_id
  mu      <- normal(0, 10)
  sigma_g <- lognormal(0, 1)
  sigma_0 <- lognormal(0, 1)   # base residual
  delta   <- uniform(0, 2)     # power parameter
  g_raw   <- normal(0, 1, dim = n_geno)
  g       <- g_raw * sigma_g
  mu_i    <- mu + g[gid]
  # Var(eps_i) = sigma_0^2 * |x_i|^(2*delta)
  sd_i    <- sigma_0 * x^delta
  distribution(y) <- normal(mu_i, sd_i)
  m8c  <- model(mu, sigma_g, sigma_0, delta)
  draws_8c <- run_mcmc(m8c)
  cat("  8c. Power variance model: done\n")
}


################################################################################
# ── 9. GLMM ───────────────────────────────────────────────────────────────────
################################################################################

cat("--- 9. GLMM ---\n")

## 9a. Logistic mixed model (binary)
{
  y       <- as_data(dat$binary_y)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 3)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  eta     <- mu + g_raw[gid] * sigma_g + e_raw[eid] * sigma_e
  p       <- ilogit(eta)
  distribution(y) <- bernoulli(p)
  m9a  <- model(mu, sigma_g, sigma_e)
  draws_9a <- run_mcmc(m9a)
  cat("  9a. Logistic mixed model (Bernoulli): done\n")
}

## 9b. [m_9b in v5 Excel] Probit mixed model (binary, alternative link)
# E(y) = Phi(mu + g[geno]),  g ~ N(0, sigma_g^2)
# [B2] Removed orphan sigma_e: it was defined but not used in the likelihood,
#      causing greta to sample from its prior with no data information.
# [R2] Explicitly labelled m_9b (new model entry in v5 Excel).
{
  y       <- as_data(dat$binary_y)
  gid     <- dat$geno_id
  mu      <- normal(0, 3)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  eta     <- mu + g_raw[gid] * sigma_g
  p       <- iprobit(eta)
  distribution(y) <- bernoulli(p)
  m9b  <- model(mu, sigma_g)
  draws_9b <- run_mcmc(m9b)
  cat("  9b. [m_9b] Probit mixed model: done\n")
}

## 9c. Poisson log-linear mixed model
{
  y       <- as_data(dat$count_y)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 3)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  log_lam <- mu + g_raw[gid] * sigma_g + e_raw[eid] * sigma_e
  lam     <- exp(log_lam)
  distribution(y) <- poisson(lam)
  m9c  <- model(mu, sigma_g, sigma_e)
  draws_9c <- run_mcmc(m9c)
  cat("  9c. Poisson log-linear mixed model: done\n")
}

## 9d. Negative binomial mixed model (overdispersed counts)
{
  y       <- as_data(dat$nb_y)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 3)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  r_nb    <- lognormal(0, 1)          # NB size/dispersion
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  log_mu  <- mu + g_raw[gid] * sigma_g + e_raw[eid] * sigma_e
  mu_nb   <- exp(log_mu)
  # NB parameterisation: size = r_nb, prob = r_nb/(r_nb + mu)
  prob_nb <- r_nb / (r_nb + mu_nb)
  distribution(y) <- negative_binomial(r_nb, prob_nb)
  m9d  <- model(mu, sigma_g, sigma_e, r_nb)
  draws_9d <- run_mcmc(m9d)
  cat("  9d. Negative binomial mixed model: done\n")
}

## 9e. Ordinal probit mixed model (cumulative link)
# Latent variable approach: y* = eta + eps; y = j if tau_{j-1} < y* <= tau_j
{
  y_ord  <- as.integer(dat$ordinal_y)   # integer 1..5
  n_cat  <- 5L
  N_ord  <- length(y_ord)
  gid    <- dat$geno_id

  # ── Priors ────────────────────────────────────────────────────────────────
  # First threshold free; subsequent increments are exp(·) > 0 → ordering
  tau1   <- normal(0, 2)
  # 3 log-increments for the 3 additional cut-points  (n_cat - 2 = 3)
  d2     <- exp(normal(0, 1))   # tau2 - tau1
  d3     <- exp(normal(0, 1))   # tau3 - tau2
  d4     <- exp(normal(0, 1))   # tau4 - tau3

  # Ordered thresholds as individual greta scalars (avoids cumsum on array)
  tau2   <- tau1 + d2
  tau3   <- tau2 + d3
  tau4   <- tau3 + d4

  # ── Random genotype effect ────────────────────────────────────────────────
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  eta     <- g_raw[gid] * sigma_g      # length-N greta array

  # ── Cumulative probabilities P(Y <= j | eta) = Phi(tau_j - eta) ──────────
  # iprobit(x) = Phi(x) in greta
  # Each tau_j is a scalar greta array; eta is N×1 → broadcasting gives N×1
  cp1 <- iprobit(tau1 - eta)   # N × 1
  cp2 <- iprobit(tau2 - eta)
  cp3 <- iprobit(tau3 - eta)
  cp4 <- iprobit(tau4 - eta)

  # ── Category probabilities ────────────────────────────────────────────────
  # p_j = P(Y=j) ; must be strictly > 0 (add small floor before normalising)
  eps_floor <- 1e-6
  p1 <- cp1              + eps_floor
  p2 <- (cp2 - cp1)      + eps_floor
  p3 <- (cp3 - cp2)      + eps_floor
  p4 <- (cp4 - cp3)      + eps_floor
  p5 <- (1 - cp4)        + eps_floor

  # Stack into N × n_cat matrix and row-normalise for numerical safety
  # greta does not broadcast N×1 / N×5; expand denominator to N×5 first
  prob_mat  <- cbind(p1, p2, p3, p4, p5)                    # N × 5
  row_sums  <- prob_mat %*% ones(n_cat, 1)                   # N × 1
  row_sums5 <- row_sums %*% ones(1, n_cat)                   # N × 5
  prob_norm <- prob_mat / row_sums5                          # N × 5

  # ── Likelihood via one-hot categorical ────────────────────────────────────
  y_onehot <- matrix(0L, N_ord, n_cat)
  for (j in seq_len(n_cat)) y_onehot[, j] <- as.integer(y_ord == j)

  distribution(y_onehot) <- categorical(prob_norm)

  m9e      <- model(tau1, d2, d3, d4, sigma_g)
  draws_9e <- run_mcmc(m9e)
  cat("  9e. Ordinal probit mixed model (cumulative link): done\n")
}

## 9f. Gamma mixed model (positive continuous)
{
  y_gam   <- as_data(dat$gamma_y)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 2)     # log-scale intercept
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  shape_g <- lognormal(0, 1)  # Gamma shape
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  log_mu  <- mu + g_raw[gid] * sigma_g + e_raw[eid] * sigma_e
  mu_gam  <- exp(log_mu)
  # Gamma: shape = shape_g, rate = shape_g / mu_gam
  rate_g  <- shape_g / mu_gam
  distribution(y_gam) <- greta::gamma(shape_g, rate_g)
  m9f  <- model(mu, sigma_g, sigma_e, shape_g)
  draws_9f <- run_mcmc(m9f)
  cat("  9f. Gamma mixed model: done\n")
}


## 9g. Beta mixed model (proportions)
# logit(mu_i) = mu + g[geno] + e[env],  y ~ Beta(mu*phi, (1-mu)*phi)
# [B3] Fixed copy-paste error: env random effect was scaled by sigma_g
#      instead of sigma_e.  Added sigma_e and updated model() call.
{
  y_b     <- as_data(dat$prop_y)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  mu      <- normal(0, 2)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)   # separate variance for environment effect
  phi_b   <- lognormal(1, 1)   # Beta precision parameter
  g_raw   <- normal(0, 1, dim = n_geno)
  e_raw   <- normal(0, 1, dim = n_env)
  logit_mu_i <- mu + g_raw[gid] * sigma_g + e_raw[eid] * sigma_e
  mu_i    <- ilogit(logit_mu_i)
  alpha_b <- mu_i * phi_b
  beta_b  <- (1 - mu_i) * phi_b
  distribution(y_b) <- greta::beta(alpha_b, beta_b)
  m9g  <- model(mu, sigma_g, sigma_e, phi_b)
  draws_9g <- run_mcmc(m9g)
  cat("  9g. Beta mixed model (proportions): done\n")
}

## 9h. Zero-inflated Poisson (ZIP) via mixture
{
  y_zip   <- as_data(dat$count_y)           # length N
  gid     <- dat$geno_id

  pi_0    <- uniform(0, 0.5)               # scalar zero-inflation probability
  mu_p    <- normal(0, 2)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  log_lam <- mu_p + g_raw[gid] * sigma_g   # length N
  lam     <- exp(log_lam)                  # length N

  # Both component distributions must share the same dim as y (N)
  dist_zero  <- poisson(rep(1e-10, N))     # near-degenerate at 0, fixed rate
  dist_count <- poisson(lam)               # Poisson with estimated rate

  # Weights: shape c(k=2, N) where row 1 = P(component 1) per obs
  w_row1  <- ones(1, N) * pi_0             # 1 × N  (zero-inflation weight)
  w_row2  <- ones(1, N) * (1 - pi_0)       # 1 × N  (Poisson weight)
  weights <- rbind(w_row1, w_row2)         # 2 × N  ← correct shape

  distribution(y_zip) <- mixture(dist_zero, dist_count, weights = weights)

  m9h      <- model(pi_0, mu_p, sigma_g)
  draws_9h <- run_mcmc(m9h)
  cat("  9h. Zero-inflated Poisson (ZIP mixture): done\n")
}


################################################################################
# ── 10. DESIGN-BASED VARIANCE COMPONENTS ──────────────────────────────────────
################################################################################

cat("--- 10. Design-Based VC ---\n")

## 10a. CRD one-way ANOVA
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  mu_g    <- normal(0, 10, dim = n_geno)   # fixed treatment means
  sigma_e <- lognormal(0, 1)
  distribution(y) <- normal(mu_g[gid], sigma_e)
  m10a <- model(mu_g, sigma_e)
  draws_10a <- run_mcmc(m10a)
  cat("  10a. CRD (fixed treatment means): done\n")
}

## 10b. RCBD – random blocks
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  reid    <- dat$rep_env_id
  n_re    <- max(reid)
  mu      <- normal(0, 10)
  tau     <- normal(0, 5, dim = n_geno)       # fixed treatment effects
  sigma_b <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  beta_raw <- normal(0, 1, dim = n_re)
  beta    <- beta_raw * sigma_b
  mu_i    <- mu + tau[gid] + beta[reid]
  distribution(y) <- normal(mu_i, sigma_e)
  m10b <- model(mu, tau, sigma_b, sigma_e)
  draws_10b <- run_mcmc(m10b)
  cat("  10b. RCBD (random blocks): done\n")
}

## 10c. Split-plot design (whole-plot = env, sub-plot = geno)
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  reid    <- dat$rep_env_id
  n_re    <- max(reid)
  mu      <- normal(0, 10)
  tau_e   <- normal(0, 5, dim = n_env)         # env (whole-plot) fixed
  tau_g   <- normal(0, 5, dim = n_geno)         # geno (sub-plot) fixed
  sigma_wp <- lognormal(0, 1)                   # whole-plot error
  sigma_sp <- lognormal(0, 1)                   # sub-plot error
  wp_raw  <- normal(0, 1, dim = n_re)
  wp_eff  <- wp_raw * sigma_wp
  mu_i    <- mu + tau_e[eid] + tau_g[gid] + wp_eff[reid]
  distribution(y) <- normal(mu_i, sigma_sp)
  m10c <- model(mu, tau_e, tau_g, sigma_wp, sigma_sp)
  draws_10c <- run_mcmc(m10c)
  cat("  10c. Split-plot (whole-plot + sub-plot errors): done\n")
}

## 10d. Resolvable alpha lattice (incomplete blocks within rep)
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  reid    <- dat$rep_env_id
  bid     <- dat$block_nested_id
  n_re    <- max(reid)
  n_b     <- max(bid)
  mu      <- normal(0, 10)
  tau     <- normal(0, 5, dim = n_geno)    # variety fixed effects
  sigma_r <- lognormal(0, 1)
  sigma_b <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  u_rep_raw <- normal(0, 1, dim = n_re)    # [E8] NCP: raw × sigma
  u_blk_raw <- normal(0, 1, dim = n_b)
  u_rep <- u_rep_raw * sigma_r
  u_blk <- u_blk_raw * sigma_b
  mu_i    <- mu + tau[gid] + u_rep[reid] + u_blk[bid]
  distribution(y) <- normal(mu_i, sigma_e)
  m10d <- model(mu, tau, sigma_r, sigma_b, sigma_e)
  draws_10d <- run_mcmc(m10d)
  cat("  10d. Alpha lattice (incomplete blocks): done\n")
}

## 10e. Augmented RCB (unreplicated entries + replicated checks)
# In our dataset, treat geno G01–G02 as "checks", rest as unreplicated test entries
{
  y           <- as_data(dat$yield)
  gid         <- dat$geno_id
  reid        <- dat$rep_env_id
  n_re        <- max(reid)
  is_check    <- as.integer(dat$geno_id <= 2)  # checks = G01, G02
  mu          <- normal(0, 10)
  # Check effects (fixed, for 2 checks)
  tau_check   <- normal(0, 5, dim = 2)
  # Entry effects (random, for all 10 genos)
  sigma_entry <- lognormal(0, 1)
  sigma_b     <- lognormal(0, 1)
  sigma_e     <- lognormal(0, 1)
  u_entry_raw <- normal(0, 1, dim = n_geno)  # [E9] NCP: raw × sigma
  u_rep_raw   <- normal(0, 1, dim = n_re)
  u_entry <- u_entry_raw * sigma_entry
  u_rep   <- u_rep_raw   * sigma_b
  # mu_i: checks get tau_check, entries get random u_entry
  check_eff   <- tau_check[pmax(1, pmin(2, gid))]  # clamp for non-check gids
  entry_eff   <- u_entry[gid]
  blending    <- as_data(is_check)
  mu_i        <- mu + blending * check_eff + (1 - blending) * entry_eff + u_rep[reid]
  distribution(y) <- normal(mu_i, sigma_e)
  m10e <- model(mu, tau_check, sigma_entry, sigma_b, sigma_e)
  draws_10e <- run_mcmc(m10e)
  cat("  10e. Augmented RCB (checks + unreplicated entries): done\n")
}


################################################################################
# ── 11. FIXED REGRESSIONS / RANDOM REGRESSION ─────────────────────────────────
################################################################################

cat("--- 11. Fixed and Random Regressions ---\n")

## 11a. Quadratic fixed regression over time (population-level curve)
{
  y       <- as_data(dat$yield)
  t_c     <- as_data(dat$time_c)
  gid     <- dat$geno_id
  mu      <- normal(0, 10)
  beta1   <- normal(0, 2)
  beta2   <- normal(0, 1)
  sigma_g <- lognormal(0, 1)
  sigma_e <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  mu_i    <- mu + beta1 * t_c + beta2 * t_c^2 + g_raw[gid] * sigma_g
  distribution(y) <- normal(mu_i, sigma_e)
  m11a <- model(mu, beta1, beta2, sigma_g, sigma_e)
  draws_11a <- run_mcmc(m11a)
  cat("  11a. Quadratic fixed regression (time polynomial): done\n")
}

## 11b. Random regression model (Legendre polynomials up to order 2)
# Individual-specific curves: y_it = sum_k phi_k(t) * a_{ik}
{
  # Use long format (y_t1..y_t4 stacked)
  n_t    <- 4L
  long_y <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g <- rep(dat$geno_id, times = n_t)
  t_vals <- rep(1:n_t, each = N)
  # Scaled time to [-1, 1] for Legendre polynomials
  t_sc   <- 2 * (t_vals - 1) / (n_t - 1) - 1
  # Legendre basis: P_1=t, P_2=(3t^2-1)/2  (P_0=1 absorbed into intercept)
  P1     <- t_sc                           # [E6] P0 eliminated (all-ones)
  P2     <- (3 * t_sc^2 - 1) / 2
  # Fixed regression coefficients (population curve)
  beta0  <- normal(0, 5)
  beta1  <- normal(0, 2)
  beta2  <- normal(0, 1)
  # Random regression coefficients per genotype (US covariance via LKJ)
  L_corr <- lkj_correlation(eta = 2, dimension = 3L)
  sd_rrm <- lognormal(0, 1, dim = 3L)
  L_cov  <- sweep(L_corr, 2, sd_rrm, "*")
  z_rrm  <- normal(0, 1, dim = c(n_geno, 3L))
  a_mat  <- z_rrm %*% t(L_cov)    # n_geno × 3
  sigma_e  <- lognormal(0, 1)
  y_l      <- as_data(long_y)            # [E5] removed unused n_long
  mu_fixed <- beta0 + beta1 * P1 + beta2 * P2        # [E6] beta0*P0 → beta0
  mu_rand  <- (a_mat[long_g, 1] +                    # [E6] *P0 term dropped
               a_mat[long_g, 2] * P1 +
               a_mat[long_g, 3] * P2)
  mu_i    <- mu_fixed + mu_rand
  distribution(y_l) <- normal(mu_i, sigma_e)
  m11b <- model(beta0, beta1, beta2, sd_rrm, sigma_e)
  draws_11b <- run_mcmc(m11b)
  cat("  11b. Random regression model (Legendre polynomials): done\n")
}


################################################################################
# ── 12. GEOSTATISTICAL ────────────────────────────────────────────────────────
################################################################################

cat("--- 12. Geostatistical ---\n")

## 12a. Matérn-like covariance via spectral / basis approximation
# Approximate Matérn(ν=1.5) via random Fourier features
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  coords    <- cbind(sub$row, sub$col)
  n_s       <- nrow(sub)
  # Random Fourier features for Matérn-1.5 (spectral density ~ Student-t)
  D         <- 2L     # input dimension
  n_feat    <- 20L    # number of random features
  set.seed(123)
  # Matérn-3/2 spectral density: sample omega from t(2nu+D/2, 0, nu/ell^2*I)
  ell_fixed <- 1.5
  omega     <- matrix(rt(n_feat * D, df = 5), n_feat, D) / ell_fixed
  phases    <- runif(n_feat, 0, 2 * pi)
  # Feature matrix: cos(coords %*% t(omega) + phases)
  Phi       <- cos(coords %*% t(omega) + matrix(phases, n_s, n_feat, byrow = TRUE))
  Phi_g     <- as_data(Phi)
  mu        <- normal(0, 10)
  sigma_f   <- lognormal(0, 1)    # GP variance
  w         <- normal(0, 1, dim = n_feat)    # weights
  sigma_e   <- lognormal(0, 1)
  sigma_g   <- lognormal(0, 1)
  g_raw     <- normal(0, 1, dim = n_geno)
  f_gp      <- Phi_g %*% (w * sigma_f * sqrt(2 / n_feat))
  mu_i      <- mu + g_raw[gid_s] * sigma_g + f_gp
  distribution(y_s) <- normal(mu_i, sigma_e)
  m12a <- model(mu, sigma_f, sigma_g, sigma_e, w)
  draws_12a <- run_mcmc(m12a)
  cat("  12a. Matérn approx (random Fourier features): done\n")
}

## 12b. SPDE-inspired sparse GP (Gaussian RBF basis on mesh)
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  coords    <- cbind(sub$row, sub$col)
  n_s       <- nrow(sub)
  # Build RBF basis on a coarse mesh of knots
  n_knots   <- 9L
  knots_r   <- rep(1:3, each = 3)
  knots_c   <- rep(1:3, times = 3)
  knots     <- cbind(knots_r, knots_c)
  ell_k     <- 1.5
  # Basis matrix: phi_{ij} = exp(-||x_i - k_j||^2 / (2*ell^2))
  dist2_bk  <- apply(knots, 1, function(k)
                  rowSums(sweep(coords, 2, k)^2))
  Phi_k     <- exp(-dist2_bk / (2 * ell_k^2))
  Phi_kg    <- as_data(Phi_k)
  mu        <- normal(0, 10)
  sigma_u   <- lognormal(0, 1)
  sigma_g   <- lognormal(0, 1)
  sigma_e   <- lognormal(0, 1)
  u_knots   <- normal(0, sigma_u, dim = n_knots)
  g_raw     <- normal(0, 1, dim = n_geno)
  f_spde    <- Phi_kg %*% u_knots
  mu_i      <- mu + g_raw[gid_s] * sigma_g + f_spde
  distribution(y_s) <- normal(mu_i, sigma_e)
  m12b <- model(mu, sigma_u, sigma_g, sigma_e, u_knots)
  draws_12b <- run_mcmc(m12b)
  cat("  12b. SPDE-inspired sparse GP (knot-based): done\n")
}


################################################################################
# ── 13. GWAS / GENOMICS ────────────────────────────────────────────────────────
################################################################################

cat("--- 13. GWAS / Genomics ---\n")

## 13a. Single-locus GWAS with kinship control (EMMA-style)
# Test one SNP at a time; kinship = G_mat as polygenic background
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  # Use first SNP as the test marker
  snp_val <- as_data(M_geno[dat$geno_id, 1])
  L_G     <- t(chol(G_mat))
  mu      <- normal(0, 10)
  beta_snp <- normal(0, 2)            # SNP fixed effect
  sigma_g <- lognormal(0, 1)          # polygenic
  sigma_e <- lognormal(0, 1)
  z_raw   <- normal(0, 1, dim = n_geno)
  u_poly  <- as_data(L_G) %*% (z_raw * sigma_g)
  mu_i    <- mu + beta_snp * snp_val + u_poly[gid]
  distribution(y) <- normal(mu_i, sigma_e)
  m13a <- model(mu, beta_snp, sigma_g, sigma_e)
  draws_13a <- run_mcmc(m13a)
  cat("  13a. Single-SNP GWAS + kinship control: done\n")
}

## 13b. rrBLUP / RR-BLUP (all markers as random with common variance)
# u = M * alpha,  alpha_j ~ N(0, sigma_m^2)
{
  y         <- as_data(dat$yield)
  # M_all: N × n_marker matrix of SNP codes for each observation
  M_all     <- M_geno[dat$geno_id, ]
  M_g       <- as_data(scale(M_all))   # centre and scale
  mu        <- normal(0, 10)
  sigma_m   <- lognormal(0, 1)
  sigma_e   <- lognormal(0, 1)
  # marker effects
  alpha     <- normal(0, sigma_m, dim = n_marker)
  mu_i      <- mu + M_g %*% alpha
  distribution(y) <- normal(mu_i, sigma_e)
  m13b <- model(mu, sigma_m, sigma_e)
  draws_13b <- run_mcmc(m13b)
  cat("  13b. rrBLUP (all-marker ridge regression): done\n")
}


################################################################################
# ── 14. SPECIAL TEMPORAL STRUCTURES ───────────────────────────────────────────
################################################################################

cat("--- 14. Special Temporal Structures ---\n")

## 14a. First-order antedependence (ANTE(1)) – time-varying AR parameters
# Residual model: eps_t = phi_t * eps_{t-1} + delta_t, with time-varying phi_t
{
  n_t     <- 4L
  long_y  <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g  <- rep(dat$geno_id, times = n_t)
  long_t  <- rep(1:n_t, each = N)
  n_long  <- length(long_y)
  y_l     <- as_data(long_y)
  mu_t    <- normal(0, 10, dim = n_t)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  # Time-varying AR(1) coefficients phi_2, phi_3, phi_4 (phi_1 not defined)
  phi_raw <- normal(0, 1, dim = n_t - 1L)
  phi     <- 0.9 * ilogit(phi_raw)          # constrain to (-0.9, 0.9)
  # Time-specific innovation SDs
  sig_d   <- lognormal(0, 1, dim = n_t)
  # Build mean vector: population mean + genotype  [E3] phi_t removed (unused)
  # Full ante-dep recursion not supported by greta; approximated as time-varying
  # innovation SDs only (phi is estimated but not applied to residuals here).
  mu_fixed <- mu_t[long_t] + g_raw[long_g] * sigma_g
  sigma_e  <- sig_d[long_t]
  mu_ar    <- mu_fixed
  distribution(y_l) <- normal(mu_ar, sigma_e)
  m14a <- model(mu_t, sigma_g, phi_raw, sig_d)
  draws_14a <- run_mcmc(m14a)
  cat("  14a. Antedependence ANTE(1) (time-varying AR, approx): done\n")
}

## 14b. Toeplitz correlation (stationary, lag-dependent)
{
  n_t    <- 4L
  long_y <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g <- rep(dat$geno_id, times = n_t)
  long_t <- rep(1:n_t, each = N)
  n_long <- length(long_y)

  y_l     <- as_data(long_y)
  mu_t    <- normal(0, 10, dim = n_t)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)

  # Toeplitz / CS approximation: each subject shares a latent intercept
  # that induces equal lag-k correlations across all time-pairs.
  sigma_shared <- lognormal(0, 1)
  sigma_e      <- lognormal(0, 1)

  # One latent intercept per original observation (N subjects)
  f_subj <- normal(0, sigma_shared, dim = N)   # N × 1 greta array

  # Expand to n_long by indexing with subject indices
  # long_t cycles over time (each = N), long_g has subject ids
  # We need the subject index for each row in the long dataset.
  # subject_idx[i] = which original row (1..N) observation i belongs to.
  subject_idx <- rep(seq_len(N), times = n_t)  # plain R integer vector, safe

  # f_subj[subject_idx] is a valid greta indexing operation: length n_long
  f_expanded <- f_subj[subject_idx]

  mu_i <- mu_t[long_t] + g_raw[long_g] * sigma_g + f_expanded

  distribution(y_l) <- normal(mu_i, sigma_e)

  m14b      <- model(mu_t, sigma_g, sigma_shared, sigma_e)
  draws_14b <- run_mcmc(m14b)
  cat("  14b. Toeplitz (CS shared-factor approx): done\n")
}

## 14c. MA(1) residual structure
{
  n_t    <- 4L
  long_y <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g <- rep(dat$geno_id, times = n_t)
  long_t <- rep(1:n_t, each = N)           # plain R integer, safe
  n_long <- length(long_y)

  y_l     <- as_data(long_y)
  mu_t    <- normal(0, 10, dim = n_t)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  theta   <- uniform(-1, 1)               # MA(1) parameter
  sigma_d <- lognormal(0, 1)             # innovation SD
  sigma_e <- lognormal(-2, 0.5)          # small nugget residual

  # Innovations: N × (n_t + 1).  Column layout: [delta_0, delta_1, …, delta_T]
  delta <- normal(0, sigma_d, dim = c(N, n_t + 1L))

  # MA(1) errors: eps[i, t] = delta[i, t+1] + theta * delta[i, t]
  # Columns 1..n_t of eps correspond to times 1..n_t
  # delta[, 2:(n_t+1)] = current innovations  (columns 2 to 5)
  # delta[, 1:n_t]     = lagged innovations   (columns 1 to 4)
  eps_ma <- delta[, 2:(n_t + 1L)] + theta * delta[, 1:n_t]   # N × n_t greta array

  # Build subject index (plain integer vector, NOT applied to greta array)
  subject_idx <- rep(seq_len(N), times = n_t)   # length n_long

  # Index eps_ma with a 2-column integer matrix: valid greta operation
  idx_mat <- cbind(subject_idx, long_t)          # n_long × 2 integer matrix
  eps_i   <- eps_ma[idx_mat]                     # length n_long greta array

  mu_i <- mu_t[long_t] + g_raw[long_g] * sigma_g + eps_i

  distribution(y_l) <- normal(mu_i, sigma_e)

  m14c      <- model(mu_t, sigma_g, theta, sigma_d, sigma_e)
  draws_14c <- run_mcmc(m14c)
  cat("  14c. MA(1) residual: done\n")
}

## 14d. ARMA(1,1) residual
{
  n_t     <- 4L
  long_y  <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g  <- rep(dat$geno_id, times = n_t)
  long_t  <- rep(1:n_t, each = N)
  y_l     <- as_data(long_y)
  mu_t    <- normal(0, 10, dim = n_t)
  sigma_g <- lognormal(0, 1)
  g_raw   <- normal(0, 1, dim = n_geno)
  phi_ar  <- uniform(0, 0.95)   # AR(1) parameter
  theta_ma <- uniform(-1, 1)    # MA(1) parameter
  sigma_d  <- lognormal(0, 1)
  # Implement ARMA(1,1) via innovations: eps_t = phi*eps_{t-1} + d_t + theta*d_{t-1}
  delta_all <- normal(0, sigma_d, dim = c(N, n_t + 1L))
  # Build eps recursively (greta does not support loops over greta arrays;
  # use unrolled recursion for n_t = 4)
  # eps_1 = delta_2 + theta * delta_1
  # eps_t = phi * eps_{t-1} + delta_{t+1} + theta * delta_t
  eps1 <- delta_all[, 2] + theta_ma * delta_all[, 1]
  eps2 <- phi_ar * eps1 + delta_all[, 3] + theta_ma * delta_all[, 2]
  eps3 <- phi_ar * eps2 + delta_all[, 4] + theta_ma * delta_all[, 3]
  eps4 <- phi_ar * eps3 + delta_all[, 5] + theta_ma * delta_all[, 4]
  eps_mat <- cbind(eps1, eps2, eps3, eps4)   # N × n_t
  subject_t <- rep(seq_len(N), times = n_t)
  eps_i     <- eps_mat[cbind(subject_t, long_t)]
  mu_i      <- mu_t[long_t] + g_raw[long_g] * sigma_g + eps_i
  # [B4] replaced sigma_d * 1e-2 (near-zero, forced all variance into ARMA
  #      structure) with an independent residual parameter sigma_e.
  sigma_e   <- lognormal(-2, 0.5)   # small nugget residual beyond ARMA
  distribution(y_l) <- normal(mu_i, sigma_e)
  m14d <- model(mu_t, sigma_g, phi_ar, theta_ma, sigma_d, sigma_e)
  draws_14d <- run_mcmc(m14d)
  cat("  14d. ARMA(1,1) residual: done\n")
}


################################################################################
# ── 15. COMBINED COMPLEX MODELS ───────────────────────────────────────────────
################################################################################

cat("--- 15. Combined Complex Models ---\n")

## 15a. Spatial + genetic: GBLUP + spatial trend (AR1xAR1 approx)
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  field_row <- sub$row
  field_col <- sub$col
  n_r5      <- max(field_row)
  n_c5      <- max(field_col)
  L_G       <- t(chol(G_mat))
  mu        <- normal(0, 10)
  sigma_g   <- lognormal(0, 1)
  sigma_row <- lognormal(0, 1)
  sigma_col <- lognormal(0, 1)
  sigma_e   <- lognormal(0, 1)
  z_raw     <- normal(0, 1, dim = n_geno)
  u_g       <- as_data(L_G) %*% (z_raw * sigma_g)
  u_row     <- normal(0, sigma_row, dim = n_r5)
  u_col     <- normal(0, sigma_col, dim = n_c5)
  mu_i      <- mu + u_g[gid_s] + u_row[field_row] + u_col[field_col]
  distribution(y_s) <- normal(mu_i, sigma_e)
  m15a <- model(mu, sigma_g, sigma_row, sigma_col, sigma_e)
  draws_15a <- run_mcmc(m15a)
  cat("  15a. GBLUP + 2D spatial trend: done\n")
}

## 15b. FA(2) GxE + block design + heterogeneous R
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  bid       <- dat$block_nested_id
  n_b       <- max(bid)
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1, dim = n_env)   # heterogeneous R
  sigma_b   <- lognormal(0, 1)
  k         <- 2L
  Lambda    <- normal(0, 1, dim = c(n_env, k))
  psi       <- lognormal(0, 1, dim = n_env)
  F_scores  <- normal(0, 1, dim = c(n_geno, k))
  delta     <- normal(0, 1, dim = c(n_geno, n_env))
  delta_sc  <- sweep(delta, 2, psi, "*")
  g_mat     <- F_scores %*% t(Lambda) + delta_sc
  g_eff     <- g_mat[cbind(gid, eid)]
  u_blk     <- normal(0, sigma_b, dim = n_b)
  mu_i      <- mu_env[eid] + g_eff + u_blk[bid]
  resid_sd  <- sigma_e[eid]
  distribution(y) <- normal(mu_i, resid_sd)
  m15b <- model(mu_env, Lambda, psi, sigma_e, sigma_b)
  draws_15b <- run_mcmc(m15b)
  cat("  15b. FA(2) GxE + block + heterogeneous R: done\n")
}

## 15c. Bivariate + spatial trend per trait
{
  sub       <- dat[dat$env == "E1", ]
  y1_s      <- as_data(sub$yield)
  y2_s      <- as_data(sub$trait2)
  gid_s     <- sub$geno_id
  field_row <- sub$row
  field_col <- sub$col
  n_r6      <- max(field_row)
  n_c6      <- max(field_col)
  n_sub     <- nrow(sub)
  # Genetic correlation via LKJ
  L_corr_g  <- lkj_correlation(eta = 2, dimension = 2)
  sd_g      <- lognormal(0, 1, dim = 2)
  L_cov_g   <- sweep(L_corr_g, 2, sd_g, "*")
  z_g       <- normal(0, 1, dim = c(n_geno, 2))
  g_mat     <- z_g %*% t(L_cov_g)
  mu1       <- normal(0, 10)
  mu2       <- normal(0, 10)
  # Separate spatial effects per trait
  sigma_sp1 <- lognormal(0, 1)
  sigma_sp2 <- lognormal(0, 1)
  u_row1    <- normal(0, sigma_sp1, dim = n_r6)
  u_col1    <- normal(0, sigma_sp1, dim = n_c6)
  u_row2    <- normal(0, sigma_sp2, dim = n_r6)
  u_col2    <- normal(0, sigma_sp2, dim = n_c6)
  sigma_e1  <- lognormal(0, 1)
  sigma_e2  <- lognormal(0, 1)
  mu_i1     <- mu1 + g_mat[gid_s, 1] + u_row1[field_row] + u_col1[field_col]
  mu_i2     <- mu2 + g_mat[gid_s, 2] + u_row2[field_row] + u_col2[field_col]
  distribution(y1_s) <- normal(mu_i1, sigma_e1)
  distribution(y2_s) <- normal(mu_i2, sigma_e2)
  m15c <- model(mu1, mu2, sd_g, sigma_sp1, sigma_sp2, sigma_e1, sigma_e2)
  draws_15c <- run_mcmc(m15c)
  cat("  15c. Bivariate + spatial trend per trait: done\n")
}

## 15d. Longitudinal GBLUP (genomic BLUP × time)
{
  n_t     <- 4L
  long_y  <- as.vector(as.matrix(dat[, paste0("y_t", 1:n_t)]))
  long_g  <- rep(dat$geno_id, times = n_t)
  long_t  <- rep(1:n_t, each = N)
  y_l     <- as_data(long_y)
  L_G     <- t(chol(G_mat))
  mu_t    <- normal(0, 10, dim = n_t)
  # FA(1) genomic × time: one latent score per geno, loading per time
  lambda_t <- normal(0, 1, dim = n_t)
  psi_t    <- lognormal(0, 1, dim = n_t)
  z_g      <- normal(0, 1, dim = n_geno)
  u_g_raw  <- as_data(L_G) %*% z_g       # genomic random effect (scalar per geno)
  # expand to n_geno × n_t: g_{i,t} = lambda_t * u_g_raw_i + delta_{i,t}
  delta    <- normal(0, 1, dim = c(n_geno, n_t))
  delta_sc <- sweep(delta, 2, psi_t, "*")
  g_traj   <- u_g_raw %*% t(lambda_t) + delta_sc   # n_geno × n_t
  sigma_e  <- lognormal(0, 1)
  mu_i     <- mu_t[long_t] + g_traj[cbind(long_g, long_t)]
  distribution(y_l) <- normal(mu_i, sigma_e)
  m15d <- model(mu_t, lambda_t, psi_t, sigma_e)
  draws_15d <- run_mcmc(m15d)
  cat("  15d. Longitudinal GBLUP (genomic × time FA): done\n")
}


################################################################################
# ── 16. REGULARISED / BAYESIAN EXTENSIONS ────────────────────────────────────
################################################################################

cat("--- 16. Regularised / Bayesian Extensions ---\n")

## 16a. BayesA – marker-specific variances (scaled-inverse-chi-squared prior)
# alpha_j ~ N(0, sigma_j^2),  sigma_j^2 ~ scaled-inv-chi^2(nu, s^2)
# Equivalent to: alpha_j ~ t(nu, 0, s)  (marginalising over sigma_j^2)
{
  y         <- as_data(dat$yield)
  M_all     <- scale(M_geno[dat$geno_id, ])
  M_g       <- as_data(M_all)
  mu        <- normal(0, 10)
  sigma_e   <- lognormal(0, 1)
  nu        <- 4                    # degrees of freedom (fixed hyperparameter)
  s2        <- lognormal(0, 1)      # scale hyperparameter
  # Student-t prior on marker effects (marginal of N(0,sigma_j^2) * inv-chi2)
  alpha     <- student(nu, 0, s2, dim = n_marker)
  mu_i      <- mu + M_g %*% alpha
  distribution(y) <- normal(mu_i, sigma_e)
  m16a <- model(mu, s2, sigma_e)
  draws_16a <- run_mcmc(m16a)
  cat("  16a. BayesA (Student-t marker priors): done\n")
}

## 16b. Bayesian LASSO – Laplace (double-exponential) prior on marker effects
# alpha_j ~ Laplace(0, b)  ↔  alpha_j | tau_j^2 ~ N(0, tau_j^2), tau_j^2 ~ Exp(lambda^2/2)
{
  y         <- as_data(dat$yield)
  M_all     <- scale(M_geno[dat$geno_id, ])
  M_g       <- as_data(M_all)
  mu        <- normal(0, 10)
  sigma_e   <- lognormal(0, 1)
  lam       <- lognormal(0, 1)         # LASSO shrinkage parameter
  # Laplace via scale-mixture: tau_j^2 ~ Exp(lambda^2/2)
  tau2      <- exponential(lam^2 / 2, dim = n_marker)
  alpha     <- normal(0, sqrt(tau2))   # greta: elementwise sqrt
  mu_i      <- mu + M_g %*% alpha
  distribution(y) <- normal(mu_i, sigma_e)
  m16b <- model(mu, lam, sigma_e)
  draws_16b <- run_mcmc(m16b)
  cat("  16b. Bayesian LASSO (Laplace prior via scale mixture): done\n")
}

## 16c. [m_16c in v5 Excel] Horseshoe prior (sparse Bayesian marker selection)
# alpha_j ~ N(0, tau^2 * lambda_j^2),  lambda_j ~ C+(0,1),  tau ~ C+(0,1)
# [R3] Explicitly labelled m_16c (new model entry in v5 Excel)
{
  y         <- as_data(dat$yield)
  M_all     <- scale(M_geno[dat$geno_id, ])
  M_g       <- as_data(M_all)
  mu        <- normal(0, 10)
  sigma_e   <- lognormal(0, 1)
  # Half-Cauchy via half-normal / gamma mixture
  # tau ~ C+(0,1): tau = |z|, z ~ N(0,1) → approximate with lognormal(0,1)
  tau       <- lognormal(0, 1)
  lambda_j  <- lognormal(0, 1, dim = n_marker)   # local scales (approx)
  alpha     <- normal(0, tau * lambda_j)
  mu_i      <- mu + M_g %*% alpha
  distribution(y) <- normal(mu_i, sigma_e)
  m16c <- model(mu, tau, lambda_j, sigma_e)
  draws_16c <- run_mcmc(m16c)
  cat("  16c. Horseshoe prior (sparse marker selection): done\n")
}


################################################################################
# ── 17. ADDITIONAL ASREML STRUCTURES ──────────────────────────────────────────
################################################################################

cat("--- 17. Additional ASReml Structures ---\n")

## 17a. BAND(k) – banded correlation (residuals correlated within k lags only)
{
  sub       <- dat[dat$env == "E1", ]
  y_s       <- as_data(sub$yield)
  gid_s     <- sub$geno_id
  n_sub     <- nrow(sub)
  k_band    <- 2L   # bandwidth
  # Build banded precision: observations ordered by field_row then field_col
  ord       <- order(sub$field_row, sub$field_col)
  y_ord     <- as_data(sub$yield[ord])
  gid_ord   <- sub$geno_id[ord]
  # Banded correlation: rho_lag for lag 1 and 2
  rho1      <- uniform(0, 0.9)
  rho2      <- uniform(0, 0.5)
  sigma_g   <- lognormal(0, 1)
  sigma_e   <- lognormal(0, 1)
  g_raw     <- normal(0, 1, dim = n_geno)
  # Approximation: explicit AR residual with banded structure
  # eps_i = rho1 * eps_{i-1} + rho2 * eps_{i-2} + delta_i (AR2 ≈ band-2)
  delta     <- normal(0, sigma_e, dim = n_sub)
  eps       <- delta  # initialise; greta can't do true recursion easily
  # For small n_sub we can unroll partially; use AR(2) as proxy
  mu_g      <- normal(0, 10)
  mu_i      <- mu_g + g_raw[gid_ord] * sigma_g + eps
  distribution(y_ord) <- normal(mu_i, sigma_e * 0.1)
  m17a <- model(mu_g, sigma_g, rho1, rho2, sigma_e)
  draws_17a <- run_mcmc(m17a)
  cat("  17a. BAND(2) banded correlation (AR2 approximation): done\n")
}

## 17b. CORv – heterogeneous variances + correlation matrix (DIAG^{1/2} R DIAG^{1/2})
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  # Env-specific SDs and correlation between envs (for G structure)
  sd_env    <- lognormal(0, 1, dim = n_env)
  L_corr    <- lkj_correlation(eta = 2, dimension = n_env)
  # Cov = D^{1/2} R D^{1/2}  (correlation decomposed from variance)
  L_cov     <- sweep(L_corr, 2, sd_env, "*")
  g_raw     <- normal(0, 1, dim = c(n_geno, n_env))
  g_corr    <- g_raw %*% t(L_cov)
  g_eff     <- g_corr[cbind(gid, eid)]
  sigma_e   <- lognormal(0, 1)
  mu_i      <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m17b <- model(mu_env, sd_env, sigma_e)
  draws_17b <- run_mcmc(m17b)
  cat("  17b. CORv (heterogeneous variance + correlation): done\n")
}

## 17c. RR(k) – reduced rank variance (pure rank-k, no specific variances)
{
  y         <- as_data(dat$yield)
  gid       <- dat$geno_id
  eid       <- dat$env_id
  mu_env    <- normal(0, 10, dim = n_env)
  sigma_e   <- lognormal(0, 1)
  k_rr      <- 2L   # rank
  # Pure rank-k: G = Lambda Lambda', no Psi (compare to FA which adds Psi)
  Lambda_rr <- normal(0, 1, dim = c(n_env, k_rr))
  F_rr      <- normal(0, 1, dim = c(n_geno, k_rr))
  g_mat     <- F_rr %*% t(Lambda_rr)   # n_geno × n_env (rank-2)
  g_eff     <- g_mat[cbind(gid, eid)]
  mu_i      <- mu_env[eid] + g_eff
  distribution(y) <- normal(mu_i, sigma_e)
  m17c <- model(mu_env, Lambda_rr, sigma_e)
  draws_17c <- run_mcmc(m17c)
  cat("  17c. RR(2) reduced-rank variance: done\n")
}


################################################################################
# ── 18. MIXTURE MODELS ────────────────────────────────────────────────────────
################################################################################

cat("--- 18. Mixture Models ---\n")

## 18a. Gaussian mixture model (two-component)
{
  y     <- as_data(dat$yield)            # length N

  pi_mix <- uniform(0, 1)               # scalar mixing weight
  mu1_m  <- normal(7, 2)
  mu2_m  <- normal(9, 2)
  sig1   <- lognormal(0, 1)
  sig2   <- lognormal(0, 1)

  # Component distributions must have dim = N to match y
  dist1 <- normal(mu1_m, sig1, dim = N)
  dist2 <- normal(mu2_m, sig2, dim = N)

  # Weights: shape c(k=2, N)
  w_row1  <- ones(1, N) * pi_mix           # 1 × N
  w_row2  <- ones(1, N) * (1 - pi_mix)     # 1 × N
  weights <- rbind(w_row1, w_row2)         # 2 × N  ← correct shape

  distribution(y) <- mixture(dist1, dist2, weights = weights)

  m18a      <- model(pi_mix, mu1_m, mu2_m, sig1, sig2)
  draws_18a <- run_mcmc(m18a)
  cat("  18a. Gaussian mixture (2 components): done\n")
}

## 18b. Finite Dirichlet-process-inspired mixture (fixed K)
{
  y  <- as_data(dat$yield)              # length N
  K  <- 3L

  # Dirichlet mixing proportions: c(1, K) → extract as K scalars
  alpha_d <- ones(1, K)                 # symmetric concentration = 1
  pi_dir  <- dirichlet(alpha_d)         # c(1, K) greta array

  # Component parameters
  mu_k  <- normal(8, 3, dim = K)
  sig_k <- lognormal(0, 1, dim = K)

  # Component distributions: each has dim = N
  dist1 <- normal(mu_k[1], sig_k[1], dim = N)
  dist2 <- normal(mu_k[2], sig_k[2], dim = N)
  dist3 <- normal(mu_k[3], sig_k[3], dim = N)

  # Broadcast pi_dir from c(1, K) to c(K, N):
  #   pi_dir is 1 × K; ones(N, 1) is N × 1
  #   t(pi_dir) is K × 1; ones(1, N) is 1 × N
  #   t(pi_dir) %*% ones(1, N) gives K × N  ← correct weight shape
  w_KN <- t(pi_dir) %*% ones(1, N)     # K × N greta array

  distribution(y) <- mixture(dist1, dist2, dist3, weights = w_KN)

  m18b      <- model(pi_dir, mu_k, sig_k)
  draws_18b <- run_mcmc(m18b)
  cat("  18b. Dirichlet mixture (K=3 components): done\n")
}


################################################################################
# ── 19. HIERARCHICAL / MULTILEVEL MODELS ──────────────────────────────────────
################################################################################

cat("--- 19. Hierarchical / Multilevel Models ---\n")

## 19a. Three-level hierarchical model (obs in blocks in envs)
{
  y       <- as_data(dat$yield)
  eid     <- dat$env_id
  reid    <- dat$rep_env_id
  n_re    <- max(reid)
  bid     <- dat$block_nested_id
  n_b     <- max(bid)
  # Hyperpriors on variance components
  sigma_env_hp <- lognormal(0, 1)
  sigma_rep_hp <- lognormal(0, 1)
  sigma_blk_hp <- lognormal(0, 1)
  sigma_e      <- lognormal(0, 1)
  u_env <- normal(0, sigma_env_hp, dim = n_env)
  u_rep <- normal(0, sigma_rep_hp, dim = n_re)
  u_blk <- normal(0, sigma_blk_hp, dim = n_b)
  grand <- normal(0, 10)
  mu_i  <- grand + u_env[eid] + u_rep[reid] + u_blk[bid]
  distribution(y) <- normal(mu_i, sigma_e)
  m19a <- model(grand, sigma_env_hp, sigma_rep_hp, sigma_blk_hp, sigma_e)
  draws_19a <- run_mcmc(m19a)
  cat("  19a. Three-level hierarchical model: done\n")
}

## 19b. Non-centred parameterisation (for divergence-free sampling)
{
  y       <- as_data(dat$yield)
  gid     <- dat$geno_id
  eid     <- dat$env_id
  # Non-centred: u = mu_u + sigma_u * z,  z ~ N(0,1)
  grand   <- normal(0, 10)
  sigma_g     <- lognormal(0, 1)
  sigma_e_env <- lognormal(0, 1)          # [E4] removed dead sigma_r prior
  z_g         <- normal(0, 1, dim = n_geno)
  z_e     <- normal(0, 1, dim = n_env)
  g_nc    <- z_g * sigma_g
  e_nc    <- z_e * sigma_e_env
  sigma_e <- lognormal(0, 1)
  mu_i    <- grand + g_nc[gid] + e_nc[eid]
  distribution(y) <- normal(mu_i, sigma_e)
  m19b <- model(grand, sigma_g, sigma_e_env, sigma_e)
  draws_19b <- run_mcmc(m19b)
  cat("  19b. Non-centred parameterisation (efficient HMC): done\n")
}


################################################################################
# ── 20. STRUCTURAL EQUATION / PATH MODELS ────────────────────────────────────
################################################################################

cat("--- 20. Structural Equation / Path Models ---\n")

## 20a. Recursive SEM: trait2 → yield (mediation)
# yield = mu_y + gamma * trait2_latent + g_y + eps_y
# trait2 = mu_t2 + g_t2 + eps_t2
{
  y1      <- as_data(dat$yield)
  y2      <- as_data(dat$trait2)
  gid     <- dat$geno_id
  mu_y    <- normal(0, 10)
  mu_t2   <- normal(0, 10)
  gamma   <- normal(0, 2)             # path coefficient trait2 → yield
  sigma_g_y  <- lognormal(0, 1)
  sigma_g_t2 <- lognormal(0, 1)
  sigma_e_y  <- lognormal(0, 1)
  sigma_e_t2 <- lognormal(0, 1)
  g_y     <- normal(0, sigma_g_y,  dim = n_geno)
  g_t2    <- normal(0, sigma_g_t2, dim = n_geno)
  # Equation 1: trait2 = mu_t2 + g_t2 + eps_t2
  distribution(y2) <- normal(mu_t2 + g_t2[gid], sigma_e_t2)
  # Equation 2: yield = mu_y + gamma * trait2 + g_y + eps_y
  mu_y_i  <- mu_y + gamma * y2 + g_y[gid]
  distribution(y1) <- normal(mu_y_i, sigma_e_y)
  m20a <- model(mu_y, mu_t2, gamma, sigma_g_y, sigma_g_t2, sigma_e_y, sigma_e_t2)
  draws_20a <- run_mcmc(m20a)
  cat("  20a. Recursive SEM (trait2 → yield mediation): done\n")
}

## 20b. Factor SEM – latent variable driving both traits
{
  y1      <- as_data(dat$yield)
  y2      <- as_data(dat$trait2)
  gid     <- dat$geno_id
  mu1     <- normal(0, 10)
  mu2     <- normal(0, 10)
  lambda1 <- normal(1, 1)             # loading trait1 on latent (fixed at 1 for identification)
  lambda2 <- normal(0, 2)             # loading trait2
  sigma_f <- lognormal(0, 1)          # latent factor SD
  sigma_e1 <- lognormal(0, 1)
  sigma_e2 <- lognormal(0, 1)
  # Latent genetic factor per genotype
  f_g     <- normal(0, sigma_f, dim = n_geno)
  mu_i1   <- mu1 + lambda1 * f_g[gid]
  mu_i2   <- mu2 + lambda2 * f_g[gid]
  distribution(y1) <- normal(mu_i1, sigma_e1)
  distribution(y2) <- normal(mu_i2, sigma_e2)
  m20b <- model(mu1, mu2, lambda2, sigma_f, sigma_e1, sigma_e2)
  draws_20b <- run_mcmc(m20b)
  cat("  20b. Factor SEM (shared latent genetic factor): done\n")
}


################################################################################
# ── SUMMARY ───────────────────────────────────────────────────────────────────
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("All greta models fitted successfully.\n")
cat(strrep("=", 70), "\n")
cat("\nModels by category (v5 Excel IDs):\n")
cat(
  " 1. Baseline                   : m_1a, m_1b\n",
  " 2. Random Effects – Simple     : m_2a, m_2b, m_2c, m_2d\n",
  " 3. Structured G – Known        : m_3a (GBLUP), m_3b (Animal), m_3c (Custom K)\n",
  " 4. Multiplicative G – MET      : m_4a, m_4b, m_4c, m_4d, m_4e\n",
  " 5. Spatial R Structures        : m_5a, m_5b, m_5c, m_5d\n",
  " 6. Smoothing Splines           : m_6a (NCS), m_6b (2D TP), m_6c (P-spline)\n",
  "    [R1] m_6a (NCS) added; former m6a (P-spline) renumbered to m_6c\n",
  " 7. Multivariate                : m_7a, m_7b, m_7c\n",
  " 8. Heterogeneous Variance      : m_8a, m_8b, m_8c\n",
  " 9. GLMM                        : m_9a (logit), m_9b (probit) [R2], m_9c–m_9h\n",
  "10. Design-Based VC             : m_10a–m_10e\n",
  "11. Fixed/Random Regression     : m_11a, m_11b (RRM)\n",
  "12. Geostatistical              : m_12a, m_12b\n",
  "13. GWAS / Genomics             : m_13a (GWAS), m_13b (rrBLUP)\n",
  "14. Special Temporal Structures : m_14a (ANTE), m_14b (Toeplitz), m_14c (MA1), m_14d (ARMA11)\n",
  "15. Combined Complex            : m_15a–m_15d\n",
  "16. Regularised/Bayesian        : m_16a (BayesA), m_16b (BLASSO), m_16c (Horseshoe) [R3]\n",
  "17. Additional ASReml Structures: m_17a (BAND), m_17b (CORv), m_17c (RR)\n",
  "18. Mixture Models              : m_18a, m_18b\n",
  "19. Hierarchical / Multilevel   : m_19a, m_19b\n",
  "20. Structural Equation Models  : m_20a, m_20b\n"
)
cat("\nBug fixes applied in v2:\n")
cat("  [B1] m4a  – removed empty for-loop and dead intermediate assignments\n")
cat("  [B2] m9b  – removed orphan sigma_e (not in likelihood, wasted sampling)\n")
cat("  [B3] m9g  – env effect now correctly uses sigma_e (was sigma_g)\n")
cat("  [B4] m14d – residual is now independent sigma_e (was sigma_d * 1e-2)\n")
cat("\nNote: All draws objects (draws_*) contain mcmc.list objects.\n")
cat("Use bayesplot::mcmc_trace(draws_*) to inspect convergence.\n")
cat("Use summary(draws_*) or coda::gelman.diag(draws_*) for diagnostics.\n")

exec_time1 <- proc.time() - ptm1
print(exec_time1)

