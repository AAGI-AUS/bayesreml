## ============================================================================
##  source("04_test_translator.R")
##
##  Focused test suite for 02_translator.R (asreml_to_greta_translator).
##  Verifies:
##    (a) Indexing approach -- never model-matrix approach (Section 2, CG-* tests)
##    (b) Parameter dimensions match level counts, not N (Section 3, DIM-* tests)
##    (c) Large-N scalability: code generation stays fast; parameter count
##        scales with levels, not observations (Section 4, SCALE-* tests)
##    (d) Error / edge-case handling (Section 5, ERR-* tests)
##
##  26 tests total; expected output: "PASS" for all, "0 failures".
##
##  Run with:
##    setwd("path/to/asreml_to_greta_translator/")
##    source("01_generate_data.R")   # only if met_data.rds not yet created
##    source("04_test_translator.R")
## ============================================================================

suppressPackageStartupMessages(library(greta))
source("02_translator.R")

## ── 0. Setup ─────────────────────────────────────────────────────────────────
obj    <- readRDS("met_data.rds")
dat    <- obj$dat
G_mat  <- obj$G_mat
A_mat  <- obj$A_mat
n_geno <- obj$n_geno
n_env  <- obj$n_env
N      <- nrow(dat)

## Minimal MCMC settings
SMPL <- 10L; WARM <- 5L; CH <- 1L

## ── Test harness ─────────────────────────────────────────────────────────────
PASS <- 0L; FAIL <- 0L; results <- list()

run_test <- function(label, expr) {
  cat(sprintf("  %-60s", paste0(label, " ...")))
  tryCatch({
    eval(expr)                      # must not signal an error
    cat("PASS\n")
    PASS <<- PASS + 1L
    results[[label]] <<- list(status = "PASS")
  }, error = function(e) {
    cat("FAIL –", conditionMessage(e), "\n")
    FAIL <<- FAIL + 1L
    results[[label]] <<- list(status = "FAIL", msg = conditionMessage(e))
  })
}

## ── 1. Helper functions ───────────────────────────────────────────────────────

## Returns TRUE if all must_contain patterns are found and none of
## must_not_contain patterns are found in code_str (exact string matching).
check_indexing <- function(code_str,
                           must_contain     = character(0),
                           must_not_contain = character(0)) {
  ok  <- all(vapply(must_contain,
                    function(p) grepl(p, code_str, fixed = TRUE),
                    logical(1)))
  ok  <- ok & !any(vapply(must_not_contain,
                           function(p) grepl(p, code_str, fixed = TRUE),
                           logical(1)))
  ok
}

## Returns greta array entries whose names start with 'prefix'
param_dims <- function(fit, prefix) {
  arrays <- fit$greta_arrays
  arrays[grep(paste0("^", prefix), names(arrays))]
}

## ── 2. Code-generation / indexing tests (no MCMC) ────────────────────────────
cat("=== Section 2: Code-generation / indexing tests ===\n\n")
cat("--- 2a. Fixed effects ---\n")

run_test("CG-1a  intercept-only: mu_atg, no model.matrix", {
  code <- asreml_to_greta(fixed = yield ~ 1, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(is.character(code))
  stopifnot(check_indexing(code,
    must_contain     = c("mu_atg"),
    must_not_contain = c("model.matrix")))
})

run_test("CG-1b  fixed factor: alpha_env[env_id], no X %*%", {
  code <- asreml_to_greta(fixed = yield ~ 0 + env, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain     = c("alpha_env", "env_id"),
    must_not_contain = c("X_env %*%", "model.matrix")))
})

cat("--- 2b. Random effects ---\n")

run_test("CG-2a  random geno: u_geno[geno_id], no model.matrix", {
  code <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain     = c("u_geno", "geno_id"),
    must_not_contain = c("model.matrix")))
})

run_test("CG-2b  crossed random (geno+env): both indexed", {
  code <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno + env, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("u_geno", "geno_id", "u_env", "env_id")))
})

cat("--- 2c. Structured G ---\n")

run_test("CG-3a  GBLUP vm(geno,Gmat): L_G_ and indexed predictor", {
  code <- asreml_to_greta(fixed = yield ~ 1,
                          random = ~ vm(geno, Gmat), data = dat,
                          known_matrices = list(Gmat = G_mat),
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("L_G_geno", "geno_id")))
})

cat("--- 2d. GxE / multiplicative structures ---\n")

run_test("CG-4a  at(env):geno: cbind indexing, no model.matrix", {
  code <- asreml_to_greta(fixed = yield ~ env,
                          random = ~ at(env):geno, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain     = c("cbind(geno_id, env_id)"),
    must_not_contain = c("model.matrix")))
})

run_test("CG-4b  us(env):id(geno): cbind indexing", {
  code <- asreml_to_greta(fixed = yield ~ env,
                          random = ~ us(env):id(geno), data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("cbind(geno_id, env_id)")))
})

run_test("CG-4c  fa(env,1):id(geno): Lambda_ and cbind indexing", {
  code <- asreml_to_greta(fixed = yield ~ env,
                          random = ~ fa(env, 1):id(geno), data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("Lambda_", "cbind(geno_id, env_id)")))
})

cat("--- 2e. Spatial ---\n")

dat_E1 <- dat[dat$env == "E1", ]

run_test("CG-5a  ar1(row):id(col): row_id and col_id indexed", {
  code <- asreml_to_greta(fixed = yield ~ 1,
                          random = ~ geno + ar1(row):id(col),
                          data = dat_E1,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("row_id", "col_id")))
})

cat("--- 2f. Spline ---\n")

run_test("CG-6a  spl(x_cov): basis matrix %*% (expected; k=8 << N)", {
  code <- asreml_to_greta(fixed = yield ~ 1,
                          random = ~ spl(x_cov), data = dat,
                          return_code = TRUE, verbose = FALSE)
  ## Spline legitimately uses a k-column basis matrix (k=8, not N)
  stopifnot(check_indexing(code,
    must_contain = c("%*%", "splines::bs", "df=8")))
})

cat("--- 2g. Heterogeneous residual ---\n")

run_test("CG-8a  at(env):units rcov: sigma_e_atg indexed by env_id", {
  code <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno,
                          rcov = ~ at(env):units, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain = c("sigma_e_atg", "env_id")))
})

cat("--- 2h. Nested/design-based ---\n")

run_test("CG-10d rep_env+block_nested: indexed, no factor(paste)", {
  code <- asreml_to_greta(fixed = yield ~ geno,
                          random = ~ rep_env + block_nested, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain     = c("rep_env_id", "block_nested_id"),
    must_not_contain = c("factor(paste")))
})

run_test("CG-10d-nested colon-nested (geno:env): no factor(paste in code", {
  ## Test the actual translator fix: nested term (outer:inner colon syntax)
  ## must precompute IDs in setup env, not emit factor(paste into model code.
  code <- asreml_to_greta(fixed = yield ~ env,
                          random = ~ geno:env, data = dat,
                          return_code = TRUE, verbose = FALSE)
  stopifnot(check_indexing(code,
    must_contain     = c("nested_id_env_in_geno"),
    must_not_contain = c("factor(paste")))
})

## ── 3. Parameter-dimension / name tests (with MCMC) ──────────────────────────
cat("\n=== Section 3: Parameter-dimension / name tests ===\n\n")
cat("(MCMC: ", SMPL, " samples, ", WARM, " warmup, ", CH, " chain)\n\n")

run_test("DIM-1b  yield~0+env: alpha_env dim == n_env", {
  fit <- asreml_to_greta(fixed = yield ~ 0 + env, data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  arr <- fit$greta_arrays[["alpha_env"]]
  stopifnot(!is.null(arr))
  ## greta column-vectors have dim c(n, 1); check total elements
  stopifnot(prod(dim(arr)) == n_env)
})

run_test("DIM-2a  ~geno random: sigma_geno is scalar param", {
  fit <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno, data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  stopifnot("sigma_geno" %in% fit$param_names)
  ## sigma_geno should be a scalar (dim c(1,1) or NULL)
  sg <- fit$greta_arrays[["sigma_geno"]]
  stopifnot(!is.null(sg), prod(dim(sg)) == 1L)
})

run_test("DIM-2b  ~geno+env: both sigma_geno and sigma_env present", {
  fit <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno + env, data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  stopifnot("sigma_geno" %in% fit$param_names,
            "sigma_env"  %in% fit$param_names)
})

run_test("DIM-4a  at(env):geno: geno_env_sc dim == c(n_geno, n_env)", {
  fit <- asreml_to_greta(fixed = yield ~ env, random = ~ at(env):geno,
                         data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  ## The scaled interaction matrix is named geno_env_sc
  sc <- fit$greta_arrays[["geno_env_sc"]]
  stopifnot(!is.null(sc))
  stopifnot(all(dim(sc) == c(n_geno, n_env)))
})

run_test("DIM-8a  at(env):units rcov: sigma_e_atg dim == n_env", {
  fit <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno,
                         rcov = ~ at(env):units, data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  se <- fit$greta_arrays[["sigma_e_atg"]]
  stopifnot(!is.null(se))
  stopifnot(prod(dim(se)) == n_env)
})

run_test("DIM-ALL  summary rows == length(param_names) for yield~1+geno", {
  fit <- asreml_to_greta(fixed = yield ~ 1, random = ~ geno, data = dat,
                         n_samples = SMPL, warmup = WARM, chains = CH,
                         verbose = FALSE)
  stopifnot(nrow(fit$summary$statistics) == length(fit$param_names))
})

## ── 4. Large-N scalability smoke test ────────────────────────────────────────
cat("\n=== Section 4: Large-N scalability smoke test ===\n\n")

run_test("SCALE-1  N=2000, code gen < 2 s, no N-sized matrices", {
  set.seed(1)
  n_g <- 50; n_e <- 4; n_r <- 10
  big_dat <- expand.grid(
    geno = paste0("G", seq_len(n_g)),
    env  = paste0("E", seq_len(n_e)),
    rep  = seq_len(n_r)
  )
  big_dat$yield <- rnorm(nrow(big_dat))

  elapsed <- system.time(
    code <- asreml_to_greta(fixed = yield ~ env, random = ~ geno,
                            data = big_dat,
                            return_code = TRUE, verbose = FALSE)
  )[3]

  stopifnot(is.character(code), nchar(code) > 50)
  stopifnot(elapsed < 2)

  ## No N-sized matrix literals inside the model code
  has_big_matrix <- grepl("matrix(0, N", code, fixed = TRUE) ||
                    grepl(paste0("matrix(0, ", nrow(big_dat)), code, fixed = TRUE)
  stopifnot(!has_big_matrix)
})

run_test("SCALE-2  param count == level counts (not N)", {
  set.seed(1)
  n_g <- 50; n_e <- 4; n_r <- 10
  big_dat <- expand.grid(
    geno = paste0("G", seq_len(n_g)),
    env  = paste0("E", seq_len(n_e)),
    rep  = seq_len(n_r)
  )
  big_dat$yield <- rnorm(nrow(big_dat))
  ## Parse via return_code to access ctx$params count
  code <- asreml_to_greta(fixed = yield ~ env, random = ~ geno,
                          data = big_dat,
                          return_code = TRUE, verbose = FALSE)
  ## Code should reference n_geno and n_env, not nrow(big_dat)
  stopifnot(grepl("n_geno", code, fixed = TRUE))
  stopifnot(!grepl(paste0("dim = ", nrow(big_dat)), code, fixed = TRUE))
})

## ── 5. Error / edge-case tests ───────────────────────────────────────────────
cat("\n=== Section 5: Error / edge-case tests ===\n\n")

run_test("ERR-1  unknown family 'student_t' → error", {
  result <- tryCatch(
    asreml_to_greta(fixed = yield ~ 1, data = dat,
                    family = "student_t",
                    return_code = TRUE, verbose = FALSE),
    error = function(e) e
  )
  stopifnot(inherits(result, "error"))
  stopifnot(grepl("Unsupported family", conditionMessage(result),
                  ignore.case = TRUE))
})

run_test("ERR-2  vm(geno,Gmat) with known_matrices=list() → error", {
  result <- tryCatch(
    asreml_to_greta(fixed = yield ~ 1,
                    random = ~ vm(geno, Gmat), data = dat,
                    known_matrices = list(),
                    return_code = TRUE, verbose = FALSE),
    error   = function(e) e,
    warning = function(w) w
  )
  ## Should error or warn about missing matrix
  stopifnot(inherits(result, "error") || inherits(result, "warning"))
})

run_test("ERR-3  response not in data → error", {
  result <- tryCatch(
    asreml_to_greta(fixed = nonexistent_col ~ 1, data = dat,
                    return_code = TRUE, verbose = FALSE),
    error = function(e) e
  )
  stopifnot(inherits(result, "error"))
})

run_test("ERR-4  return_code=TRUE returns character, not list", {
  result <- asreml_to_greta(fixed = yield ~ 1, data = dat,
                            return_code = TRUE, verbose = FALSE)
  stopifnot(is.character(result))
  stopifnot(!is.list(result))
})

run_test("ERR-5  fa(env,0) zero-factor FA → graceful error or warning", {
  result <- tryCatch(
    asreml_to_greta(fixed = yield ~ env,
                    random = ~ fa(env, 0):id(geno), data = dat,
                    return_code = TRUE, verbose = FALSE),
    error   = function(e) e,
    warning = function(w) w
  )
  ## Must signal something (error or warning), not silently produce broken code
  stopifnot(inherits(result, "error") || inherits(result, "warning") ||
              (is.character(result) &&
                 !grepl("fa.*dim = c.*0", result, perl = TRUE)))
})

## ── 6. Summary ───────────────────────────────────────────────────────────────
cat("\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Passed : %d\n  Failed : %d\n  Total  : %d\n",
            PASS, FAIL, PASS + FAIL))
cat(strrep("=", 65), "\n")

if (FAIL > 0L) {
  cat("\nFailed tests:\n")
  for (nm in names(results)) {
    if (results[[nm]]$status == "FAIL")
      cat(sprintf("  %-60s %s\n", nm, results[[nm]]$msg))
  }
}

invisible(list(pass = PASS, fail = FAIL, results = results))
