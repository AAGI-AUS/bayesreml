## ============================================================================
##  bayesreml -- Comprehensive Model Gallery
##
##  This script demonstrates every model type, family, random-effect structure,
##  and residual structure supported by bayesreml. Each model is generated
##  using return_code = TRUE so the generated greta code can be inspected
##  without requiring TensorFlow. Section 6 shows how to actually fit a model
##  and use all downstream methods.
##
##  Usage:
##    source("all_models.R")
##    # Or run interactively section by section
## ============================================================================

library(bayesreml)

## -- Setup: load bundled dataset ----------------------------------------------

data(met_example)
dat    <- met_example$dat
G_mat  <- met_example$G_mat
A_mat  <- met_example$A_mat

cat("Dataset:", nrow(dat), "obs |",
    met_example$n_geno, "genotypes |",
    met_example$n_env, "environments\n")
cat("Columns:", paste(names(dat), collapse = ", "), "\n\n")

## Helper to show a model's generated code with a header
show_model <- function(title, ...) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(" ", title, "\n")
  cat(strrep("=", 70), "\n")
  code <- bayesreml(..., return_code = TRUE, verbose = FALSE,
                    n_samples = 500, warmup = 200, chains = 2)
  cat(code, "\n")
  invisible(code)
}


## ============================================================================
##  SECTION 1: FIXED-EFFECT STRUCTURES
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 1: Fixed-Effect Structures\n")
cat(strrep("#", 70), "\n")

## 1a. Intercept-only
show_model("1a. Intercept-only (yield ~ 1)",
  fixed = yield ~ 1, data = dat)

## 1b. Single factor (no intercept -- cell-means parameterisation)
show_model("1b. Factor no-intercept (yield ~ 0 + env)",
  fixed = yield ~ 0 + env, data = dat)

## 1c. Factor with intercept (effects parameterisation)
show_model("1c. Factor with intercept (yield ~ env)",
  fixed = yield ~ env, data = dat)

## 1d. Continuous covariate
show_model("1d. Continuous covariate (yield ~ x_cov)",
  fixed = yield ~ x_cov, data = dat)

## 1e. Factor + continuous
show_model("1e. Factor + continuous (yield ~ env + x_cov)",
  fixed = yield ~ env + x_cov, data = dat)

## 1f. Factor interaction in fixed effects
show_model("1f. Factor interaction (yield ~ env:geno)",
  fixed = yield ~ env:geno, data = dat)

## 1g. I() expression (polynomial)
show_model("1g. Polynomial via I() (yield ~ x_cov + I(x_cov^2))",
  fixed = yield ~ x_cov + I(x_cov^2), data = dat)

## 1h. Continuous interaction
dat$x2 <- rnorm(nrow(dat))
show_model("1h. Continuous interaction (yield ~ x_cov:x2)",
  fixed = yield ~ x_cov:x2, data = dat)
dat$x2 <- NULL  # clean up


## ============================================================================
##  SECTION 2: RANDOM-EFFECT STRUCTURES (Gaussian family)
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 2: Random-Effect Structures\n")
cat(strrep("#", 70), "\n")

## 2a. Simple IID random intercept
show_model("2a. Simple random intercept (~ geno)",
  fixed = yield ~ env, random = ~ geno, data = dat)

## 2b. Crossed random effects (two IID terms)
show_model("2b. Crossed random effects (~ geno + rep)",
  fixed = yield ~ env, random = ~ geno + rep, data = dat)

## 2c. Nested random effect (block nested within env)
show_model("2c. Nested random (~ env:block)",
  fixed = yield ~ 1, random = ~ env:block, data = dat)

## 2d. Three-way combination (geno:env:rep)
show_model("2d. Three-way nesting (~ geno:env:rep)",
  fixed = yield ~ 1, random = ~ geno:env:rep, data = dat)

## 2e. GBLUP -- vm(geno, Gmat) with genomic relationship matrix
show_model("2e. GBLUP with vm() (~ vm(geno, Gmat))",
  fixed = yield ~ env, random = ~ vm(geno, Gmat),
  data = dat, known_matrices = list(Gmat = G_mat))

## 2f. Pedigree-based animal model -- ped(geno, Amat)
show_model("2f. Animal model with ped() (~ ped(geno, Amat))",
  fixed = yield ~ env, random = ~ ped(geno, Amat),
  data = dat, known_matrices = list(Amat = A_mat))

## 2g. Diagonal GxE -- at(env):geno
show_model("2g. Diagonal GxE (~ at(env):geno)",
  fixed = yield ~ env, random = ~ at(env):geno, data = dat)

## 2h. Unstructured GxE -- us(env):id(geno)
show_model("2h. Unstructured GxE (~ us(env):id(geno))",
  fixed = yield ~ env, random = ~ us(env):id(geno), data = dat)

## 2i. Factor-analytic GxE -- fa(env,2):id(geno)
show_model("2i. Factor-analytic GxE, k=2 (~ fa(env,2):id(geno))",
  fixed = yield ~ env, random = ~ fa(env, 2):id(geno), data = dat)

## 2j. Factor-analytic GxE -- fa(env,1):id(geno) (single factor)
show_model("2j. Factor-analytic GxE, k=1 (~ fa(env,1):id(geno))",
  fixed = yield ~ env, random = ~ fa(env, 1):id(geno), data = dat)

## 2k. Spatial AR1 approximation -- ar1(row):id(col)
show_model("2k. Spatial AR1 (~ ar1(row):id(col))",
  fixed = yield ~ env, random = ~ ar1(row):id(col), data = dat)

## 2l. P-spline smooth -- spl(x_cov)
show_model("2l. P-spline smooth (~ spl(x_cov))",
  fixed = yield ~ env, random = ~ spl(x_cov), data = dat)

## 2m. Genomic GxE -- vm(geno, Gmat):id(env)
show_model("2m. Genomic GxE (~ vm(geno, Gmat):id(env))",
  fixed = yield ~ 1, random = ~ vm(geno, Gmat):id(env),
  data = dat, known_matrices = list(Gmat = G_mat))


## ============================================================================
##  SECTION 3: RESIDUAL (RCOV) STRUCTURES
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 3: Residual (rcov) Structures\n")
cat(strrep("#", 70), "\n")

## 3a. Default IID residuals (~ units, implicit)
show_model("3a. Homogeneous residual (default ~ units)",
  fixed = yield ~ env, random = ~ geno, data = dat)

## 3b. Heterogeneous residual variance by environment
show_model("3b. Heterogeneous residual (rcov = ~ at(env):units)",
  fixed = yield ~ env, random = ~ geno,
  rcov = ~ at(env):units, data = dat)


## ============================================================================
##  SECTION 4: FAMILY / LINK FUNCTIONS
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 4: Family / Link Functions\n")
cat(strrep("#", 70), "\n")

## 4a. Gaussian (default)
show_model("4a. Gaussian / identity (default)",
  fixed = yield ~ env, random = ~ geno, data = dat,
  family = "gaussian")

## 4b. Binomial (logit link)
show_model("4b. Binomial / logit",
  fixed = bin_y ~ env, random = ~ geno, data = dat,
  family = "binomial")

## 4c. Binomial (probit link)
show_model("4c. Binomial / probit",
  fixed = bin_y ~ env, random = ~ geno, data = dat,
  family = "binomial", link = "probit")

## 4d. Poisson (log link)
show_model("4d. Poisson / log",
  fixed = count_y ~ env, random = ~ geno, data = dat,
  family = "poisson")

## 4e. Negative binomial (log link)
show_model("4e. Negative binomial / log",
  fixed = count_y ~ env, random = ~ geno, data = dat,
  family = "negative_binomial")

## 4f. Gamma (log link)
dat$gamma_y <- rgamma(nrow(dat), shape = 2, rate = 0.5)
show_model("4f. Gamma / log",
  fixed = gamma_y ~ env, random = ~ geno, data = dat,
  family = "gamma")

## 4g. Beta (logit link)
dat$beta_y <- rbeta(nrow(dat), shape1 = 2, shape2 = 5)
show_model("4g. Beta / logit",
  fixed = beta_y ~ env, random = ~ geno, data = dat,
  family = "beta")


## ============================================================================
##  SECTION 5: COMBINED / COMPLEX MODELS
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 5: Combined / Complex Models\n")
cat(strrep("#", 70), "\n")

## 5a. Fixed + multiple random + heterogeneous residual
show_model("5a. Multi-component (env + x_cov | geno + env:block | at(env):units)",
  fixed = yield ~ env + x_cov,
  random = ~ geno + env:block,
  rcov = ~ at(env):units,
  data = dat)

## 5b. GBLUP + environment-specific residual
show_model("5b. GBLUP + heterogeneous residual",
  fixed = yield ~ env,
  random = ~ vm(geno, Gmat),
  rcov = ~ at(env):units,
  data = dat,
  known_matrices = list(Gmat = G_mat))

## 5c. Diagonal GxE + spline + nested block
show_model("5c. Diagonal GxE + spline + nested block",
  fixed = yield ~ env + x_cov,
  random = ~ at(env):geno + spl(x_cov) + env:block,
  data = dat)

## 5d. GBLUP + diagonal GxE in same model
show_model("5d. Additive genetics (vm) + diagonal GxE (at)",
  fixed = yield ~ env,
  random = ~ vm(geno, Gmat) + at(env):geno,
  data = dat,
  known_matrices = list(Gmat = G_mat))

## 5e. Pedigree model + crossed random + custom priors
show_model("5e. Pedigree + rep + custom priors",
  fixed = yield ~ env,
  random = ~ ped(geno, Amat) + rep,
  data = dat,
  known_matrices = list(Amat = A_mat),
  prior_fixed_sd = 20, prior_vc_sd = 2)

## 5f. Binomial GLMM with random intercept
show_model("5f. Binomial GLMM (bin_y ~ env | geno)",
  fixed = bin_y ~ env, random = ~ geno, data = dat,
  family = "binomial")

## 5g. Poisson GLMM with nested random
show_model("5g. Poisson GLMM + nested random",
  fixed = count_y ~ env, random = ~ geno + env:block, data = dat,
  family = "poisson")

## 5h. Weighted Gaussian model
wts <- runif(nrow(dat), 0.5, 2.0)
show_model("5h. Weighted Gaussian model",
  fixed = yield ~ env, random = ~ geno, data = dat,
  weights = wts)

## 5i. Spatial model + random genotype
show_model("5i. Spatial AR1 + genotype random",
  fixed = yield ~ env,
  random = ~ geno + ar1(row):id(col),
  data = dat)


## ============================================================================
##  SECTION 6: FITTING A MODEL AND DOWNSTREAM ANALYSIS
##  (requires greta + TensorFlow -- skipped gracefully if unavailable)
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 6: Fitting & Downstream Methods\n")
cat(strrep("#", 70), "\n")

greta_ok <- tryCatch({
  requireNamespace("greta", quietly = TRUE) &&
  requireNamespace("tensorflow", quietly = TRUE) &&
  !is.null(tensorflow::tf_version()) &&
  tensorflow::tf_version() >= "2.0" &&
  # Verify greta can actually create arrays (TF backend works)
  { x <- greta::normal(0, 1); TRUE }
}, error = function(e) FALSE)

if (greta_ok) {
  cat("\ngreta + TensorFlow detected -- fitting a live model...\n\n")

  fit <- bayesreml(
    fixed     = yield ~ env,
    random    = ~ geno,
    data      = dat,
    n_samples = 500, warmup = 200, chains = 2,
    verbose   = TRUE
  )

  ## -- S3 methods on bayesreml object --
  cat("\n--- print() ---\n")
  print(fit)

  cat("\n--- summary() ---\n")
  summary(fit)

  cat("\n--- coef() ---\n")
  print(coef(fit))

  cat("\n--- vcov() ---\n")
  print(round(vcov(fit), 4))

  cat("\n--- confint() ---\n")
  print(confint(fit))

  cat("\n--- fitted() (first 10) ---\n")
  print(head(fitted(fit), 10))

  cat("\n--- residuals() (first 10) ---\n")
  print(head(residuals(fit), 10))

  cat("\n--- nobs() ---\n")
  print(nobs(fit))

  cat("\n--- family() ---\n")
  print(family(fit))

  cat("\n--- formula() ---\n")
  print(formula(fit))

  cat("\n--- logLik() ---\n")
  print(logLik(fit))

  cat("\n--- model.matrix() (first 5 rows) ---\n")
  print(head(model.matrix(fit), 5))

  cat("\n--- predict() (first 10) ---\n")
  print(head(predict(fit), 10))

  cat("\n--- predict(se.fit = TRUE) ---\n")
  p <- predict(fit, se.fit = TRUE)
  cat("  fit:", head(p$fit, 5), "\n")
  cat("  se: ", head(p$se.fit, 5), "\n")

  ## -- S3 methods on bayesreml_glm sub-object --
  cat("\n--- summary(fit$glm) ---\n")
  summary(fit$glm)

  cat("\n--- confint(fit$glm) ---\n")
  print(confint(fit$glm))

  ## -- Broom-style tidiers --
  cat("\n--- tidy() ---\n")
  print(tidy.bayesreml(fit))

  cat("\n--- tidy(effects = 'random') ---\n")
  print(tidy.bayesreml(fit, effects = "random"))

  cat("\n--- glance() ---\n")
  print(glance.bayesreml(fit))

  cat("\n--- augment() (first 5 rows) ---\n")
  print(head(augment.bayesreml(fit), 5))

  ## -- Variance components & BLUPs --
  cat("\n--- Variance components ---\n")
  print(fit$extras$variance_comps)

  cat("\n--- BLUPs (genotype) ---\n")
  if (length(fit$extras$blups) > 0) {
    for (nm in names(fit$extras$blups)) {
      cat("  ", nm, ": ", paste(round(as.numeric(fit$extras$blups[[nm]]), 3),
                                collapse = ", "), "\n")
    }
  } else {
    cat("  (none)\n")
  }

  ## -- Convergence diagnostics --
  cat("\n--- Convergence: effective sample sizes ---\n")
  print(fit$extras$convergence$n_eff)

  ## -- Posterior draws (native greta output) --
  cat("\n--- Posterior draws (class) ---\n")
  cat("Class:", class(fit$greta$draws), "\n")
  cat("Dimensions:", dim(fit$greta$draws[[1]]), "\n")

  ## -- Generated code --
  cat("\n--- Generated greta code ---\n")
  cat(fit$extras$code, "\n")

  ## -- Diagnostic plots --
  cat("\n--- Diagnostic plots ---\n")
  plot(fit, type = "residuals")
  plot(fit, type = "effects")
  plot(fit, type = "variance")
  plot(fit, type = "blups")
  plot(fit, type = "pp_check")

  ## -- Model comparison --
  cat("\n--- Model comparison (anova) ---\n")
  fit2 <- bayesreml(
    fixed     = yield ~ env,
    random    = ~ geno + env:block,
    data      = dat,
    n_samples = 500, warmup = 200, chains = 2,
    verbose   = FALSE
  )
  anova(fit, fit2)

  ## -- emmeans integration --
  if (requireNamespace("emmeans", quietly = TRUE)) {
    cat("\n--- emmeans ---\n")
    tryCatch({
      em <- emmeans::emmeans(fit, ~ env)
      print(em)
    }, error = function(e) cat("  emmeans error:", conditionMessage(e), "\n"))
  }

  cat("\n--- Model info ---\n")
  print(fit$extras$model_info)

} else {
  cat("\ngreta/TensorFlow not available -- skipping live model fitting.\n")
  cat("All", 30, "model code-generation examples above completed successfully.\n")
  cat("To fit models, install greta:\n")
  cat("  install.packages('greta')\n")
  cat("  greta::install_greta_deps()\n")
}


## ============================================================================
##  SECTION 7: METHODS DEMONSTRATION WITH MOCK OBJECT
##  (always runs -- no greta required)
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 7: S3 Methods (mock object, no greta needed)\n")
cat(strrep("#", 70), "\n")

## Build a mock bayesreml object to demonstrate all methods
n <- nrow(dat)
n_samples <- 100
draws_mat <- matrix(rnorm(n_samples * 3), ncol = 3)
draws_mat[, 1] <- draws_mat[, 1] + 50
draws_mat[, 2] <- abs(draws_mat[, 2])
draws_mat[, 3] <- abs(draws_mat[, 3]) + 1
colnames(draws_mat) <- c("mu_atg", "sigma_geno", "sigma_e_atg")
draws <- coda::mcmc.list(coda::mcmc(draws_mat))

beta <- c("(Intercept)" = 50.1)
V <- matrix(0.5, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))
fitted_vals <- rep(50.1, n)
resid_vals  <- dat$yield - fitted_vals

glm_obj <- list(
  coefficients      = beta,
  residuals         = resid_vals,
  fitted.values     = fitted_vals,
  family            = gaussian(),
  formula           = yield ~ 1,
  terms             = terms(yield ~ 1, data = dat),
  model             = model.frame(yield ~ 1, data = dat),
  data              = dat,
  df.residual       = n - 1L,
  rank              = 1L,
  y                 = dat$yield,
  linear.predictors = fitted_vals,
  call              = quote(bayesreml(fixed = yield ~ 1, data = dat)),
  qr                = qr(matrix(1, n, 1))
)
attr(glm_obj, "posterior_vcov") <- V
class(glm_obj) <- c("bayesreml_glm", "glm", "lm")

greta_out <- structure(
  list(model = NULL, draws = draws,
       greta_arrays = list(), env = new.env(parent = emptyenv())),
  class = "bayesreml_greta"
)

post_summary <- summary(draws)
vc_table <- data.frame(
  component = c("sigma_geno", "sigma_e_atg"),
  estimate  = c(1.0, 2.0), sd = c(0.3, 0.2),
  q2.5 = c(0.5, 1.5), q50 = c(1.0, 2.0), q97.5 = c(1.5, 2.5),
  stringsAsFactors = FALSE
)

extras <- structure(list(
  summary        = post_summary,
  convergence    = list(n_eff = coda::effectiveSize(draws), gelman = NULL),
  variance_comps = vc_table,
  blups          = list(u_geno = setNames(rnorm(5), paste0("G", 1:5))),
  predictions    = data.frame(obs = 1:n, observed = dat$yield,
                              fitted = fitted_vals, residual = resid_vals),
  code           = "mu_atg <- normal(0, 10)\nsigma_e_atg <- lognormal(0, 1)",
  param_names    = c("mu_atg", "sigma_geno", "sigma_e_atg"),
  parse_info     = list(
    fixed  = list(response = "yield", intercept = TRUE, terms = list()),
    random = list(), rcov = list(list(type = "units")),
    family = list(family = "gaussian", link = "identity")
  ),
  call_info = list(
    fixed = yield ~ 1, random = NULL, rcov = NULL,
    data_name = "dat", family = "gaussian", link = NULL,
    known_matrices = list(), weights = NULL,
    n_samples = 100, warmup = 50, chains = 1,
    prior_fixed_sd = 10, prior_vc_sd = 1
  ),
  run_time   = 1.5,
  model_info = list(
    n_obs = n, n_fixed = 1L, n_random = 0L, n_params = 3L,
    family = "gaussian", link = "identity"
  )
), class = "bayesreml_extras")

mock_fit <- structure(
  list(glm = glm_obj, greta = greta_out, extras = extras),
  class = "bayesreml"
)

## Now demonstrate every S3 method
cat("\n--- print ---\n")
print(mock_fit)

cat("\n--- summary ---\n")
summary(mock_fit)

cat("\n--- coef ---\n")
print(coef(mock_fit))

cat("\n--- vcov ---\n")
print(vcov(mock_fit))

cat("\n--- confint ---\n")
print(confint(mock_fit))

cat("\n--- fitted (first 6) ---\n")
print(head(fitted(mock_fit)))

cat("\n--- residuals (first 6) ---\n")
print(head(residuals(mock_fit)))

cat("\n--- nobs ---\n")
print(nobs(mock_fit))

cat("\n--- family ---\n")
print(family(mock_fit))

cat("\n--- formula ---\n")
print(formula(mock_fit))

cat("\n--- logLik ---\n")
print(logLik(mock_fit))

cat("\n--- model.matrix (first 3 rows) ---\n")
print(head(model.matrix(mock_fit), 3))

cat("\n--- predict (first 6) ---\n")
print(head(predict(mock_fit)))

cat("\n--- predict with se.fit ---\n")
p <- predict(mock_fit, se.fit = TRUE)
cat("  fit:", head(p$fit, 3), "\n")
cat("  se: ", head(p$se.fit, 3), "\n")

cat("\n--- summary(glm sub-object) ---\n")
summary(mock_fit$glm)

cat("\n--- confint(glm sub-object) ---\n")
print(confint(mock_fit$glm))

cat("\n--- tidy ---\n")
print(tidy.bayesreml(mock_fit))

cat("\n--- tidy (random effects) ---\n")
print(tidy.bayesreml(mock_fit, effects = "random"))

cat("\n--- glance ---\n")
print(glance.bayesreml(mock_fit))

cat("\n--- augment (first 5 rows) ---\n")
print(head(augment.bayesreml(mock_fit), 5))

cat("\n--- plot: residuals ---\n")
plot(mock_fit, type = "residuals")

cat("\n--- plot: effects ---\n")
plot(mock_fit, type = "effects")

cat("\n--- plot: variance ---\n")
plot(mock_fit, type = "variance")

cat("\n--- plot: blups ---\n")
plot(mock_fit, type = "blups")

cat("\n--- plot: pp_check ---\n")
plot(mock_fit, type = "pp_check")

cat("\n--- anova (model comparison) ---\n")
anova(mock_fit, mock_fit)


## ============================================================================
##  SECTION 8: SUMMARY TABLE
## ============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION 8: Model Catalogue Summary\n")
cat(strrep("#", 70), "\n\n")

models <- data.frame(
  Section  = c(rep("Fixed", 8), rep("Random", 13), rep("Rcov", 2),
               rep("Family", 7), rep("Complex", 9)),
  ID       = c(paste0("1", letters[1:8]),
               paste0("2", letters[1:13]),
               paste0("3", letters[1:2]),
               paste0("4", letters[1:7]),
               paste0("5", letters[1:9])),
  Model    = c(
    # Fixed
    "yield ~ 1", "yield ~ 0 + env", "yield ~ env",
    "yield ~ x_cov", "yield ~ env + x_cov", "yield ~ env:geno",
    "yield ~ x_cov + I(x_cov^2)", "yield ~ x_cov:x2",
    # Random
    "~ geno (IID)", "~ geno + rep (crossed)", "~ env:block (nested)",
    "~ geno:env:rep (3-way)", "~ vm(geno, G) (GBLUP)",
    "~ ped(geno, A) (pedigree)", "~ at(env):geno (diagonal GxE)",
    "~ us(env):id(geno) (unstructured)", "~ fa(env,2):id(geno) (FA k=2)",
    "~ fa(env,1):id(geno) (FA k=1)", "~ ar1(row):id(col) (spatial)",
    "~ spl(x_cov) (P-spline)", "~ vm(geno,G):id(env) (genomic GxE)",
    # Rcov
    "~ units (homogeneous)", "~ at(env):units (heterogeneous)",
    # Families
    "gaussian/identity", "binomial/logit", "binomial/probit",
    "poisson/log", "negative_binomial/log", "gamma/log", "beta/logit",
    # Complex
    "Multi-component", "GBLUP + het. resid", "Diag GxE + spline + nested",
    "vm + at GxE", "Pedigree + crossed + priors", "Binomial GLMM",
    "Poisson GLMM + nested", "Weighted Gaussian", "Spatial + genotype"
  ),
  stringsAsFactors = FALSE
)

print(models, right = FALSE, row.names = FALSE)

cat("\n\nTotal models demonstrated:", nrow(models), "\n")
cat("\n=== All models gallery complete ===\n")
