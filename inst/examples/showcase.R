## ============================================================================
##  bayesreml showcase — selected models demonstrating the package
##
##  This script runs several bayesreml models using the bundled met_example
##  dataset and demonstrates the full analysis workflow.
##
##  Prerequisites:
##    library(bayesreml)
##    greta::install_greta_deps()  # one-time TensorFlow setup
## ============================================================================

library(bayesreml)
data(met_example)
dat    <- met_example$dat
G_mat  <- met_example$G_mat
A_mat  <- met_example$A_mat

cat("Dataset:", nrow(dat), "observations,",
    met_example$n_geno, "genotypes x",
    met_example$n_env, "environments\n\n")

## ── 1. Baseline: Intercept-only ─────────────────────────────────────────────

cat("=== Model 1: Intercept-only ===\n")
fit1 <- bayesreml(
  fixed     = yield ~ 1,
  data      = dat,
  n_samples = 500, warmup = 200, chains = 2
)
summary(fit1)

## ── 2. Fixed effects + random intercept ──────────────────────────────────────

cat("\n=== Model 2: Fixed env + random geno ===\n")
fit2 <- bayesreml(
  fixed     = yield ~ env,
  random    = ~ geno,
  data      = dat,
  n_samples = 500, warmup = 200, chains = 2
)
summary(fit2)
cat("\nBLUPs:\n")
print(round(fit2$extras$blups$u_geno, 3))

## ── 3. GBLUP with genomic relationship matrix ────────────────────────────────

cat("\n=== Model 3: GBLUP ===\n")
fit3 <- bayesreml(
  fixed          = yield ~ env,
  random         = ~ vm(geno, Gmat),
  data           = dat,
  known_matrices = list(Gmat = G_mat),
  n_samples      = 500, warmup = 200, chains = 2
)
summary(fit3)

## ── 4. GxE with diagonal structure ───────────────────────────────────────────

cat("\n=== Model 4: at(env):geno (diagonal GxE) ===\n")
fit4 <- bayesreml(
  fixed     = yield ~ env,
  random    = ~ at(env):geno,
  data      = dat,
  n_samples = 500, warmup = 200, chains = 2
)
summary(fit4)

## ── 5. Heterogeneous residual variance ───────────────────────────────────────

cat("\n=== Model 5: Heterogeneous residual ===\n")
fit5 <- bayesreml(
  fixed     = yield ~ env,
  random    = ~ geno,
  rcov      = ~ at(env):units,
  data      = dat,
  n_samples = 500, warmup = 200, chains = 2
)
summary(fit5)

## ── 6. Model comparison ──────────────────────────────────────────────────────

cat("\n=== Model Comparison ===\n")
anova(fit1, fit2, fit5)

## ── 7. Downstream analysis ──────────────────────────────────────────────────

cat("\n=== Downstream: tidy output ===\n")
print(tidy(fit2))
print(glance(fit2))

## ── 8. View generated code ──────────────────────────────────────────────────

cat("\n=== Generated greta code for Model 2 ===\n")
cat(fit2$extras$code, "\n")

## ── 9. Diagnostic plots ──────────────────────────────────────────────────────

cat("\n=== Diagnostic plots ===\n")
plot(fit2, type = "residuals")
plot(fit2, type = "effects")
plot(fit2, type = "variance")

cat("\n=== Showcase complete ===\n")
