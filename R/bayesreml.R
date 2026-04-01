#' Bayesian Mixed Models with ASReml Syntax
#'
#' Specify mixed models using ASReml formula syntax and estimate them via
#' Bayesian MCMC using greta. Returns a three-part result: a GLM-compatible
#' object for use with standard R packages (emmeans, marginaleffects, etc.),
#' native greta output for Bayesian diagnostics, and extras for secondary
#' analyses.
#'
#' @param fixed Two-sided formula: `response ~ fixed_effects`.
#' @param random One-sided formula: `~ random_terms` using ASReml syntax.
#'   Supports `vm()`, `at()`, `us()`, `fa()`, `ar1()`, `spl()`, `ped()`,
#'   `dsum()`, `id()`, and nested colon terms.
#' @param rcov One-sided formula: `~ residual_structure`. Default `~ units`
#'   (iid residuals). Use `~ at(env):units` for heterogeneous variance.
#' @param data A data.frame containing all variables referenced in the formulas.
#' @param family Character: `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"negative_binomial"`, `"gamma"`, or `"beta"`.
#' @param link Character or NULL: override the default link function
#'   (e.g., `"probit"` for binomial).
#' @param known_matrices Named list of matrices referred to in the random
#'   formula (e.g., `list(Gmat = G_mat, Amat = A_mat)`).
#' @param weights Optional numeric weight vector (length N); sets
#'   Var = sigma^2 / w.
#' @param n_samples Integer: number of posterior samples per chain.
#' @param warmup Integer: number of warmup (burn-in) iterations per chain.
#' @param chains Integer: number of MCMC chains.
#' @param prior_fixed_sd Numeric: SD for fixed-effect normal priors.
#' @param prior_vc_sd Numeric: hyperparameter for lognormal variance-component
#'   priors.
#' @param verbose Logical: print generated greta code to console.
#' @param mcmc_verbose Logical: show MCMC progress bar from greta.
#' @param return_code Logical: if TRUE, return the generated code string
#'   without fitting the model.
#'
#' @return An object of class `"bayesreml"` — a list with three components:
#' \describe{
#'   \item{`$glm`}{A GLM-compatible object (class `c("bayesreml_glm", "glm",
#'     "lm")`) with posterior mean coefficients, vcov, residuals, fitted values,
#'     etc. Works with `summary()`, `emmeans()`, `marginaleffects()`,
#'     `effectsize()`.}
#'   \item{`$greta`}{Native greta output: `model`, `draws` (mcmc.list),
#'     `greta_arrays`, and `env`. Use with `bayesplot`, `greta::calculate()`,
#'     `posterior::as_draws()`.}
#'   \item{`$extras`}{Additional outputs: posterior `summary`, `convergence`
#'     diagnostics, `variance_comps`, `blups`, `predictions`, generated `code`,
#'     `param_names`, `parse_info`, `call_info`, `run_time`, `model_info`.}
#' }
#'
#' If `return_code = TRUE`, returns a character string of greta code instead.
#'
#' @examples
#' \dontrun{
#' # Simple random intercept model
#' fit <- bayesreml(
#'   fixed  = yield ~ env,
#'   random = ~ geno,
#'   data   = met_example,
#'   n_samples = 1000, warmup = 500, chains = 4
#' )
#'
#' # Inspect results
#' summary(fit)
#' coef(fit)
#' confint(fit)
#'
#' # Use with emmeans
#' library(emmeans)
#' emmeans(fit, ~ env)
#'
#' # Use with marginaleffects
#' library(marginaleffects)
#' avg_slopes(fit)
#'
#' # Bayesian diagnostics
#' library(bayesplot)
#' mcmc_trace(fit$greta$draws)
#'
#' # View generated greta code
#' cat(fit$extras$code)
#' }
#'
#' @importFrom stats terms model.frame model.matrix model.response
#'   formula family gaussian binomial poisson Gamma
#'   fitted residuals predict confint coef vcov logLik
#'   nobs update anova na.omit var cov quantile median
#'   qnorm pnorm dnorm setNames df.residual deviance
#'   .getXlevels dbinom density dpois lowess
#'   printCoefmat qqline qqnorm
#' @importFrom graphics abline axis hist legend lines par segments
#' @importFrom methods is
#' @importFrom coda effectiveSize gelman.diag
#' @importFrom splines bs
#' @export
bayesreml <- function(
    fixed,
    random           = NULL,
    rcov             = NULL,
    data,
    family           = "gaussian",
    link             = NULL,
    known_matrices   = list(),
    weights          = NULL,
    n_samples        = 1000,
    warmup           = 500,
    chains           = 4,
    prior_fixed_sd   = 10,
    prior_vc_sd      = 1,
    verbose          = TRUE,
    mcmc_verbose     = TRUE,
    return_code      = FALSE
) {
  if (!requireNamespace("greta", quietly = TRUE))
    stop("Package 'greta' is required. Install with:\n",
         "  install.packages('greta')\n",
         "  greta::install_greta_deps()")

  the_call <- match.call()

  # 1. Parse formulas
  fixed_info   <- .parse_fixed(fixed, data)
  random_terms <- if (!is.null(random)) .parse_formula(random, data) else list()
  rcov_terms   <- if (!is.null(rcov))   .parse_formula(rcov,   data) else
                  list(list(type = "units"))

  fam_link <- .resolve_family(family, link)

  # 2. Build evaluation environment
  #    Parent must include greta's namespace so bare calls like normal(),

  #    as_data(), distribution() resolve correctly.
  greta_ns <- asNamespace("greta")
  ev <- new.env(parent = greta_ns)
  .setup_env(ev, fixed_info, random_terms, rcov_terms,
             data, known_matrices, weights)

  # 3. Generate code
  ctx <- list(
    params    = character(0),
    predictor = character(0),
    code      = character(0),
    fam_link  = fam_link,
    prior_vc  = prior_vc_sd,
    prior_fx  = prior_fixed_sd
  )

  ctx <- .code_fixed(ctx, fixed_info)
  ctx <- .code_random(ctx, random_terms, data, known_matrices)
  ctx <- .code_rcov(ctx, rcov_terms, data)
  ctx <- .code_predictor(ctx, fixed_info)
  ctx <- .code_likelihood(ctx, fixed_info, rcov_terms, data, weights)
  ctx <- .code_model(ctx, n_samples, warmup, chains, mcmc_verbose)

  code_str <- paste(ctx$code, collapse = "\n")

  if (verbose) {
    cat("\n-- bayesreml: generated greta code ",
        paste(rep("-", 40), collapse = ""), "\n", sep = "")
    cat(code_str, "\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
  }

  if (return_code) return(invisible(code_str))

  # 4. Evaluate
  t0 <- proc.time()
  tryCatch(
    eval(parse(text = code_str), envir = ev),
    error = function(e)
      stop("bayesreml evaluation error: ", conditionMessage(e))
  )
  elapsed <- unname((proc.time() - t0)["elapsed"])

  # 5. Post-processing
  # Extract greta arrays
  all_obj_names <- ls(ev)
  is_greta <- vapply(all_obj_names, function(nm) {
    tryCatch(inherits(get(nm, envir = ev, inherits = FALSE), "greta_array"),
             error = function(e) FALSE)
  }, logical(1))
  greta_arrays <- if (any(is_greta)) {
    nms <- all_obj_names[is_greta]
    setNames(lapply(nms, function(nm) get(nm, envir = ev, inherits = FALSE)),
             nms)
  } else {
    list()
  }

  # Posterior summary
  post_summary <- tryCatch(summary(ev$atg_draws), error = function(e) NULL)

  # Convergence diagnostics
  n_eff  <- tryCatch(coda::effectiveSize(ev$atg_draws),
                     error = function(e) NULL)
  gelman <- if (chains >= 2 && length(ev$atg_draws) >= 2) {
    tryCatch(coda::gelman.diag(ev$atg_draws, multivariate = FALSE),
             error = function(e) NULL)
  } else NULL
  conv_diag <- list(n_eff = n_eff, gelman = gelman)

  # 6. Build the three-part output

  # Part 1: GLM-compatible output
  glm_obj <- .build_glm(
    fixed_info  = fixed_info,
    data        = data,
    draws       = ev$atg_draws,
    ev          = ev,
    fam_link    = fam_link,
    the_call    = the_call,
    fixed       = fixed
  )

  # Part 2: Native greta output
  greta_out <- structure(
    list(
      model        = ev$atg_model,
      draws        = ev$atg_draws,
      greta_arrays = greta_arrays,
      env          = ev
    ),
    class = "bayesreml_greta"
  )

  # Part 3: Extras
  # Variance components table
  vc_table <- .build_variance_comps(post_summary, ctx$params)
  # BLUPs
  blups <- .build_blups(post_summary, random_terms)
  # Predictions
  preds <- .build_predictions(data, fixed_info, post_summary, ev, fam_link)

  extras <- structure(
    list(
      summary       = post_summary,
      convergence   = conv_diag,
      variance_comps = vc_table,
      blups         = blups,
      predictions   = preds,
      code          = code_str,
      param_names   = unique(ctx$params),
      parse_info    = list(
        fixed  = fixed_info,
        random = random_terms,
        rcov   = rcov_terms,
        family = fam_link
      ),
      call_info     = list(
        fixed          = fixed,
        random         = random,
        rcov           = rcov,
        data_name      = deparse(substitute(data)),
        family         = family,
        link           = link,
        known_matrices = known_matrices,
        weights        = weights,
        n_samples      = n_samples,
        warmup         = warmup,
        chains         = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd    = prior_vc_sd
      ),
      run_time      = elapsed,
      model_info    = list(
        n_obs     = nrow(data),
        n_fixed   = length(fixed_info$terms) + fixed_info$intercept,
        n_random  = length(random_terms),
        n_params  = length(unique(ctx$params)),
        family    = fam_link$family,
        link      = fam_link$link
      )
    ),
    class = "bayesreml_extras"
  )

  # 7. Assemble final object
  result <- structure(
    list(
      glm    = glm_obj,
      greta  = greta_out,
      extras = extras
    ),
    class = "bayesreml"
  )

  invisible(result)
}
