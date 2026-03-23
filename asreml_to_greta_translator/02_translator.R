## ============================================================================
##  source("02_translator.R")
##
##  Translates ASReml model specifications to efficient greta implementations.
##  v2: enriched return value — richer output for downstream editing / analysis.
##
##  Usage:
##    fit <- asreml_to_greta(
##      fixed  = yield ~ env + x_cov,
##      random = ~ vm(geno, Gmat) + at(env):geno,
##      rcov   = ~ at(env):units,
##      data   = dat,
##      known_matrices = list(Gmat = G_mat)
##    )
##
##  Returns an object of class "atg_fit" — a list with:
##
##  Core greta objects (ready for further use):
##    $model        greta_model — re-run with mcmc($model, n_samples=...)
##    $draws        mcmc.list   — raw posterior samples, works with all coda fns
##    $greta_arrays named list of ALL greta arrays in the eval env, keyed by name
##                  use with calculate() for derived quantities & posterior pred.
##
##  Posterior summaries:
##    $summary      output of summary($draws)  — means, SDs, quantiles
##    $convergence  list(n_eff, gelman)         — ESS and Gelman-Rubin R-hat
##
##  Metadata for editing / reproducibility:
##    $code         complete generated greta R code string — paste and edit
##    $param_names  character vector of parameters passed to model()
##    $parse_info   list(fixed, random, rcov, family) — parsed formula structure
##    $call_info    list of all original arguments to asreml_to_greta()
##    $run_time     elapsed seconds for mcmc()
##
##  Environment:
##    $env          the evaluation environment containing all objects created
##                  by the generated code (data arrays, design matrices, etc.)
##
##  S3 methods:  print(), summary(), coef()
##
##  Unchanged from v1:  all .atg_* internal parsers and code generators.
## ============================================================================

## ── Main exported function ──────────────────────────────────────────────────

#' Translate an ASReml model specification to a fitted greta model.
#'
#' @param fixed     Formula: response ~ fixed_effects  (standard R formula)
#' @param random    Formula: ~ random_terms  (ASReml syntax)
#' @param rcov      Formula: ~ residual_structure  (ASReml syntax)
#' @param data      data.frame
#' @param family    "gaussian" | "binomial" | "poisson" | "negative_binomial" |
#'                  "gamma" | "beta"
#' @param link      Override default link (e.g. "probit" for binomial)
#' @param known_matrices  Named list of matrices referred to in random formula
#'                        e.g. list(Gmat = G_mat, Amat = A_mat)
#' @param weights   Optional numeric weight vector (length N); Var = sigma^2/w
#' @param n_samples MCMC samples
#' @param warmup    MCMC warmup
#' @param chains    Number of MCMC chains
#' @param prior_fixed_sd  SD for fixed-effect normal priors (default 10)
#' @param prior_vc_sd     Hyperparameter for lognormal variance-component priors
#' @param verbose      Print generated code to console
#' @param mcmc_verbose Pass verbose = TRUE/FALSE to greta::mcmc() (default TRUE — shows MCMC progress)
#' @param return_code  If TRUE, return code string without fitting
#'
#' @return An "atg_fit" object (list).  Key components:
#'   \item{model}{greta_model — re-run MCMC with mcmc(fit$model, ...)}
#'   \item{draws}{mcmc.list  — posterior samples (coda-compatible)}
#'   \item{greta_arrays}{named list of all greta arrays for calculate()}
#'   \item{summary}{summary(draws) — posterior statistics}
#'   \item{convergence}{list(n_eff, gelman) — diagnostics}
#'   \item{code}{generated greta R code string}
#'   \item{param_names}{parameters included in model()}
#'   \item{parse_info}{parsed formula structure}
#'   \item{call_info}{original call arguments}
#'   \item{run_time}{elapsed seconds for mcmc()}
#'   \item{env}{evaluation environment}
asreml_to_greta <- function(
    fixed,
    random           = NULL,
    rcov             = NULL,
    data,
    family           = "gaussian",
    link             = NULL,
    known_matrices   = list(),
    weights          = NULL,
    n_samples        = 500,
    warmup           = 200,
    chains           = 2,
    prior_fixed_sd   = 10,
    prior_vc_sd      = 1,
    verbose          = TRUE,
    mcmc_verbose     = TRUE,
    return_code      = FALSE
) {
  if (!requireNamespace("greta", quietly = TRUE))
    stop("Package 'greta' is required.")

  # ── 1. Parse formulas ──────────────────────────────────────────────────────
  fixed_info    <- .atg_parse_fixed(fixed, data)
  random_terms  <- if (!is.null(random)) .atg_parse_formula(random, data) else list()
  rcov_terms    <- if (!is.null(rcov))   .atg_parse_formula(rcov,   data) else
                   list(list(type = "units"))

  fam_link <- .atg_resolve_family(family, link)

  # ── 2. Build evaluation environment ────────────────────────────────────────
  ev <- new.env(parent = parent.frame())
  .atg_setup_env(ev, fixed_info, random_terms, rcov_terms,
                 data, known_matrices, weights)

  # ── 3. Generate code ───────────────────────────────────────────────────────
  ctx <- list(
    params     = character(0),   # names to include in model()
    predictor  = character(0),   # additive terms for mu_i
    code       = character(0),   # accumulated code lines
    fam_link   = fam_link,
    prior_vc   = prior_vc_sd,
    prior_fx   = prior_fixed_sd
  )

  ctx <- .atg_code_fixed(ctx, fixed_info)
  ctx <- .atg_code_random(ctx, random_terms, data, known_matrices)
  ctx <- .atg_code_rcov(ctx, rcov_terms, data)
  ctx <- .atg_code_predictor(ctx, fixed_info)
  ctx <- .atg_code_likelihood(ctx, fixed_info, rcov_terms, data, weights)
  ctx <- .atg_code_model(ctx, n_samples, warmup, chains, mcmc_verbose)

  code_str <- paste(ctx$code, collapse = "\n")

  if (verbose) {
    cat("\n── asreml_to_greta: generated greta code ──────────────────────────\n")
    cat(code_str, "\n")
    cat("───────────────────────────────────────────────────────────────────\n\n")
  }

  if (return_code) return(invisible(code_str))

  # ── 4. Evaluate ─────────────────────────────────────────────────────────────
  t0 <- proc.time()
  tryCatch(
    eval(parse(text = code_str), envir = ev),
    error = function(e) stop("asreml_to_greta evaluation error: ", conditionMessage(e))
  )
  elapsed <- unname((proc.time() - t0)["elapsed"])

  # ── 5. Post-processing: extract greta arrays, summaries, diagnostics ────────

  ## Named list of ALL greta_array objects created in ev
  ## (covers params + intermediate quantities like u_*, g_mat_*, f_spl_* etc.)
  ## Use with: calculate(fit$greta_arrays$u_geno, values = fit$draws)
  all_obj_names <- ls(ev)
  is_greta <- vapply(all_obj_names, function(nm) {
    tryCatch(inherits(get(nm, envir = ev, inherits = FALSE), "greta_array"),
             error = function(e) FALSE)
  }, logical(1))
  greta_arrays <- if (any(is_greta)) {
    nms <- all_obj_names[is_greta]
    setNames(lapply(nms, function(nm) get(nm, envir = ev, inherits = FALSE)), nms)
  } else {
    list()
  }

  ## Posterior summary (coda summary.mcmc.list)
  post_summary <- tryCatch(summary(ev$atg_draws), error = function(e) NULL)

  ## Convergence diagnostics
  n_eff  <- tryCatch(coda::effectiveSize(ev$atg_draws), error = function(e) NULL)
  gelman <- if (chains >= 2 && length(ev$atg_draws) >= 2) {
    tryCatch(coda::gelman.diag(ev$atg_draws, multivariate = FALSE),
             error = function(e) NULL)
  } else {
    NULL
  }
  conv_diag <- list(n_eff = n_eff, gelman = gelman)

  ## ── 6. Assemble return value ────────────────────────────────────────────────
  result <- structure(
    list(
      ## ── Core greta objects ─────────────────────────────────────────────────
      ## Re-run MCMC:  mcmc(fit$model, n_samples = 2000, warmup = 500, chains = 4)
      ## Compute derived quantities:  calculate(fit$greta_arrays$u_geno, values = fit$draws)
      model        = ev$atg_model,
      draws        = ev$atg_draws,
      greta_arrays = greta_arrays,

      ## ── Posterior summaries ────────────────────────────────────────────────
      ## Posterior statistics:       fit$summary$statistics  (mean, SD, naive SE, TS-SE)
      ## Posterior quantiles:        fit$summary$quantiles   (2.5%, 25%, 50%, 75%, 97.5%)
      ## Effective sample sizes:     fit$convergence$n_eff
      ## Gelman-Rubin R-hat:         fit$convergence$gelman$psrf  (if chains >= 2)
      summary      = post_summary,
      convergence  = conv_diag,

      ## ── Generated code ─────────────────────────────────────────────────────
      ## View / edit:  cat(fit$code)
      ## Source directly:  eval(parse(text = fit$code), envir = new.env(...))
      code         = code_str,

      ## ── Metadata ───────────────────────────────────────────────────────────
      param_names  = unique(ctx$params),

      ## parse_info: programmatic access to what was translated
      ##   fit$parse_info$fixed   — response, intercept flag, term list (type/var/n_levels)
      ##   fit$parse_info$random  — list of random term descriptors
      ##   fit$parse_info$rcov    — list of rcov term descriptors
      ##   fit$parse_info$family  — list(family, link)
      parse_info   = list(
        fixed   = fixed_info,
        random  = random_terms,
        rcov    = rcov_terms,
        family  = fam_link
      ),

      ## call_info: original arguments — for reproducibility and re-running
      ##   Refit with more samples:
      ##     do.call(asreml_to_greta, modifyList(fit$call_info, list(n_samples=2000)))
      call_info    = list(
        fixed          = fixed,
        random         = random,
        rcov           = rcov,
        data           = data,
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

      run_time     = elapsed,   # elapsed seconds for mcmc()

      ## ── Full evaluation environment ────────────────────────────────────────
      ## Access any object by name:  fit$env$sigma_geno,  fit$env$u_geno, etc.
      ## List all objects:  ls(fit$env)
      env          = ev
    ),
    class = "atg_fit"
  )

  invisible(result)
}

## ── S3 methods ───────────────────────────────────────────────────────────────

#' Print an atg_fit object
#' @export
print.atg_fit <- function(x, ...) {
  ci <- x$call_info
  pi <- x$parse_info

  cat("ASReml-to-greta fit  [atg_fit]\n")
  cat(strrep("-", 50), "\n")

  # Model specification
  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) cat("  Random :", deparse(ci$random), "\n")
  if (!is.null(ci$rcov))   cat("  Rcov   :", deparse(ci$rcov),   "\n")
  cat("  Family :", pi$family$family,
      "(", pi$family$link, "link )\n")

  # MCMC info
  nch <- length(x$draws)
  ns  <- nrow(x$draws[[1]])
  cat("  MCMC   :", nch, "chain(s) x", ns, "samples",
      "(warmup =", ci$warmup, ") —",
      round(x$run_time, 1), "sec\n")

  # Parameters
  np <- length(x$param_names)
  nall <- length(x$greta_arrays)
  cat("  Params :", np, "in model(); ",
      nall - np, "intermediate greta arrays available via $greta_arrays\n")

  # Quick convergence
  if (!is.null(x$convergence$gelman)) {
    rhat <- x$convergence$gelman$psrf[, "Point est."]
    max_rhat <- max(rhat, na.rm = TRUE)
    flag <- if (max_rhat < 1.05) " [OK]" else if (max_rhat < 1.10) " [borderline]" else " [!]"
    cat("  Max R-hat:", round(max_rhat, 3), flag, "\n")
  }
  min_eff <- if (!is.null(x$convergence$n_eff))
               min(x$convergence$n_eff, na.rm = TRUE) else NA
  if (!is.na(min_eff))
    cat("  Min n_eff:", round(min_eff, 0), "\n")

  cat(strrep("-", 50), "\n")
  cat("  $model        — re-run:  mcmc($model, n_samples=...)\n")
  cat("  $draws        — mcmc.list: use with coda, bayesplot, etc.\n")
  cat("  $greta_arrays — for calculate() on derived quantities\n")
  cat("  $code         — generated greta code (paste & edit)\n")
  cat("  Use summary(fit) for full posterior table.\n")

  invisible(x)
}

#' Summarise posterior of an atg_fit object
#' @export
summary.atg_fit <- function(object, ...) {
  ci <- object$call_info

  cat("Posterior summary — atg_fit\n")
  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) cat("  Random :", deparse(ci$random), "\n")
  cat("  Family :", object$parse_info$family$family, "/",
      object$parse_info$family$link, "\n\n")

  cat("── Parameter statistics ──────────────────────────────────────\n")
  print(object$summary$statistics)

  cat("\n── Posterior quantiles ───────────────────────────────────────\n")
  print(object$summary$quantiles)

  cat("\n── Effective sample sizes ────────────────────────────────────\n")
  if (!is.null(object$convergence$n_eff)) {
    print(round(object$convergence$n_eff, 0))
  } else {
    cat("  (not available)\n")
  }

  if (!is.null(object$convergence$gelman)) {
    cat("\n── Gelman-Rubin R-hat ────────────────────────────────────────\n")
    print(round(object$convergence$gelman$psrf, 3))
    mpsrf <- object$convergence$gelman$mpsrf
    if (!is.null(mpsrf))
      cat("  Multivariate PSRF:", round(mpsrf, 3), "\n")
  }

  invisible(object$summary)
}

#' Extract posterior means of model parameters from an atg_fit object
#'
#' Returns a named numeric vector of posterior means for all parameters
#' that were included in model() (i.e. fit$param_names).
#' @export
coef.atg_fit <- function(object, ...) {
  stats <- object$summary$statistics
  if (is.null(stats)) {
    warning("No summary statistics available.")
    return(NULL)
  }
  # statistics has rownames = parameter names (may include [i] indices)
  means <- stats[, "Mean"]
  means
}

## ── Internal: formula parsers ────────────────────────────────────────────────
## All code below is identical to v1 — no changes.

# Parse a standard R fixed formula → list(response, intercept, terms)
.atg_parse_fixed <- function(fixed, data) {
  if (length(fixed) < 3) stop("fixed must be a two-sided formula: response ~ ...")
  response  <- deparse(fixed[[2]])
  if (!response %in% names(data))
    stop("Response variable '", response, "' not found in data.")
  rhs       <- fixed[[3]]

  # Use terms() on a version with 1 for intercept detection
  tt <- tryCatch(terms(fixed, data = data), error = function(e) terms(fixed))
  intercept   <- attr(tt, "intercept") == 1
  term_labels <- attr(tt, "term.labels")

  term_info <- lapply(term_labels, function(lbl) .atg_classify_fixed_term(lbl, data))

  list(response = response, intercept = intercept, terms = term_info)
}

# Classify a single fixed-formula term label
.atg_classify_fixed_term <- function(lbl, data) {
  if (grepl("^I\\(", lbl)) {
    return(list(type = "expression", label = lbl))
  }
  if (grepl(":", lbl, fixed = TRUE)) {
    parts <- strsplit(lbl, ":")[[1]]
    all_factors <- all(vapply(parts, function(p) {
      p %in% names(data) && (is.factor(data[[p]]) || is.character(data[[p]]))
    }, logical(1)))
    if (all_factors)
      return(list(type = "factor_interaction", vars = parts, label = lbl))
    return(list(type = "interaction", vars = parts, label = lbl))
  }
  if (lbl %in% names(data)) {
    col <- data[[lbl]]
    if (is.factor(col) || is.character(col)) {
      f <- factor(col)
      return(list(type = "factor", var = lbl,
                  n_levels = nlevels(f), levels = levels(f), label = lbl))
    } else {
      return(list(type = "continuous", var = lbl, label = lbl))
    }
  }
  list(type = "expression", label = lbl)
}

# Parse an ASReml one-sided formula (random or rcov) via AST walking
.atg_parse_formula <- function(formula, data) {
  rhs   <- if (length(formula) == 3) formula[[3]] else formula[[2]]
  terms <- .atg_walk(rhs)
  lapply(terms, .atg_enrich, data = data)
}

# ── Recursive AST walker ─────────────────────────────────────────────────────
.atg_walk <- function(expr) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm == "units") return(list(list(type = "units")))
    return(list(list(type = "simple", var = nm)))
  }
  if (!is.call(expr)) return(list())

  fn <- as.character(expr[[1]])

  # Binary operators
  if (fn == "+") return(c(.atg_walk(expr[[2]]), .atg_walk(expr[[3]])))

  if (fn == ":") {
    left  <- .atg_walk(expr[[2]])
    right <- .atg_walk(expr[[3]])
    l <- if (length(left))  left[[1]]  else list(type = "unknown")
    r <- if (length(right)) right[[1]] else list(type = "unknown")
    return(.atg_classify_interaction(l, r))
  }

  # ASReml-specific functions
  switch(fn,
    "vm"   = list(list(type = "vm",  var = .dep(expr, 2),
                       mat  = .dep(expr, 3))),
    "ide"  = ,
    "id"   = list(list(type = "ide", var = .dep(expr, 2))),
    "ped"  = list(list(type = "ped", var = .dep(expr, 2),
                       mat  = .dep(expr, 3))),
    "fa"   = list(list(type = "fa",  var = .dep(expr, 2),
                       k    = as.integer(expr[[3]]))),
    "us"   = list(list(type = "us",  var = .dep(expr, 2))),
    "at"   = list(list(type = "at",  var = .dep(expr, 2),
                       level = if (length(expr) > 2) .dep(expr, 3) else NULL)),
    "ar1"  = list(list(type = "ar1", var = .dep(expr, 2))),
    "ar2"  = list(list(type = "ar2", var = .dep(expr, 2))),
    "cor"  = list(list(type = "cor", var = .dep(expr, 2))),
    "str"  = list(list(type = "str", formula = expr[[2]],
                       structure = if (length(expr) > 2) expr[[3]] else NULL)),
    "spl"  = list(list(type = "spline", var = .dep(expr, 2))),
    "dsum" = {
      # dsum(~ units|grp) — heterogeneous residual by grouping variable
      # Equivalent to at(grp):units
      # The argument is a language call `~(units|grp)`, not a formula object
      arg <- expr[[2]]
      grp <- tryCatch({
        # Navigate past the tilde to reach `units|grp`
        body_expr <- if (is.call(arg) && as.character(arg[[1]]) == "~")
          arg[[length(arg)]]
        else
          arg
        if (is.call(body_expr) && as.character(body_expr[[1]]) == "|")
          deparse(body_expr[[3]])
        else
          NULL
      }, error = function(e) NULL)
      if (!is.null(grp))
        list(list(type = "at_units", var = grp, level = NULL))
      else
        list(list(type = "units"))
    },
    "lin"  = list(list(type = "continuous", var = .dep(expr, 2))),
    "pol"  = list(list(type = "polynomial", var = .dep(expr, 2),
                       degree = if (length(expr) > 2) as.integer(expr[[3]]) else 2L)),
    # Default: treat as a simple named variable
    list(list(type = "simple", var = deparse(expr)))
  )
}

.dep <- function(expr, i) if (length(expr) >= i) deparse(expr[[i]]) else NA_character_

# Classify what an `a:b` interaction means in ASReml context
.atg_classify_interaction <- function(l, r) {
  lt <- l$type; rt <- r$type
  # FA(k) GxE:  fa(env,k) : id(geno)
  if (lt == "fa" && rt %in% c("id","ide","simple")) {
    if (!is.null(l$k) && l$k < 1L)
      stop("fa() requires k >= 1; got k = ", l$k, ".")
    return(list(list(type="fa_gxe",  outer=l$var, k=l$k,
                     inner=r$var)))
  }
  # US GxE:     us(env) : id(geno)
  if (lt == "us" && rt %in% c("id","ide","simple"))
    return(list(list(type="us_gxe",  outer=l$var, inner=r$var)))
  # at(env):units  – heterogeneous residual
  if (lt == "at" && rt == "units")
    return(list(list(type="at_units", var=l$var, level=l$level)))
  # at(env):geno   – DIAG(env) x I(geno)
  if (lt == "at" && rt %in% c("simple","ide","id"))
    return(list(list(type="at_simple", outer=l$var, level=l$level,
                     inner=r$var)))
  # AR1 spatial
  if (lt == "ar1" && rt %in% c("ar1","id","ide","simple","units"))
    return(list(list(type="ar1_spatial", row_var=l$var, col_var=r$var,
                     col_ar1 = (rt == "ar1"))))
  # vm(geno,G) : id(env)  – structured GxE with known G
  if (lt == "vm" && rt %in% c("id","ide","simple"))
    return(list(list(type="vm_gxe", inner=l$var, mat=l$mat, outer=r$var)))
  # simple:simple  – nested random effect
  if (lt %in% c("simple","ide","id") && rt %in% c("simple","ide","id"))
    return(list(list(type="nested", outer=l$var, inner=r$var)))
  # nested:simple  – three-way combination (e.g. rep:loc:yearf)
  if (lt == "nested" && rt %in% c("simple","ide","id"))
    return(list(list(type="combo", vars=c(l$outer, l$inner, r$var))))
  # combo:simple   – four-way+ combination (e.g. A:B:C:D)
  if (lt == "combo" && rt %in% c("simple","ide","id"))
    return(list(list(type="combo", vars=c(l$vars, r$var))))
  # fallback
  list(list(type="interaction_generic", left=l, right=r))
}

# Enrich a term descriptor with data-derived information
.atg_enrich <- function(term, data) {
  .add_var_info <- function(term, var_field = "var") {
    vname <- term[[var_field]]
    if (!is.null(vname) && !is.na(vname) && vname %in% names(data)) {
      col <- data[[vname]]
      f   <- factor(col)
      term[[paste0(var_field, "_n")]]      <- nlevels(f)
      term[[paste0(var_field, "_levels")]] <- levels(f)
    }
    term
  }

  switch(term$type,
    "simple" = ,
    "ide"    = ,
    "id"     = .add_var_info(term),
    "vm"     = .add_var_info(term),
    "ped"    = .add_var_info(term),
    "fa"     = { term <- .add_var_info(term); term },
    "fa_gxe" = {
      if (term$outer %in% names(data)) {
        f <- factor(data[[term$outer]])
        term$n_outer <- nlevels(f)
      }
      if (term$inner %in% names(data)) {
        f <- factor(data[[term$inner]])
        term$n_inner <- nlevels(f)
      }
      term
    },
    "us_gxe" = {
      if (term$outer %in% names(data)) term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data)) term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    "at_simple" = {
      if (term$outer %in% names(data)) term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data)) term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    "at_units" = {
      if (term$var %in% names(data)) term$n_var <- nlevels(factor(data[[term$var]]))
      term
    },
    "ar1_spatial" = {
      if (!is.null(term$row_var) && term$row_var %in% names(data))
        term$n_row <- length(unique(data[[term$row_var]]))
      if (!is.null(term$col_var) && term$col_var %in% names(data))
        term$n_col <- length(unique(data[[term$col_var]]))
      term
    },
    "nested" = {
      if (term$outer %in% names(data)) term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data)) term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    term  # default: return unchanged
  )
}

## ── Internal: environment setup ──────────────────────────────────────────────

# Populate ev with all R objects the generated code will reference
.atg_setup_env <- function(ev, fixed_info, random_terms, rcov_terms,
                            data, known_matrices, weights) {
  N <- nrow(data)
  ev$N_atg   <- N
  ev$dat_atg <- data

  # Response
  resp <- fixed_info$response
  ev$y_atg <- data[[resp]]

  # Fixed effect id vectors
  for (term in fixed_info$terms) {
    if (term$type == "factor") {
      f <- factor(data[[term$var]])
      assign(paste0(term$var, "_id"),  as.integer(f),  envir = ev)
      assign(paste0("n_", term$var),   nlevels(f),      envir = ev)
    } else if (term$type == "continuous") {
      if (term$var %in% names(data))
        assign(term$var, as.numeric(data[[term$var]]), envir = ev)
    } else if (term$type == "interaction") {
      for (v in term$vars) {
        if (v %in% names(data)) {
          col <- data[[v]]
          if (is.factor(col) || is.character(col)) {
            f <- factor(col)
            assign(paste0(v, "_id"), as.integer(f), envir = ev)
            assign(paste0("n_", v),  nlevels(f),    envir = ev)
          } else {
            assign(v, as.numeric(col), envir = ev)
          }
        }
      }
    } else if (term$type == "factor_interaction") {
      tag   <- paste(term$vars, collapse = "_x_")
      combo <- do.call(paste, c(lapply(term$vars, function(v) data[[v]]), list(sep = "_")))
      f     <- factor(combo)
      assign(paste0(tag, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", tag),  nlevels(f),    envir = ev)
    }
  }

  # Random and rcov id vectors / level counts / matrices
  all_terms <- c(random_terms, rcov_terms)
  for (term in all_terms) {
    .atg_setup_term_env(ev, term, data, known_matrices)
  }

  # Weights
  if (!is.null(weights)) ev$wt_atg <- weights
}

.atg_setup_term_env <- function(ev, term, data, known_matrices) {
  mk_id <- function(vname) {
    if (vname %in% names(data)) {
      f <- factor(data[[vname]])
      assign(paste0(vname, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", vname),  nlevels(f),    envir = ev)
    } else if (paste0(vname, "_id") %in% names(data)) {
      # Column already stored with _id suffix (pre-computed integer IDs)
      id_col <- data[[paste0(vname, "_id")]]
      assign(paste0(vname, "_id"), id_col,         envir = ev)
      assign(paste0("n_", vname),  max(id_col),    envir = ev)
    } else {
      warning("Variable '", vname, "' not found in data — id setup skipped.")
    }
  }
  mk_num <- function(vname) {
    if (!vname %in% names(data)) return()
    assign(vname, as.numeric(data[[vname]]), envir = ev)
  }
  mk_mat <- function(mat_name) {
    if (is.null(mat_name) || is.na(mat_name)) return()
    if (!mat_name %in% names(known_matrices))
      stop("Matrix '", mat_name, "' not found in known_matrices. ",
           "Provide it via known_matrices = list(", mat_name, " = <your matrix>).")
    assign(mat_name, known_matrices[[mat_name]], envir = ev)
  }

  switch(term$type,
    "simple" = , "ide" = , "id" = mk_id(term$var),
    "vm"     = { mk_id(term$var); mk_mat(term$mat) },
    "ped"    = { mk_id(term$var); mk_mat(term$mat) },
    "fa"     = mk_id(term$var),
    "fa_gxe" = { mk_id(term$outer); mk_id(term$inner) },
    "us_gxe" = { mk_id(term$outer); mk_id(term$inner) },
    "at_simple" = { mk_id(term$outer); mk_id(term$inner) },
    "at_units"  = mk_id(term$var),
    "nested" = {
      tag <- paste0(term$inner, "_in_", term$outer)
      nid <- as.integer(factor(paste(data[[term$outer]], data[[term$inner]], sep = "_")))
      ev[[paste0("nested_id_", tag)]] <- nid
      ev[[paste0("n_",         tag)]] <- length(unique(nid))
    },
    "combo" = {
      tag <- paste(rev(term$vars), collapse = "_in_")
      nid <- as.integer(factor(do.call(paste, c(lapply(term$vars, function(v) data[[v]]), list(sep = "_")))))
      ev[[paste0("combo_id_", tag)]] <- nid
      ev[[paste0("n_",        tag)]] <- length(unique(nid))
    },
    "ar1_spatial" = {
      if (!is.null(term$row_var) && term$row_var %in% names(data)) {
        row_vals <- data[[term$row_var]]
        assign(paste0(term$row_var, "_id"), as.integer(factor(row_vals)), envir = ev)
        assign(paste0("n_", term$row_var),  length(unique(row_vals)),     envir = ev)
      }
      if (!is.null(term$col_var) && term$col_var %in% names(data)) {
        col_vals <- data[[term$col_var]]
        assign(paste0(term$col_var, "_id"), as.integer(factor(col_vals)), envir = ev)
        assign(paste0("n_", term$col_var),  length(unique(col_vals)),     envir = ev)
      }
    }
  )
}

## ── Internal: resolve family/link ────────────────────────────────────────────

.atg_resolve_family <- function(family, link) {
  defaults <- list(
    gaussian          = "identity",
    binomial          = "logit",
    binary            = "logit",
    poisson           = "log",
    negative_binomial = "log",
    negbinom          = "log",
    gamma             = "log",
    beta              = "logit"
  )
  fam <- tolower(family)
  if (!fam %in% names(defaults)) stop("Unsupported family: ", family)
  lnk <- if (!is.null(link)) tolower(link) else defaults[[fam]]
  list(family = fam, link = lnk)
}

## ── Internal: code generators ────────────────────────────────────────────────

# Add a line of code to ctx
.add <- function(ctx, ...) {
  ctx$code <- c(ctx$code, paste0(...))
  ctx
}
.add_param <- function(ctx, ...) {
  ctx$params <- c(ctx$params, c(...))
  ctx
}
.add_pred <- function(ctx, term_str) {
  ctx$predictor <- c(ctx$predictor, term_str)
  ctx
}

# ── Fixed effects ─────────────────────────────────────────────────────────────
.atg_code_fixed <- function(ctx, fixed_info) {
  psd <- ctx$prior_fx
  ctx <- .add(ctx, "# ── Fixed effects ───────────────────────────────────────")

  # Intercept
  if (fixed_info$intercept && length(fixed_info$terms) == 0) {
    ctx <- .add(ctx, "mu_atg <- normal(0, ", psd, ")")
    ctx <- .add_param(ctx, "mu_atg")
    ctx <- .add_pred(ctx, "mu_atg")
  } else if (fixed_info$intercept) {
    ctx <- .add(ctx, "mu_atg <- normal(0, ", psd, ")")
    ctx <- .add_param(ctx, "mu_atg")
    ctx <- .add_pred(ctx, "mu_atg")
  }

  for (term in fixed_info$terms) {
    ctx <- switch(term$type,
      "factor" = {
        tag <- term$var
        if (fixed_info$intercept) {
          # Deviation contrasts: tau[k] with mu_atg already in predictor
          ctx <- .add(ctx, "tau_", tag, " <- normal(0, ", psd, ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("tau_", tag))
          ctx <- .add_pred(ctx, paste0("tau_", tag, "[", tag, "_id]"))
        } else {
          # Means model: alpha[k] absorbs intercept
          ctx <- .add(ctx, "alpha_", tag, " <- normal(0, ", psd, ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("alpha_", tag))
          ctx <- .add_pred(ctx, paste0("alpha_", tag, "[", tag, "_id]"))
        }
        ctx
      },
      "continuous" = {
        tag <- term$var
        ctx <- .add(ctx, "beta_", tag, " <- normal(0, 2)")
        ctx <- .add_param(ctx, paste0("beta_", tag))
        ctx <- .add_pred(ctx, paste0("beta_", tag, " * as_data(", tag, ")"))
        ctx
      },
      "factor_interaction" = {
        tag <- paste(term$vars, collapse = "_x_")
        if (fixed_info$intercept) {
          ctx <- .add(ctx, "tau_", tag, " <- normal(0, ", psd, ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("tau_", tag))
          ctx <- .add_pred(ctx, paste0("tau_", tag, "[", tag, "_id]"))
        } else {
          ctx <- .add(ctx, "alpha_", tag, " <- normal(0, ", psd, ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("alpha_", tag))
          ctx <- .add_pred(ctx, paste0("alpha_", tag, "[", tag, "_id]"))
        }
        ctx
      },
      "interaction" = {
        v1 <- term$vars[1]; v2 <- term$vars[2]
        tag <- paste0(v1, "_x_", v2)
        ctx <- .add(ctx, "beta_", tag, " <- normal(0, 2)")
        ctx <- .add_param(ctx, paste0("beta_", tag))
        ctx <- .add_pred(ctx,
          paste0("beta_", tag, " * as_data(", v1, ") * as_data(", v2, ")"))
        ctx
      },
      "expression" = {
        raw_lbl  <- term$label
        tag      <- gsub("[^A-Za-z0-9_]", "_", raw_lbl)
        # Strip I() wrapper so as_data() receives plain numeric, not AsIs
        inner    <- gsub("^I\\((.*)\\)$", "\\1", trimws(raw_lbl))
        ctx <- .add(ctx, "beta_", tag, " <- normal(0, 2)")
        ctx <- .add(ctx, "x_", tag, " <- as_data(as.numeric(with(dat_atg, ", inner, ")))")
        ctx <- .add_param(ctx, paste0("beta_", tag))
        ctx <- .add_pred(ctx, paste0("beta_", tag, " * x_", tag))
        ctx
      },
      ctx
    )
  }
  ctx
}

# ── Random effects ─────────────────────────────────────────────────────────────
.atg_code_random <- function(ctx, random_terms, data, known_matrices) {
  if (length(random_terms) == 0) return(ctx)
  ctx <- .add(ctx, "# ── Random effects ──────────────────────────────────────")

  for (i in seq_along(random_terms)) {
    term <- random_terms[[i]]
    pvc  <- ctx$prior_vc

    ctx <- switch(term$type,

      # ── iid (simple, ide, id) ────────────────────────────────────────────
      "simple" = , "ide" = , "id" = {
        tag <- term$var
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- ", tag, "_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # ── GBLUP: vm(geno, G) ───────────────────────────────────────────────
      "vm" = {
        tag <- term$var; mat <- term$mat
        ctx <- .add(ctx, "L_G_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- as_data(L_G_", tag, ") %*% (z_", tag, " * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # ── Animal model: ped(animal, A) ─────────────────────────────────────
      "ped" = {
        tag <- term$var; mat <- term$mat
        ctx <- .add(ctx, "L_A_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- as_data(L_A_", tag, ") %*% (z_", tag, " * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # ── DIAG(outer) x I(inner): at(env):geno ────────────────────────────
      "at_simple" = {
        outer <- term$outer; inner <- term$inner
        tag   <- paste0(inner, "_", outer)
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ", dim = n_", outer, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = c(n_", inner, ", n_", outer, "))")
        ctx <- .add(ctx, tag, "_sc  <- sweep(", tag, "_raw, 2, sigma_", tag, ", \"*\")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0(tag, "_sc[cbind(", inner, "_id, ", outer, "_id)]"))
        ctx
      },

      # ── US(outer) x I(inner): us(env):id(geno) ───────────────────────────
      "us_gxe" = {
        outer <- term$outer; inner <- term$inner
        tag   <- paste0(inner, "_", outer, "_us")
        ctx <- .add(ctx, "L_corr_", tag, " <- lkj_correlation(eta = 2, dimension = n_", outer, ")")
        ctx <- .add(ctx, "sd_", tag, " <- greta::lognormal(0, ", pvc, ", dim = n_", outer, ")")
        ctx <- .add(ctx, "L_", tag, " <- sweep(L_corr_", tag, ", 2, sd_", tag, ", \"*\")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = c(n_", inner, ", n_", outer, "))")
        ctx <- .add(ctx, tag, "_mat <- ", tag, "_raw %*% t(L_", tag, ")")
        ctx <- .add_param(ctx, paste0("sd_", tag))
        ctx <- .add_pred(ctx, paste0(tag, "_mat[cbind(", inner, "_id, ", outer, "_id)]"))
        ctx
      },

      # ── FA(k)(outer) x I(inner): fa(env,k):id(geno) ──────────────────────
      "fa_gxe" = {
        outer <- term$outer; inner <- term$inner; k <- term$k
        tag   <- paste0(inner, "_", outer, "_fa", k)
        ctx <- .add(ctx, "Lambda_", tag, " <- normal(0, 1, dim = c(n_", outer, ", ", k, "L))")
        ctx <- .add(ctx, "psi_", tag, "    <- greta::lognormal(0, ", pvc, ", dim = n_", outer, ")")
        ctx <- .add(ctx, "F_", tag, "      <- normal(0, 1, dim = c(n_", inner, ", ", k, "L))")
        ctx <- .add(ctx, "delta_", tag, "  <- normal(0, 1, dim = c(n_", inner, ", n_", outer, "))")
        ctx <- .add(ctx, "delta_sc_", tag, " <- sweep(delta_", tag, ", 2, psi_", tag, ", \"*\")")
        ctx <- .add(ctx, "g_mat_", tag, "  <- F_", tag, " %*% t(Lambda_", tag, ") + delta_sc_", tag)
        ctx <- .add_param(ctx, paste0("Lambda_", tag), paste0("psi_", tag))
        ctx <- .add_pred(ctx, paste0("g_mat_", tag, "[cbind(", inner, "_id, ", outer, "_id)]"))
        ctx
      },

      # ── US(outer) x Struct(inner): vm(geno,G):id(env) ────────────────────
      "vm_gxe" = {
        inner <- term$inner; outer <- term$outer; mat <- term$mat
        tag   <- paste0(inner, "_", outer, "_vm_gxe")
        ctx <- .add(ctx, "L_G_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "L_corr_", tag, " <- lkj_correlation(eta = 2, dimension = n_", outer, ")")
        ctx <- .add(ctx, "sd_", tag, " <- greta::lognormal(0, ", pvc, ", dim = n_", outer, ")")
        ctx <- .add(ctx, "L_", tag, " <- sweep(L_corr_", tag, ", 2, sd_", tag, ", \"*\")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = c(n_", inner, ", n_", outer, "))")
        ctx <- .add(ctx, "uA_", tag, " <- as_data(L_G_", tag, ") %*% z_", tag)
        ctx <- .add(ctx, "g_mat_", tag, " <- uA_", tag, " %*% t(L_", tag, ")")
        ctx <- .add_param(ctx, paste0("sd_", tag))
        ctx <- .add_pred(ctx, paste0("g_mat_", tag, "[cbind(", inner, "_id, ", outer, "_id)]"))
        ctx
      },

      # ── Nested random effect: outer:inner ────────────────────────────────
      "nested" = {
        outer <- term$outer; inner <- term$inner
        tag   <- paste0(inner, "_in_", outer)
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- ", tag, "_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[nested_id_", tag, "]"))
        ctx
      },

      # ── Multi-way combination: A:B:C (three-way or higher) ───────────────
      "combo" = {
        tag <- paste(rev(term$vars), collapse = "_in_")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- ", tag, "_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[combo_id_", tag, "]"))
        ctx
      },

      # ── Spatial AR1 approximation ─────────────────────────────────────────
      "ar1_spatial" = {
        rv <- term$row_var; cv <- term$col_var
        tag <- paste0("sp_", rv, "_", cv)
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_row_raw <- normal(0, 1, dim = n_", rv, ")")
        ctx <- .add(ctx, tag, "_col_raw <- normal(0, 1, dim = n_", cv, ")")
        ctx <- .add(ctx, "u_", tag, "_row <- ", tag, "_row_raw * sigma_", tag)
        ctx <- .add(ctx, "u_", tag, "_col <- ", tag, "_col_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "_row[", rv, "_id] + u_", tag, "_col[", cv, "_id]"))
        ctx
      },

      # ── P-spline: spl(x) ─────────────────────────────────────────────────
      "spline" = {
        vname <- term$var
        tag   <- paste0("spl_", vname)
        ctx <- .add(ctx, "x_std_", tag, " <- as.numeric(scale(dat_atg$", vname, "))")
        ctx <- .add(ctx, "B_", tag, " <- splines::bs(x_std_", tag, ", df=8, degree=3, intercept=FALSE)")
        ctx <- .add(ctx, "B_g_", tag, " <- as_data(B_", tag, ")")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = ncol(B_", tag, "))")
        ctx <- .add(ctx, "f_", tag, " <- B_g_", tag, " %*% (", tag, "_raw * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("f_", tag))
        ctx
      },

      ctx  # default: skip unknown terms
    )
  }
  ctx
}

# ── Residual (rcov) code ───────────────────────────────────────────────────────
.atg_code_rcov <- function(ctx, rcov_terms, data) {
  ctx <- .add(ctx, "# ── Residual ────────────────────────────────────────────")
  pvc <- ctx$prior_vc

  # Find the dominant rcov structure
  has_at_units  <- any(vapply(rcov_terms, function(t) t$type == "at_units",  logical(1)))
  has_ar1       <- any(vapply(rcov_terms, function(t) t$type == "ar1_spatial",logical(1)))

  if (has_at_units) {
    term <- Filter(function(t) t$type == "at_units", rcov_terms)[[1]]
    tag  <- term$var
    ctx  <- .add(ctx, "sigma_e_atg <- greta::lognormal(0, ", pvc, ", dim = n_", tag, ")")
    ctx  <- .add_param(ctx, "sigma_e_atg")
    ctx$resid_expr <- paste0("sigma_e_atg[", tag, "_id]")
  } else {
    ctx  <- .add(ctx, "sigma_e_atg <- greta::lognormal(0, ", pvc, ")")
    ctx  <- .add_param(ctx, "sigma_e_atg")
    ctx$resid_expr <- "sigma_e_atg"
  }

  ctx
}

# ── Linear predictor assembly ──────────────────────────────────────────────────
.atg_code_predictor <- function(ctx, fixed_info) {
  ctx <- .add(ctx, "# ── Linear predictor ────────────────────────────────────")
  if (length(ctx$predictor) == 0) {
    ctx <- .add(ctx, "mu_i_atg <- zeros(N_atg, 1)")
  } else {
    parts <- ctx$predictor
    lhs   <- "mu_i_atg"
    ctx   <- .add(ctx, lhs, " <- ", paste(parts, collapse = " + "))
  }
  ctx
}

# ── Likelihood ────────────────────────────────────────────────────────────────
.atg_code_likelihood <- function(ctx, fixed_info, rcov_terms, data, weights) {
  ctx  <- .add(ctx, "# ── Likelihood ──────────────────────────────────────────")
  resp <- fixed_info$response
  fl   <- ctx$fam_link
  resid_sd <- if (!is.null(ctx$resid_expr)) ctx$resid_expr else "sigma_e_atg"

  ctx <- .add(ctx, "y_atg_obs <- as_data(y_atg)")

  if (!is.null(weights)) {
    ctx <- .add(ctx, "wt_atg_data <- as_data(wt_atg)")
    resid_sd <- paste0(resid_sd, " / sqrt(wt_atg_data)")
  }

  if (fl$family == "gaussian") {
    ctx <- .add(ctx, "distribution(y_atg_obs) <- normal(mu_i_atg, ", resid_sd, ")")

  } else if (fl$family %in% c("binomial","binary")) {
    if (fl$link == "logit") {
      ctx <- .add(ctx, "p_atg <- ilogit(mu_i_atg)")
    } else {
      ctx <- .add(ctx, "p_atg <- iprobit(mu_i_atg)")
    }
    y_vals <- data[[resp]]
    if (all(y_vals %in% c(0, 1))) {
      ctx <- .add(ctx, "distribution(y_atg_obs) <- bernoulli(p_atg)")
    } else {
      ctx <- .add(ctx, "# Binomial with trials — adapt n_trials as needed")
      ctx <- .add(ctx, "distribution(y_atg_obs) <- bernoulli(p_atg)")
    }

  } else if (fl$family == "poisson") {
    ctx <- .add(ctx, "lam_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "distribution(y_atg_obs) <- poisson(lam_atg)")

  } else if (fl$family %in% c("negative_binomial","negbinom")) {
    ctx <- .add(ctx, "r_nb_atg <- greta::lognormal(0, 1)")
    ctx <- .add_param(ctx, "r_nb_atg")
    ctx <- .add(ctx, "mu_nb_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "prob_nb_atg <- r_nb_atg / (r_nb_atg + mu_nb_atg)")
    ctx <- .add(ctx, "distribution(y_atg_obs) <- negative_binomial(r_nb_atg, prob_nb_atg)")

  } else if (fl$family == "gamma") {
    ctx <- .add(ctx, "shape_atg <- greta::lognormal(0, 1)")
    ctx <- .add_param(ctx, "shape_atg")
    ctx <- .add(ctx, "mu_gam_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "rate_atg <- shape_atg / mu_gam_atg")
    ctx <- .add(ctx, "distribution(y_atg_obs) <- greta::gamma(shape_atg, rate_atg)")

  } else if (fl$family == "beta") {
    ctx <- .add(ctx, "phi_b_atg <- greta::lognormal(1, 1)")
    ctx <- .add_param(ctx, "phi_b_atg")
    ctx <- .add(ctx, "mu_b_atg <- ilogit(mu_i_atg)")
    ctx <- .add(ctx, "alpha_b_atg <- mu_b_atg * phi_b_atg")
    ctx <- .add(ctx, "beta_b_atg  <- (1 - mu_b_atg) * phi_b_atg")
    ctx <- .add(ctx, "distribution(y_atg_obs) <- greta::beta(alpha_b_atg, beta_b_atg)")
  }

  ctx
}

# ── model() + mcmc() ──────────────────────────────────────────────────────────
.atg_code_model <- function(ctx, n_samples, warmup, chains, mcmc_verbose = TRUE) {
  ctx <- .add(ctx, "# ── Model and MCMC ───────────────────────────────────────")
  params_str <- paste(unique(ctx$params), collapse = ", ")
  ctx <- .add(ctx, "atg_model <- greta::model(", params_str, ")")
  ctx <- .add(ctx, "atg_draws <- greta::mcmc(atg_model, n_samples = ", n_samples,
              ", warmup = ", warmup, ", chains = ", chains,
              ", verbose = ", mcmc_verbose, ")")
  ctx
}
