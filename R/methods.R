# S3 methods for the bayesreml class

#' Print a bayesreml object
#' @param x A bayesreml object
#' @param ... Additional arguments (ignored)
#' @export
print.bayesreml <- function(x, ...) {
  ci <- x$extras$call_info
  mi <- x$extras$model_info

  cat("Bayesian mixed model  [bayesreml]\n")
  cat(strrep("-", 55), "\n")

  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) cat("  Random :", deparse(ci$random), "\n")
  if (!is.null(ci$rcov) && !identical(deparse(ci$rcov), "~units"))
    cat("  Rcov   :", deparse(ci$rcov), "\n")
  cat("  Family :", mi$family, "(", mi$link, "link )\n")

  nch <- ci$chains
  ns  <- ci$n_samples
  cat("  MCMC   :", nch, "chain(s) x", ns, "samples",
      "(warmup =", ci$warmup, ") --",
      round(x$extras$run_time, 1), "sec\n")

  cat("  Params :", mi$n_params, "monitored;",
      mi$n_fixed, "fixed,", mi$n_random, "random terms\n")

  # Quick convergence
  if (!is.null(x$extras$convergence$gelman)) {
    rhat <- x$extras$convergence$gelman$psrf[, "Point est."]
    max_rhat <- max(rhat, na.rm = TRUE)
    flag <- if (max_rhat < 1.05) " [OK]" else
            if (max_rhat < 1.10) " [borderline]" else " [!]"
    cat("  Max Rhat:", round(max_rhat, 3), flag, "\n")
  }
  min_eff <- if (!is.null(x$extras$convergence$n_eff))
               min(x$extras$convergence$n_eff, na.rm = TRUE) else NA
  if (!is.na(min_eff))
    cat("  Min ESS:", round(min_eff, 0), "\n")

  cat(strrep("-", 55), "\n")
  cat("  $glm    -- GLM-compatible (summary, emmeans, etc.)\n")
  cat("  $greta  -- native greta (draws, model, calculate)\n")
  cat("  $extras -- diagnostics, BLUPs, variance components\n")

  invisible(x)
}

#' Summarise a bayesreml object
#'
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the posterior summary
#' @export
summary.bayesreml <- function(object, ...) {
  ci <- object$extras$call_info
  mi <- object$extras$model_info

  cat("Bayesian mixed model summary  [bayesreml]\n")
  cat(strrep("=", 60), "\n")
  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) cat("  Random :", deparse(ci$random), "\n")
  cat("  Family :", mi$family, "/", mi$link, "\n")
  cat("  N =", mi$n_obs, ", chains =", ci$chains,
      ", samples =", ci$n_samples, "\n\n")

  # Fixed effects
  cat("-- Fixed effects (posterior) ", strrep("-", 33), "\n")
  beta <- coef(object)
  if (length(beta) > 0) {
    ci_mat <- confint(object)
    fx_df <- data.frame(
      Estimate   = beta,
      Post.SD    = sqrt(diag(vcov(object))),
      `2.5%`     = ci_mat[, 1],
      `97.5%`    = ci_mat[, 2],
      check.names = FALSE
    )
    print(round(fx_df, 4))
  } else {
    cat("  (none)\n")
  }

  # Variance components
  cat("\n-- Variance components ", strrep("-", 38), "\n")
  vc <- object$extras$variance_comps
  if (!is.null(vc) && nrow(vc) > 0) {
    print(data.frame(
      Component = vc$component,
      Estimate  = round(vc$estimate, 4),
      SD        = round(vc$sd, 4),
      `2.5%`    = round(vc$q2.5, 4),
      `97.5%`   = round(vc$q97.5, 4),
      check.names = FALSE,
      row.names = NULL
    ))
  } else {
    cat("  (none available)\n")
  }

  # Convergence
  cat("\n-- Convergence ", strrep("-", 45), "\n")
  conv <- object$extras$convergence
  if (!is.null(conv$gelman)) {
    rhat <- conv$gelman$psrf[, "Point est."]
    cat("  Rhat range:", round(min(rhat, na.rm = TRUE), 3), "-",
        round(max(rhat, na.rm = TRUE), 3), "\n")
  }
  if (!is.null(conv$n_eff)) {
    cat("  ESS  range:", round(min(conv$n_eff, na.rm = TRUE), 0), "-",
        round(max(conv$n_eff, na.rm = TRUE), 0), "\n")
  }
  cat("  Run time  :", round(object$extras$run_time, 1), "sec\n")

  invisible(object$extras$summary)
}

#' Extract fixed effect coefficients
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Named numeric vector of posterior mean fixed effects
#' @export
coef.bayesreml <- function(object, ...) {
  object$glm$coefficients
}

#' Extract variance-covariance matrix of fixed effects
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Posterior covariance matrix of fixed effect coefficients
#' @export
vcov.bayesreml <- function(object, ...) {
  attr(object$glm, "posterior_vcov")
}

#' Credible intervals for fixed effects
#'
#' Returns posterior quantile-based credible intervals, not
#' frequentist confidence intervals.
#'
#' @param object A bayesreml object
#' @param parm Parameter names (NULL for all fixed effects)
#' @param level Credible level (default 0.95)
#' @param ... Additional arguments (ignored)
#' @return Matrix with lower and upper credible bounds
#' @export
confint.bayesreml <- function(object, parm = NULL, level = 0.95, ...) {
  draws <- object$greta$draws
  all_draws <- do.call(rbind, lapply(draws, as.matrix))

  # Get fixed effect column names
  beta_names <- names(coef(object))
  if (length(beta_names) == 0) return(matrix(nrow = 0, ncol = 2))

  # Map coefficient names back to draw column names
  coef_info <- .extract_fixed_coefs(object$extras$parse_info$fixed, draws)
  # Use the raw draw columns
  col_names <- colnames(all_draws)

  # Find the draws columns for fixed effects
  fixed_info <- object$extras$parse_info$fixed
  fixed_draw_cols <- .get_fixed_draw_columns(fixed_info, col_names)

  if (length(fixed_draw_cols) == 0)
    return(matrix(nrow = 0, ncol = 2))

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  avail <- fixed_draw_cols[fixed_draw_cols %in% col_names]
  ci_mat <- t(apply(all_draws[, avail, drop = FALSE], 2, quantile,
                    probs = probs))
  rownames(ci_mat) <- beta_names[seq_along(avail)]
  colnames(ci_mat) <- paste0(round(probs * 100, 1), "%")

  if (!is.null(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }

  ci_mat
}

# Helper: get draw column names for fixed effects
.get_fixed_draw_columns <- function(fixed_info, col_names) {
  result <- character(0)
  if (fixed_info$intercept && "mu_atg" %in% col_names) {
    result <- c(result, "mu_atg")
  }
  for (term in fixed_info$terms) {
    if (term$type == "factor") {
      prefix <- if (fixed_info$intercept) paste0("tau_", term$var) else
                paste0("alpha_", term$var)
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      result <- c(result, matches)
    } else if (term$type == "continuous") {
      nm <- paste0("beta_", term$var)
      if (nm %in% col_names) result <- c(result, nm)
    } else if (term$type == "factor_interaction") {
      tag <- paste(term$vars, collapse = "_x_")
      prefix <- if (fixed_info$intercept) paste0("tau_", tag) else
                paste0("alpha_", tag)
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      result <- c(result, matches)
    } else if (term$type %in% c("interaction", "expression")) {
      tag <- if (term$type == "interaction") {
        paste0(term$vars[1], "_x_", term$vars[2])
      } else {
        gsub("[^A-Za-z0-9_]", "_", term$label)
      }
      nm <- paste0("beta_", tag)
      if (nm %in% col_names) result <- c(result, nm)
    }
  }
  result
}

#' Extract fitted values
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of posterior mean fitted values
#' @export
fitted.bayesreml <- function(object, ...) {
  object$glm$fitted.values
}

#' Extract residuals
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of residuals (observed - fitted)
#' @export
residuals.bayesreml <- function(object, ...) {
  object$glm$residuals
}

#' Predict from a bayesreml model
#'
#' @param object A bayesreml object
#' @param newdata Optional new data frame for prediction. If NULL, returns
#'   fitted values from the original data.
#' @param type `"link"` for linear predictor, `"response"` for response scale.
#' @param se.fit Logical: return standard errors?
#' @param ... Additional arguments (ignored)
#' @return If `se.fit = FALSE`, a numeric vector. If `se.fit = TRUE`, a list
#'   with `fit` and `se.fit`.
#' @export
predict.bayesreml <- function(object, newdata = NULL,
                               type = c("response", "link"),
                               se.fit = FALSE, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "link") {
      pred <- object$glm$linear.predictors
    } else {
      pred <- object$glm$fitted.values
    }
    if (se.fit) {
      return(list(fit = pred, se.fit = rep(NA_real_, length(pred))))
    }
    return(pred)
  }

  # For new data, use the model matrix approach with posterior means
  tryCatch({
    mm <- model.matrix(object$glm$formula, data = newdata)
    beta <- coef(object)

    # Match columns
    common <- intersect(colnames(mm), names(beta))
    if (length(common) == 0) {
      warning("No matching coefficients for newdata prediction.")
      pred <- rep(NA_real_, nrow(newdata))
    } else {
      pred <- as.numeric(mm[, common, drop = FALSE] %*% beta[common])
    }

    if (type == "response") {
      fam_link <- object$extras$parse_info$family
      pred <- switch(fam_link$link,
        "identity" = pred,
        "log"      = exp(pred),
        "logit"    = 1 / (1 + exp(-pred)),
        "probit"   = pnorm(pred),
        pred
      )
    }

    if (se.fit) {
      V <- vcov(object)
      if (length(common) > 0 && nrow(V) > 0) {
        se <- sqrt(rowSums(
          (mm[, common, drop = FALSE] %*% V[common, common]) *
           mm[, common, drop = FALSE]
        ))
      } else {
        se <- rep(NA_real_, nrow(newdata))
      }
      return(list(fit = pred, se.fit = se))
    }

    pred
  }, error = function(e) {
    warning("Prediction with newdata failed: ", conditionMessage(e))
    rep(NA_real_, nrow(newdata))
  })
}

#' Log-likelihood (approximate)
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return logLik object
#' @export
logLik.bayesreml <- function(object, ...) {
  y <- object$glm$y
  fitted <- fitted(object)
  fam_link <- object$extras$parse_info$family

  ll <- tryCatch({
    if (fam_link$family == "gaussian") {
      # Get residual SD from posterior
      all_draws <- do.call(rbind, lapply(object$greta$draws, as.matrix))
      sigma_cols <- grep("^sigma_e_atg", colnames(all_draws), value = TRUE)
      if (length(sigma_cols) > 0) {
        sigma_e <- mean(all_draws[, sigma_cols[1]])
      } else {
        sigma_e <- sqrt(mean((y - fitted)^2))
      }
      sum(dnorm(y, mean = fitted, sd = sigma_e, log = TRUE))
    } else if (fam_link$family %in% c("binomial", "binary")) {
      p <- pmax(pmin(fitted, 1 - 1e-10), 1e-10)
      sum(dbinom(y, size = 1, prob = p, log = TRUE))
    } else if (fam_link$family == "poisson") {
      sum(dpois(y, lambda = pmax(fitted, 1e-10), log = TRUE))
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  structure(ll,
    df = object$extras$model_info$n_params,
    nobs = object$extras$model_info$n_obs,
    class = "logLik"
  )
}

#' Number of observations
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return Integer
#' @export
nobs.bayesreml <- function(object, ...) {
  object$extras$model_info$n_obs
}

#' Extract model family
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return A family object
#' @export
family.bayesreml <- function(object, ...) {
  object$glm$family
}

#' Extract model formula
#' @param x A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return The fixed formula
#' @export
formula.bayesreml <- function(x, ...) {
  x$glm$formula
}

#' Extract model matrix
#' @param object A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return The fixed effects model matrix
#' @export
model.matrix.bayesreml <- function(object, ...) {
  model.matrix(object$glm$formula, data = object$glm$data)
}

#' Update a bayesreml model
#'
#' Re-fit the model with modified arguments.
#'
#' @param object A bayesreml object
#' @param ... Arguments to override (e.g., `n_samples = 2000`)
#' @return A new bayesreml object
#' @export
update.bayesreml <- function(object, ...) {
  cl <- object$extras$call_info
  dots <- list(...)

  args <- list(
    fixed          = cl$fixed,
    random         = cl$random,
    rcov           = cl$rcov,
    data           = object$glm$data,
    family         = cl$family,
    link           = cl$link,
    known_matrices = cl$known_matrices,
    weights        = cl$weights,
    n_samples      = cl$n_samples,
    warmup         = cl$warmup,
    chains         = cl$chains,
    prior_fixed_sd = cl$prior_fixed_sd,
    prior_vc_sd    = cl$prior_vc_sd
  )

  # Override with user-supplied arguments
  for (nm in names(dots)) {
    args[[nm]] <- dots[[nm]]
  }

  do.call(bayesreml, args)
}

#' Compare bayesreml models
#'
#' Bayesian model comparison using WAIC or DIC.
#'
#' @param object A bayesreml object
#' @param ... Additional bayesreml objects to compare
#' @return A data frame with model comparison statistics
#' @export
anova.bayesreml <- function(object, ...) {
  models <- c(list(object), list(...))
  n_models <- length(models)

  # Compute log-likelihoods
  lls <- vapply(models, function(m) {
    ll <- logLik(m)
    as.numeric(ll)
  }, numeric(1))

  n_params <- vapply(models, function(m) {
    m$extras$model_info$n_params
  }, integer(1))

  # Approximate DIC-like comparison
  dic <- -2 * lls + 2 * n_params

  model_names <- paste0("Model ", seq_len(n_models))

  result <- data.frame(
    Model   = model_names,
    logLik  = round(lls, 2),
    npar    = n_params,
    DIC     = round(dic, 2),
    delta   = round(dic - min(dic), 2),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  cat("Bayesian model comparison\n")
  cat(strrep("-", 50), "\n")
  print(result)
  cat("\nNote: DIC is approximate. For rigorous comparison,\n",
      "use loo::loo() on the greta draws.\n")

  invisible(result)
}
