# broom-style tidy, glance, augment methods for bayesreml

#' Tidy a bayesreml model
#'
#' Returns a tidy data frame of fixed effect estimates with posterior
#' statistics.
#'
#' @param x A bayesreml object
#' @param conf.int Logical: include credible intervals?
#' @param conf.level Credible level (default 0.95)
#' @param effects Which effects to return: `"fixed"` or `"random"`
#' @param ... Additional arguments (ignored)
#' @return A data frame with columns: term, estimate, std.error,
#'   conf.low, conf.high
#' @export
tidy.bayesreml <- function(x, conf.int = TRUE, conf.level = 0.95,
                            effects = c("fixed", "random"), ...) {
  effects <- match.arg(effects)

  if (effects == "fixed") {
    beta <- coef(x)
    if (length(beta) == 0) return(data.frame())

    V <- vcov(x)
    se <- sqrt(diag(V))

    result <- data.frame(
      term      = names(beta),
      estimate  = unname(beta),
      std.error = unname(se),
      stringsAsFactors = FALSE
    )

    if (conf.int) {
      ci <- confint(x, level = conf.level)
      result$conf.low  <- ci[, 1]
      result$conf.high <- ci[, 2]
    }

  } else {
    # Random effects (variance components)
    vc <- x$extras$variance_comps
    if (is.null(vc)) return(data.frame())

    result <- data.frame(
      term      = vc$component,
      estimate  = vc$estimate,
      std.error = vc$sd,
      stringsAsFactors = FALSE
    )

    if (conf.int) {
      result$conf.low  <- vc$q2.5
      result$conf.high <- vc$q97.5
    }
  }

  result
}

#' Glance at a bayesreml model
#'
#' Returns a one-row summary of model-level statistics.
#'
#' @param x A bayesreml object
#' @param ... Additional arguments (ignored)
#' @return A one-row data frame
#' @export
glance.bayesreml <- function(x, ...) {
  mi <- x$extras$model_info
  ll <- logLik(x)

  conv <- x$extras$convergence
  max_rhat <- if (!is.null(conv$gelman)) {
    max(conv$gelman$psrf[, "Point est."], na.rm = TRUE)
  } else NA_real_

  min_ess <- if (!is.null(conv$n_eff)) {
    min(conv$n_eff, na.rm = TRUE)
  } else NA_real_

  data.frame(
    nobs      = mi$n_obs,
    npar      = mi$n_params,
    logLik    = as.numeric(ll),
    family    = mi$family,
    link      = mi$link,
    chains    = x$extras$call_info$chains,
    samples   = x$extras$call_info$n_samples,
    max_rhat  = round(max_rhat, 3),
    min_ess   = round(min_ess, 0),
    run_time  = round(x$extras$run_time, 1),
    stringsAsFactors = FALSE
  )
}

#' Augment a bayesreml model
#'
#' Returns the original data with fitted values, residuals, and other
#' observation-level statistics.
#'
#' @param x A bayesreml object
#' @param data Optional data frame (defaults to original data)
#' @param ... Additional arguments (ignored)
#' @return A data frame with added columns: .fitted, .resid
#' @export
augment.bayesreml <- function(x, data = NULL, ...) {
  if (is.null(data)) data <- x$glm$data

  result <- data
  result$.fitted <- fitted(x)
  result$.resid  <- residuals(x)

  result
}
