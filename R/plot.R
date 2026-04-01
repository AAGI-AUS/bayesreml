# Plot methods for bayesreml

#' Plot diagnostics for a bayesreml model
#'
#' @param x A bayesreml object
#' @param type Character: type of plot to produce.
#'   - `"diagnostics"`: trace plots + density (requires bayesplot)
#'   - `"residuals"`: residuals vs fitted + QQ plot
#'   - `"effects"`: forest plot of fixed effects with credible intervals
#'   - `"variance"`: bar chart of variance components with credible intervals
#'   - `"blups"`: caterpillar plot of BLUPs
#'   - `"pp_check"`: posterior predictive check (observed vs replicated)
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.bayesreml <- function(x, type = c("diagnostics", "residuals", "effects",
                                        "variance", "blups", "pp_check"),
                            ...) {
  type <- match.arg(type)

  switch(type,
    "diagnostics" = .plot_diagnostics(x, ...),
    "residuals"   = .plot_residuals(x, ...),
    "effects"     = .plot_effects(x, ...),
    "variance"    = .plot_variance(x, ...),
    "blups"       = .plot_blups(x, ...),
    "pp_check"    = .plot_pp_check(x, ...)
  )
}

# MCMC diagnostics: trace + density
.plot_diagnostics <- function(x, ...) {
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    draws <- x$greta$draws
    p1 <- bayesplot::mcmc_trace(draws, ...)
    print(p1)
  } else {
    # Base R fallback
    draws <- x$greta$draws
    all_draws <- do.call(rbind, lapply(draws, as.matrix))
    n_params <- min(ncol(all_draws), 6)
    par_names <- colnames(all_draws)[seq_len(n_params)]

    old_par <- par(mfrow = c(n_params, 2), mar = c(3, 3, 2, 1))
    on.exit(par(old_par))

    for (i in seq_len(n_params)) {
      nm <- par_names[i]
      # Trace
      for (ch in seq_along(draws)) {
        vals <- as.matrix(draws[[ch]])[, nm]
        if (ch == 1) {
          plot(vals, type = "l", main = paste("Trace:", nm),
               xlab = "", ylab = "", col = ch)
        } else {
          lines(vals, col = ch)
        }
      }
      # Density
      vals <- all_draws[, nm]
      plot(density(vals), main = paste("Density:", nm),
           xlab = "", ylab = "")
    }
  }
  invisible(NULL)
}

# Residual diagnostics
.plot_residuals <- function(x, ...) {
  fitted_vals <- fitted(x)
  resid_vals  <- residuals(x)

  if (any(is.na(fitted_vals)) || any(is.na(resid_vals))) {
    message("Cannot plot residuals: fitted values contain NA.")
    return(invisible(NULL))
  }

  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  # 1. Residuals vs Fitted
  plot(fitted_vals, resid_vals,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted", pch = 16, col = "#00000060")
  abline(h = 0, lty = 2, col = "red")
  lines(lowess(fitted_vals, resid_vals), col = "blue")

  # 2. QQ plot
  qqnorm(resid_vals, main = "Normal Q-Q", pch = 16, col = "#00000060")
  qqline(resid_vals, col = "red")

  # 3. Scale-Location
  plot(fitted_vals, sqrt(abs(resid_vals)),
       xlab = "Fitted values", ylab = "sqrt(|Residuals|)",
       main = "Scale-Location", pch = 16, col = "#00000060")
  lines(lowess(fitted_vals, sqrt(abs(resid_vals))), col = "blue")

  # 4. Histogram of residuals
  hist(resid_vals, breaks = 30, main = "Residual Distribution",
       xlab = "Residuals", col = "lightblue", border = "white")

  invisible(NULL)
}

# Forest plot of fixed effects
.plot_effects <- function(x, level = 0.95, ...) {
  beta <- coef(x)
  if (length(beta) == 0) {
    message("No fixed effects to plot.")
    return(invisible(NULL))
  }

  ci <- confint(x, level = level)

  n <- length(beta)
  y_pos <- seq_len(n)

  old_par <- par(mar = c(4, max(nchar(names(beta))) * 0.6 + 2, 2, 1))
  on.exit(par(old_par))

  xlim <- range(c(ci, 0)) * 1.1

  plot(beta, y_pos, xlim = xlim, yaxt = "n",
       xlab = "Estimate", ylab = "",
       main = paste0("Fixed effects (", round(level * 100), "% CrI)"),
       pch = 16, cex = 1.2)
  segments(ci[, 1], y_pos, ci[, 2], y_pos, lwd = 2)
  abline(v = 0, lty = 2, col = "grey50")
  axis(2, at = y_pos, labels = names(beta), las = 1, cex.axis = 0.8)

  invisible(NULL)
}

# Variance components plot
.plot_variance <- function(x, ...) {
  vc <- x$extras$variance_comps
  if (is.null(vc) || nrow(vc) == 0) {
    message("No variance components to plot.")
    return(invisible(NULL))
  }

  n <- nrow(vc)
  y_pos <- seq_len(n)

  old_par <- par(mar = c(4, max(nchar(vc$component)) * 0.5 + 2, 2, 1))
  on.exit(par(old_par))

  xlim <- c(0, max(vc$q97.5, na.rm = TRUE) * 1.1)

  plot(vc$estimate, y_pos, xlim = xlim, yaxt = "n",
       xlab = "Estimate", ylab = "",
       main = "Variance components (95% CrI)",
       pch = 16, cex = 1.2)
  segments(vc$q2.5, y_pos, vc$q97.5, y_pos, lwd = 2)
  axis(2, at = y_pos, labels = vc$component, las = 1, cex.axis = 0.7)

  invisible(NULL)
}

# Caterpillar plot of BLUPs
.plot_blups <- function(x, ...) {
  blups <- x$extras$blups
  if (length(blups) == 0) {
    message("No BLUPs available.")
    return(invisible(NULL))
  }

  n_terms <- length(blups)
  old_par <- par(mfrow = c(1, n_terms), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  for (nm in names(blups)) {
    vals <- blups[[nm]]
    n <- length(vals)
    ord <- order(vals)

    plot(vals[ord], seq_len(n), pch = 16, cex = 0.6,
         xlab = "BLUP", ylab = "Rank",
         main = nm, col = "#00000080")
    abline(v = 0, lty = 2, col = "red")
  }

  invisible(NULL)
}

# Posterior predictive check
.plot_pp_check <- function(x, ...) {
  y <- x$glm$y
  fitted_vals <- fitted(x)

  if (any(is.na(fitted_vals))) {
    message("Cannot produce posterior predictive check: fitted values contain NA.")
    return(invisible(NULL))
  }

  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  # 1. Observed vs predicted
  lims <- range(c(y, fitted_vals))
  plot(fitted_vals, y, pch = 16, col = "#00000060",
       xlab = "Predicted", ylab = "Observed",
       main = "Observed vs Predicted",
       xlim = lims, ylim = lims)
  abline(0, 1, col = "red", lty = 2)

  # 2. Density overlay
  plot(density(y), main = "Posterior Predictive",
       xlab = "Value", col = "black", lwd = 2)
  lines(density(fitted_vals), col = "blue", lwd = 2, lty = 2)
  legend("topright", c("Observed", "Predicted"),
         col = c("black", "blue"), lty = c(1, 2), lwd = 2, bty = "n")

  invisible(NULL)
}
