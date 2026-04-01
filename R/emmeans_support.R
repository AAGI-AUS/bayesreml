# emmeans integration for bayesreml
#
# Provides recover_data() and emm_basis() methods so that
# emmeans::emmeans(fit, ~ factor) works out of the box.

#' Recover data for emmeans
#'
#' @param object A bayesreml object
#' @param data Optional data frame override (default uses stored data)
#' @param ... Additional arguments passed to emmeans
#' @return A data frame suitable for emmeans
#' @export
recover_data.bayesreml <- function(object, data = NULL, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE))
    stop("Package 'emmeans' is required for this method.")

  glm_obj <- object$glm
  trms    <- glm_obj$terms
  if (is.null(data)) data <- glm_obj$data

  # Build a call that recover_data.call can use
  cl <- glm_obj$call
  emmeans::recover_data(cl, trms, na.action = NULL, data = data, ...)
}

#' Compute emmeans basis
#'
#' @param object A bayesreml object
#' @param trms Terms object
#' @param xlev Factor levels
#' @param grid Reference grid
#' @param ... Additional arguments
#' @return A list with X, bhat, V, dffun, dfargs
#' @export
emm_basis.bayesreml <- function(object, trms, xlev, grid, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE))
    stop("Package 'emmeans' is required for this method.")

  # Delegate to the GLM method
  emm_basis_result <- tryCatch(
    emmeans::emm_basis(object$glm, trms, xlev, grid, ...),
    error = function(e) {
      # Fallback: build manually
      mm <- model.matrix(trms, data = grid, xlev = xlev)
      bhat <- coef(object)
      V <- vcov(object)

      # Match columns
      common <- intersect(colnames(mm), names(bhat))
      if (length(common) == 0) {
        common <- seq_len(min(ncol(mm), length(bhat)))
      }

      list(
        X      = mm[, common, drop = FALSE],
        bhat   = bhat[common],
        V      = V[common, common, drop = FALSE],
        dffun  = function(k, dfargs) Inf,
        dfargs = list()
      )
    }
  )

  emm_basis_result
}
