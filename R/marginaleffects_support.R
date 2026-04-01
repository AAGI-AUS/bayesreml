# marginaleffects integration for bayesreml
#
# marginaleffects works via coef(), vcov(), predict() S3 methods,
# which are already defined in methods.R. This file provides any
# additional support needed.

# No additional methods needed - marginaleffects dispatches via:
#   coef.bayesreml()    -> methods.R
#   vcov.bayesreml()    -> methods.R
#   predict.bayesreml() -> methods.R (handles newdata + type)
#   formula.bayesreml() -> methods.R
#   family.bayesreml()  -> methods.R
#
# For full Bayesian posterior propagation, users should work with
# fit$greta$draws directly:
#
#   library(marginaleffects)
#   avg_slopes(fit)                         # uses point estimates + vcov
#   predictions(fit, newdata = nd)          # uses predict.bayesreml

# Register bayesreml as supported by marginaleffects
# This is called in .onLoad
.register_marginaleffects <- function() {
  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    # marginaleffects >= 0.18 auto-discovers via S3 methods
    # No explicit registration needed
    TRUE
  } else {
    FALSE
  }
}
