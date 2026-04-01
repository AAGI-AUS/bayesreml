# Package load hooks

.onLoad <- function(libname, pkgname) {
  # Register emmeans methods if emmeans is available
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("bayesreml", pkgname)
  }

  # Register broom S3 methods if generics is available
  if (requireNamespace("generics", quietly = TRUE)) {
    tryCatch({
      registerS3method("tidy", "bayesreml", tidy.bayesreml,
                        envir = asNamespace("generics"))
      registerS3method("glance", "bayesreml", glance.bayesreml,
                        envir = asNamespace("generics"))
      registerS3method("augment", "bayesreml", augment.bayesreml,
                        envir = asNamespace("generics"))
    }, error = function(e) NULL)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "bayesreml ", utils::packageVersion(pkgname),
    " -- Bayesian mixed models with ASReml syntax"
  )

  # Check greta availability
  if (!requireNamespace("greta", quietly = TRUE)) {
    packageStartupMessage(
      "  Note: greta is not installed. Install with:\n",
      "    install.packages('greta')\n",
      "    greta::install_greta_deps()"
    )
  }
}
