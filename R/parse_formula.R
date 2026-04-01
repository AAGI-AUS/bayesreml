# Formula parsing for bayesreml
# Translates ASReml-style formulas into internal term descriptors.
# Not exported.

# Parse a standard R fixed formula into response + intercept + terms list
#
# @param fixed Two-sided formula: response ~ predictors
# @param data data.frame
# @return List with response, intercept (logical), terms (list of term descriptors)
.parse_fixed <- function(fixed, data) {
  if (length(fixed) < 3)
    stop("fixed must be a two-sided formula: response ~ ...")
  response <- deparse(fixed[[2]])
  if (!response %in% names(data))
    stop("Response variable '", response, "' not found in data.")

  tt <- tryCatch(terms(fixed, data = data), error = function(e) terms(fixed))
  intercept   <- attr(tt, "intercept") == 1
  term_labels <- attr(tt, "term.labels")

  term_info <- lapply(term_labels, function(lbl) .classify_fixed_term(lbl, data))

  list(response = response, intercept = intercept, terms = term_info)
}

# Classify a single fixed-formula term label
.classify_fixed_term <- function(lbl, data) {
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
#
# @param formula One-sided formula: ~ random_terms
# @param data data.frame
# @return List of enriched term descriptors
.parse_formula <- function(formula, data) {
  rhs   <- if (length(formula) == 3) formula[[3]] else formula[[2]]
  terms <- .walk(rhs)
  lapply(terms, .enrich, data = data)
}

# Recursive AST walker — converts formula expression tree to term descriptors
.walk <- function(expr) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm == "units") return(list(list(type = "units")))
    return(list(list(type = "simple", var = nm)))
  }
  if (!is.call(expr)) return(list())

  fn <- as.character(expr[[1]])

  # Binary operators
  if (fn == "+") return(c(.walk(expr[[2]]), .walk(expr[[3]])))

  if (fn == ":") {
    left  <- .walk(expr[[2]])
    right <- .walk(expr[[3]])
    l <- if (length(left))  left[[1]]  else list(type = "unknown")
    r <- if (length(right)) right[[1]] else list(type = "unknown")
    return(.classify_interaction(l, r))
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
      arg <- expr[[2]]
      grp <- tryCatch({
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

# Classify what an a:b interaction means in ASReml context
.classify_interaction <- function(l, r) {
  lt <- l$type; rt <- r$type
  # FA(k) GxE:  fa(env,k) : id(geno)
  if (lt == "fa" && rt %in% c("id", "ide", "simple")) {
    if (!is.null(l$k) && l$k < 1L)
      stop("fa() requires k >= 1; got k = ", l$k, ".")
    return(list(list(type = "fa_gxe", outer = l$var, k = l$k,
                     inner = r$var)))
  }
  # US GxE:     us(env) : id(geno)
  if (lt == "us" && rt %in% c("id", "ide", "simple"))
    return(list(list(type = "us_gxe", outer = l$var, inner = r$var)))
  # at(env):units  - heterogeneous residual
  if (lt == "at" && rt == "units")
    return(list(list(type = "at_units", var = l$var, level = l$level)))
  # at(env):geno   - DIAG(env) x I(geno)
  if (lt == "at" && rt %in% c("simple", "ide", "id"))
    return(list(list(type = "at_simple", outer = l$var, level = l$level,
                     inner = r$var)))
  # AR1 spatial
  if (lt == "ar1" && rt %in% c("ar1", "id", "ide", "simple", "units"))
    return(list(list(type = "ar1_spatial", row_var = l$var, col_var = r$var,
                     col_ar1 = (rt == "ar1"))))
  # vm(geno,G) : id(env)  - structured GxE with known G
  if (lt == "vm" && rt %in% c("id", "ide", "simple"))
    return(list(list(type = "vm_gxe", inner = l$var, mat = l$mat,
                     outer = r$var)))
  # simple:simple  - nested random effect
  if (lt %in% c("simple", "ide", "id") && rt %in% c("simple", "ide", "id"))
    return(list(list(type = "nested", outer = l$var, inner = r$var)))
  # nested:simple  - three-way combination
  if (lt == "nested" && rt %in% c("simple", "ide", "id"))
    return(list(list(type = "combo", vars = c(l$outer, l$inner, r$var))))
  # combo:simple   - four-way+ combination
  if (lt == "combo" && rt %in% c("simple", "ide", "id"))
    return(list(list(type = "combo", vars = c(l$vars, r$var))))
  # fallback
  list(list(type = "interaction_generic", left = l, right = r))
}

# Enrich a term descriptor with data-derived information
.enrich <- function(term, data) {
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
    "fa"     = .add_var_info(term),
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
      if (term$outer %in% names(data))
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data))
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    "at_simple" = {
      if (term$outer %in% names(data))
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data))
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    "at_units" = {
      if (term$var %in% names(data))
        term$n_var <- nlevels(factor(data[[term$var]]))
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
      if (term$outer %in% names(data))
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      if (term$inner %in% names(data))
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      term
    },
    term
  )
}
