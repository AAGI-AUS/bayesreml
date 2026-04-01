# Environment setup for bayesreml
# Populates the evaluation environment with all R objects the generated code
# will reference (ID vectors, level counts, known matrices, etc.)
# Not exported.

# Set up the evaluation environment with data objects
#
# @param ev Environment to populate
# @param fixed_info Parsed fixed formula info
# @param random_terms List of random term descriptors
# @param rcov_terms List of rcov term descriptors
# @param data data.frame
# @param known_matrices Named list of matrices
# @param weights Optional numeric weight vector
.setup_env <- function(ev, fixed_info, random_terms, rcov_terms,
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
      assign(paste0("n_", term$var),   nlevels(f),     envir = ev)
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
      combo <- do.call(paste, c(lapply(term$vars, function(v) data[[v]]),
                                list(sep = "_")))
      f     <- factor(combo)
      assign(paste0(tag, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", tag),  nlevels(f),    envir = ev)
    }
  }

  # Random and rcov id vectors / level counts / matrices
  all_terms <- c(random_terms, rcov_terms)
  for (term in all_terms) {
    .setup_term_env(ev, term, data, known_matrices)
  }

  # Weights
  if (!is.null(weights)) ev$wt_atg <- weights
}

# Set up environment objects for a single random/rcov term
.setup_term_env <- function(ev, term, data, known_matrices) {
  mk_id <- function(vname) {
    if (vname %in% names(data)) {
      f <- factor(data[[vname]])
      assign(paste0(vname, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", vname),  nlevels(f),    envir = ev)
    } else if (paste0(vname, "_id") %in% names(data)) {
      id_col <- data[[paste0(vname, "_id")]]
      assign(paste0(vname, "_id"), id_col,      envir = ev)
      assign(paste0("n_", vname),  max(id_col), envir = ev)
    } else {
      warning("Variable '", vname, "' not found in data.")
    }
  }
  mk_mat <- function(mat_name) {
    if (is.null(mat_name) || is.na(mat_name)) return()
    if (!mat_name %in% names(known_matrices))
      stop("Matrix '", mat_name, "' not found in known_matrices. ",
           "Provide it via known_matrices = list(", mat_name,
           " = <your matrix>).")
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
      nid <- as.integer(factor(paste(data[[term$outer]],
                                     data[[term$inner]], sep = "_")))
      ev[[paste0("nested_id_", tag)]] <- nid
      ev[[paste0("n_",         tag)]] <- length(unique(nid))
    },
    "combo" = {
      tag <- paste(rev(term$vars), collapse = "_in_")
      nid <- as.integer(factor(do.call(paste,
        c(lapply(term$vars, function(v) data[[v]]), list(sep = "_")))))
      ev[[paste0("combo_id_", tag)]] <- nid
      ev[[paste0("n_",        tag)]] <- length(unique(nid))
    },
    "ar1_spatial" = {
      if (!is.null(term$row_var) && term$row_var %in% names(data)) {
        row_vals <- data[[term$row_var]]
        assign(paste0(term$row_var, "_id"),
               as.integer(factor(row_vals)), envir = ev)
        assign(paste0("n_", term$row_var),
               length(unique(row_vals)), envir = ev)
      }
      if (!is.null(term$col_var) && term$col_var %in% names(data)) {
        col_vals <- data[[term$col_var]]
        assign(paste0(term$col_var, "_id"),
               as.integer(factor(col_vals)), envir = ev)
        assign(paste0("n_", term$col_var),
               length(unique(col_vals)), envir = ev)
      }
    }
  )
}
