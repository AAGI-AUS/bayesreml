# Code generation for bayesreml
# Generates greta code strings from parsed formula descriptors.
# Not exported.

# Generate code for fixed effects
.code_fixed <- function(ctx, fixed_info) {
  psd <- ctx$prior_fx
  ctx <- .add(ctx, "# -- Fixed effects -------------------------------------------")

  # Intercept
  if (fixed_info$intercept) {
    ctx <- .add(ctx, "mu_atg <- normal(0, ", psd, ")")
    ctx <- .add_param(ctx, "mu_atg")
    ctx <- .add_pred(ctx, "mu_atg")
  }

  for (term in fixed_info$terms) {
    ctx <- switch(term$type,
      "factor" = {
        tag <- term$var
        if (fixed_info$intercept) {
          ctx <- .add(ctx, "tau_", tag, " <- normal(0, ", psd,
                      ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("tau_", tag))
          ctx <- .add_pred(ctx, paste0("tau_", tag, "[", tag, "_id]"))
        } else {
          ctx <- .add(ctx, "alpha_", tag, " <- normal(0, ", psd,
                      ", dim = n_", tag, ")")
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
          ctx <- .add(ctx, "tau_", tag, " <- normal(0, ", psd,
                      ", dim = n_", tag, ")")
          ctx <- .add_param(ctx, paste0("tau_", tag))
          ctx <- .add_pred(ctx, paste0("tau_", tag, "[", tag, "_id]"))
        } else {
          ctx <- .add(ctx, "alpha_", tag, " <- normal(0, ", psd,
                      ", dim = n_", tag, ")")
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
        inner    <- gsub("^I\\((.*)\\)$", "\\1", trimws(raw_lbl))
        ctx <- .add(ctx, "beta_", tag, " <- normal(0, 2)")
        ctx <- .add(ctx, "x_", tag,
                    " <- as_data(as.numeric(with(dat_atg, ", inner, ")))")
        ctx <- .add_param(ctx, paste0("beta_", tag))
        ctx <- .add_pred(ctx, paste0("beta_", tag, " * x_", tag))
        ctx
      },
      ctx
    )
  }
  ctx
}

# Generate code for random effects
.code_random <- function(ctx, random_terms, data, known_matrices) {
  if (length(random_terms) == 0) return(ctx)
  ctx <- .add(ctx, "# -- Random effects ------------------------------------------")

  for (i in seq_along(random_terms)) {
    term <- random_terms[[i]]
    pvc  <- ctx$prior_vc

    ctx <- switch(term$type,

      # iid (simple, ide, id)
      "simple" = , "ide" = , "id" = {
        tag <- term$var
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- ", tag, "_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # GBLUP: vm(geno, G)
      "vm" = {
        tag <- term$var; mat <- term$mat
        ctx <- .add(ctx, "L_G_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- as_data(L_G_", tag,
                    ") %*% (z_", tag, " * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # Animal model: ped(animal, A)
      "ped" = {
        tag <- term$var; mat <- term$mat
        ctx <- .add(ctx, "L_A_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- as_data(L_A_", tag,
                    ") %*% (z_", tag, " * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[", tag, "_id]"))
        ctx
      },

      # DIAG(outer) x I(inner): at(env):geno
      "at_simple" = {
        outer <- term$outer; inner <- term$inner
        tag   <- paste0(inner, "_", outer)
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc,
                    ", dim = n_", outer, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = c(n_", inner,
                    ", n_", outer, "))")
        ctx <- .add(ctx, tag, "_sc  <- sweep(", tag, "_raw, 2, sigma_",
                    tag, ", \"*\")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0(tag, "_sc[cbind(", inner, "_id, ",
                                     outer, "_id)]"))
        ctx
      },

      # US(outer) x I(inner): us(env):id(geno)
      "us_gxe" = {
        outer <- term$outer; inner <- term$inner
        tag   <- paste0(inner, "_", outer, "_us")
        ctx <- .add(ctx, "L_corr_", tag,
                    " <- lkj_correlation(eta = 2, dimension = n_", outer, ")")
        ctx <- .add(ctx, "sd_", tag, " <- greta::lognormal(0, ", pvc,
                    ", dim = n_", outer, ")")
        ctx <- .add(ctx, "L_", tag, " <- sweep(L_corr_", tag, ", 2, sd_",
                    tag, ", \"*\")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = c(n_", inner,
                    ", n_", outer, "))")
        ctx <- .add(ctx, tag, "_mat <- ", tag, "_raw %*% t(L_", tag, ")")
        ctx <- .add_param(ctx, paste0("sd_", tag))
        ctx <- .add_pred(ctx, paste0(tag, "_mat[cbind(", inner, "_id, ",
                                     outer, "_id)]"))
        ctx
      },

      # FA(k)(outer) x I(inner): fa(env,k):id(geno)
      "fa_gxe" = {
        outer <- term$outer; inner <- term$inner; k <- term$k
        tag   <- paste0(inner, "_", outer, "_fa", k)
        ctx <- .add(ctx, "Lambda_", tag, " <- normal(0, 1, dim = c(n_",
                    outer, ", ", k, "L))")
        ctx <- .add(ctx, "psi_", tag, "    <- greta::lognormal(0, ", pvc,
                    ", dim = n_", outer, ")")
        ctx <- .add(ctx, "F_", tag, "      <- normal(0, 1, dim = c(n_",
                    inner, ", ", k, "L))")
        ctx <- .add(ctx, "delta_", tag,
                    "  <- normal(0, 1, dim = c(n_", inner, ", n_", outer, "))")
        ctx <- .add(ctx, "delta_sc_", tag, " <- sweep(delta_", tag,
                    ", 2, psi_", tag, ", \"*\")")
        ctx <- .add(ctx, "g_mat_", tag, "  <- F_", tag,
                    " %*% t(Lambda_", tag, ") + delta_sc_", tag)
        ctx <- .add_param(ctx, paste0("Lambda_", tag), paste0("psi_", tag))
        ctx <- .add_pred(ctx, paste0("g_mat_", tag, "[cbind(", inner,
                                     "_id, ", outer, "_id)]"))
        ctx
      },

      # US(outer) x Struct(inner): vm(geno,G):id(env)
      "vm_gxe" = {
        inner <- term$inner; outer <- term$outer; mat <- term$mat
        tag   <- paste0(inner, "_", outer, "_vm_gxe")
        ctx <- .add(ctx, "L_G_", tag, " <- t(chol(", mat, "))")
        ctx <- .add(ctx, "L_corr_", tag,
                    " <- lkj_correlation(eta = 2, dimension = n_", outer, ")")
        ctx <- .add(ctx, "sd_", tag, " <- greta::lognormal(0, ", pvc,
                    ", dim = n_", outer, ")")
        ctx <- .add(ctx, "L_", tag, " <- sweep(L_corr_", tag, ", 2, sd_",
                    tag, ", \"*\")")
        ctx <- .add(ctx, "z_", tag, " <- normal(0, 1, dim = c(n_", inner,
                    ", n_", outer, "))")
        ctx <- .add(ctx, "uA_", tag, " <- as_data(L_G_", tag,
                    ") %*% z_", tag)
        ctx <- .add(ctx, "g_mat_", tag, " <- uA_", tag,
                    " %*% t(L_", tag, ")")
        ctx <- .add_param(ctx, paste0("sd_", tag))
        ctx <- .add_pred(ctx, paste0("g_mat_", tag, "[cbind(", inner,
                                     "_id, ", outer, "_id)]"))
        ctx
      },

      # Nested random effect: outer:inner
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

      # Multi-way combination: A:B:C
      "combo" = {
        tag <- paste(rev(term$vars), collapse = "_in_")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = n_", tag, ")")
        ctx <- .add(ctx, "u_", tag, " <- ", tag, "_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "[combo_id_", tag, "]"))
        ctx
      },

      # Spatial AR1 approximation
      "ar1_spatial" = {
        rv <- term$row_var; cv <- term$col_var
        tag <- paste0("sp_", rv, "_", cv)
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_row_raw <- normal(0, 1, dim = n_", rv, ")")
        ctx <- .add(ctx, tag, "_col_raw <- normal(0, 1, dim = n_", cv, ")")
        ctx <- .add(ctx, "u_", tag, "_row <- ", tag,
                    "_row_raw * sigma_", tag)
        ctx <- .add(ctx, "u_", tag, "_col <- ", tag,
                    "_col_raw * sigma_", tag)
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("u_", tag, "_row[", rv, "_id] + u_",
                                     tag, "_col[", cv, "_id]"))
        ctx
      },

      # P-spline: spl(x)
      "spline" = {
        vname <- term$var
        tag   <- paste0("spl_", vname)
        ctx <- .add(ctx, "x_std_", tag,
                    " <- as.numeric(scale(dat_atg$", vname, "))")
        ctx <- .add(ctx, "B_", tag, " <- splines::bs(x_std_", tag,
                    ", df=8, degree=3, intercept=FALSE)")
        ctx <- .add(ctx, "B_g_", tag, " <- as_data(B_", tag, ")")
        ctx <- .add(ctx, "sigma_", tag, " <- greta::lognormal(0, ", pvc, ")")
        ctx <- .add(ctx, tag, "_raw <- normal(0, 1, dim = ncol(B_", tag, "))")
        ctx <- .add(ctx, "f_", tag, " <- B_g_", tag, " %*% (", tag,
                    "_raw * sigma_", tag, ")")
        ctx <- .add_param(ctx, paste0("sigma_", tag))
        ctx <- .add_pred(ctx, paste0("f_", tag))
        ctx
      },

      ctx  # default: skip unknown terms
    )
  }
  ctx
}

# Generate code for residual (rcov) structure
.code_rcov <- function(ctx, rcov_terms, data) {
  ctx <- .add(ctx, "# -- Residual ------------------------------------------------")
  pvc <- ctx$prior_vc

  has_at_units <- any(vapply(rcov_terms,
    function(t) t$type == "at_units", logical(1)))

  if (has_at_units) {
    term <- Filter(function(t) t$type == "at_units", rcov_terms)[[1]]
    tag  <- term$var
    ctx  <- .add(ctx, "sigma_e_atg <- greta::lognormal(0, ", pvc,
                 ", dim = n_", tag, ")")
    ctx  <- .add_param(ctx, "sigma_e_atg")
    ctx$resid_expr <- paste0("sigma_e_atg[", tag, "_id]")
  } else {
    ctx  <- .add(ctx, "sigma_e_atg <- greta::lognormal(0, ", pvc, ")")
    ctx  <- .add_param(ctx, "sigma_e_atg")
    ctx$resid_expr <- "sigma_e_atg"
  }

  ctx
}

# Generate code for the linear predictor assembly
.code_predictor <- function(ctx, fixed_info) {
  ctx <- .add(ctx, "# -- Linear predictor ----------------------------------------")
  if (length(ctx$predictor) == 0) {
    ctx <- .add(ctx, "mu_i_atg <- zeros(N_atg, 1)")
  } else {
    parts <- ctx$predictor
    ctx   <- .add(ctx, "mu_i_atg <- ", paste(parts, collapse = " + "))
  }
  ctx
}

# Generate code for the likelihood
.code_likelihood <- function(ctx, fixed_info, rcov_terms, data, weights) {
  ctx  <- .add(ctx, "# -- Likelihood ----------------------------------------------")
  resp <- fixed_info$response
  fl   <- ctx$fam_link
  resid_sd <- if (!is.null(ctx$resid_expr)) ctx$resid_expr else "sigma_e_atg"

  ctx <- .add(ctx, "y_atg_obs <- as_data(y_atg)")

  if (!is.null(weights)) {
    ctx <- .add(ctx, "wt_atg_data <- as_data(wt_atg)")
    resid_sd <- paste0(resid_sd, " / sqrt(wt_atg_data)")
  }

  if (fl$family == "gaussian") {
    ctx <- .add(ctx, "distribution(y_atg_obs) <- normal(mu_i_atg, ",
                resid_sd, ")")

  } else if (fl$family %in% c("binomial", "binary")) {
    if (fl$link == "logit") {
      ctx <- .add(ctx, "p_atg <- ilogit(mu_i_atg)")
    } else {
      ctx <- .add(ctx, "p_atg <- iprobit(mu_i_atg)")
    }
    y_vals <- data[[resp]]
    if (all(y_vals %in% c(0, 1))) {
      ctx <- .add(ctx, "distribution(y_atg_obs) <- bernoulli(p_atg)")
    } else {
      ctx <- .add(ctx, "distribution(y_atg_obs) <- bernoulli(p_atg)")
    }

  } else if (fl$family == "poisson") {
    ctx <- .add(ctx, "lam_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "distribution(y_atg_obs) <- poisson(lam_atg)")

  } else if (fl$family %in% c("negative_binomial", "negbinom")) {
    ctx <- .add(ctx, "r_nb_atg <- greta::lognormal(0, 1)")
    ctx <- .add_param(ctx, "r_nb_atg")
    ctx <- .add(ctx, "mu_nb_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "prob_nb_atg <- r_nb_atg / (r_nb_atg + mu_nb_atg)")
    ctx <- .add(ctx,
      "distribution(y_atg_obs) <- negative_binomial(r_nb_atg, prob_nb_atg)")

  } else if (fl$family == "gamma") {
    ctx <- .add(ctx, "shape_atg <- greta::lognormal(0, 1)")
    ctx <- .add_param(ctx, "shape_atg")
    ctx <- .add(ctx, "mu_gam_atg <- exp(mu_i_atg)")
    ctx <- .add(ctx, "rate_atg <- shape_atg / mu_gam_atg")
    ctx <- .add(ctx,
      "distribution(y_atg_obs) <- greta::gamma(shape_atg, rate_atg)")

  } else if (fl$family == "beta") {
    ctx <- .add(ctx, "phi_b_atg <- greta::lognormal(1, 1)")
    ctx <- .add_param(ctx, "phi_b_atg")
    ctx <- .add(ctx, "mu_b_atg <- ilogit(mu_i_atg)")
    ctx <- .add(ctx, "alpha_b_atg <- mu_b_atg * phi_b_atg")
    ctx <- .add(ctx, "beta_b_atg  <- (1 - mu_b_atg) * phi_b_atg")
    ctx <- .add(ctx,
      "distribution(y_atg_obs) <- greta::beta(alpha_b_atg, beta_b_atg)")
  }

  ctx
}

# Generate code for model() + mcmc()
.code_model <- function(ctx, n_samples, warmup, chains, mcmc_verbose = TRUE) {
  ctx <- .add(ctx,
    "# -- Model and MCMC -----------------------------------------")
  params_str <- paste(unique(ctx$params), collapse = ", ")
  ctx <- .add(ctx, "atg_model <- greta::model(", params_str, ")")
  ctx <- .add(ctx, "atg_draws <- greta::mcmc(atg_model, n_samples = ",
              n_samples, ", warmup = ", warmup, ", chains = ", chains,
              ", verbose = ", mcmc_verbose, ")")
  ctx
}
