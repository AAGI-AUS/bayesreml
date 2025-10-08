# This is the code updated based on comments from Francis (email Friday, 19 September 2025), coding fixed effects in adifferent way.
# 
# One more minor improvement made is replacing `as.numeric` with `as.integer` for factors.
# 
# From eyeballing outputs, `summary(m1_greta)` estimates still differ from the initial "model.matrix" approach.

# source("r_code_04a_greta_test_v4.R")

##############
# Author: Max Moldovan
# max.moldovan@adelaide.edu.au
# version 0.0.1, 12 Sep 2025
# version 0.0.2, 17 Sep 2025, GitHub Francis
# version 0.0.3, 18 Sep 2025, attempting to further improve
# version 0.0.4, 07 Oct 2025, further fixed effects improvement, see Francis email Friday, 19 September 2025
##############


# The code is written with the purpose of testing how different Bayesian modelling packages
# (greta, INLA, brms, nimble and rjags) can handle the same (relatively complex) 
# linear mixed model over increasing sample sizes. The initial asreml model is as follows (as per R agridat::barrero.maize):

# the asreml model, see R agridat::barrero.maize
##################
# m1 <- asreml::asreml(yield ~ env, data=dat,
#                random = ~ gen + gen:env,
#                residual = ~ dsum( ~ units|env),
#                workspace="500mb")
##################


rm(list = ls())
library(tidyverse)
library(agridat)
library(greta)
library(Matrix)
library(bayesplot)
library(coda)

# organise the dataset and create a directory for outputs
##################
data(barrero.maize)
data_initial <- barrero.maize
n_initial <- nrow(data_initial)
rm(barrero.maize)

# only retain necessary columns
data_initial <- data_initial %>% 
    select(-c("yor", "daystoflower", "plantheight", "earheight", "population", "lodged", "moisture", "testweight"))

# # if needed, expend the initial table x times
# ###############################
# dataset_expend <- FALSE
# if (dataset_expend) {
#     
#     x <- 100
#     
#     list_all_tabs <- list()
#     list_all_tabs[[1]] <- data_initial
#     ii0 <- 1
#     set.seed(55)
#     for (i0 in 2:x) {
# 
#         tab_new <- copy(data_initial)
#         tab_new[, year := year + 11*ii0]
#         tab_new[, env := paste0(yield, loc)]
#         tab_new[, yield := sample(yield, n_initial, replace = T)]
#         list_all_tabs[[i0]] <- tab_new
#         ii0 <- ii0 + 1
#     }
#     
#     data_initial_expended <- rbindlist(list_all_tabs)
#     data_initial <- data_initial_expended
#     
#     }


package_name <- "greta_index"

dir_out <- paste0("dir_packages_tesing/", package_name, "/")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)


# increase the memory limites to maximum
mem.maxVSize(vsize = Inf)
mem.maxNSize(nsize = Inf)

print(paste0("maximum sample size: ", nrow(data_initial)))

# 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010
# we will be taking the datasets 2000:2001, 2000:2002, ..., 2000:2010
years_vec <- sort(data_initial$year %>% unique)
list_computational_times <- list()
list_draws <- list()
ii1 <- 1

set.seed(55)
#rand_seeds_vec <- sample(10:1000, length(years_vec), replace = FALSE)
sample_sizes_vec <- 1:length(years_vec)*0

# run through 2000:2001, 2000:2002, ..., 2000:length(years_vec)
# for (i1 in 2:length(years_vec)) {
i1 = 2
#' FKCH: Using a fixed subset of years as an example to keep things simple. You can revert this back to a for loop and save results later...
    
    years_needed <- years_vec[1:i1]
    #seed_local <- rand_seeds_vec[i1]
    
    print(paste0("number of years to be processed: ", length(years_needed)))
    
    # drop missing outcomes, as required for greta (to be investigated later)
    data_work <- data_initial %>%
        filter(year %in% years_needed) %>% 
        na.omit()
    sample_sizes_vec[i1] <- nrow(data_work)
    
    # set the factors
    data_work <- data_work %>% 
        arrange(env) %>% 
        mutate(yearf = factor(year),
               env = factor(env),
               loc = factor(loc),
               gen = factor(gen),
               rep = factor(rep))

    data_work <- data_work %>% 
        mutate(loc_numeric = loc %>% factor %>% as.numeric,
               yearf_numeric = yearf %>% factor %>% as.numeric,
               env_numeric = env %>% factor %>% as.numeric,
               gen_numeric = gen %>% factor %>% as.numeric,
               rep_numeric = rep %>% factor %>% as.numeric,
               repenv_numeric = interaction(rep, env) %>% factor %>% as.numeric,
               genyear_numeric = interaction(gen, yearf) %>% factor %>% as.numeric,
               genloc_numeric = interaction(gen, loc) %>% factor %>% as.numeric,
               genenv_numeric = interaction(gen, env) %>% factor %>% as.numeric) 
    
    str(data_work)
    gc()
    
    ##--------------------
    # Model in asreml (benchmark)
    #' FKCH: Since I do not won asreml, then I asked Emi to run the (simpler) model and save the results on GitHub. Thanks to Emi!
    ##--------------------
    # m1_asreml <- asreml::asreml(yield ~ env, data = data_work,
    #                      random = ~ gen + gen:env,
    #                      residual = ~ dsum( ~ units|env),
    #                      workspace="500mb")
    
    m1_asreml_fixed <- readRDS(file = "../bayesreml/outputs/01-asreml-barrero-maize-m2-fixed.rds")
    m1_asreml_random <- readRDS(file = "../bayesreml/outputs/01-asreml-barrero-maize-m2-random.rds")
    m1_asreml_varcomp <- readRDS(file = "../bayesreml/outputs/01-asreml-barrero-maize-m2-vcomp.rds")
    
    
    ##--------------------
    # Start greta model
    ## Fixed effects component
    ##--------------------
    # computational time
    ptm1 <- proc.time() 
    
    # Define priors for fixed effects coefficients using index approach
    n_env <- table(data_work$env_numeric) %>% length

    # Intercept
    beta_intercept <- normal(0, 10)
    
    # Fixed effect for env (first level set as a reference)
    beta_env <- greta_array(dim = n_env)
    beta_env[1] <- 0
    beta_env[2:n_env] <- normal(0, 10, dim = n_env - 1)
    
    #' # Fixed effect for env, which is effectively the interaction between loc and yearf
    #' beta_env <- greta_array(dim = c(n_loc, n_yearf))
    #' beta_env[1,] <- 0
    #' beta_env[,1] <- 0
    #' beta_env[2:n_loc, 2:n_yearf] <- normal(0, 10, dim = (n_loc-1)*(n_yearf-1))
    #' 
    #' ## Check if any of interactions are missing in the data
    #' #' There are probably more precise ways to do this, but for now I will cheat and use the estimated lm model to work this out
    #' m1 <- lm(yield ~ loc*yearf, data = data_work)
    #' get_env_coefs <- m1$coefficients[-(1:(1+n_loc-1+n_yearf-1))] %>% 
    #'     matrix(ncol = n_yearf-1, byrow = FALSE)
    #' find_missing_interactions <- which(is.na(get_env_coefs), arr.ind = TRUE)
    #' if (nrow(find_missing_interactions) > 0) {
    #'     for (i2 in 1:nrow(find_missing_interactions)) {
    #'         beta_env[find_missing_interactions[i2,1]+1, find_missing_interactions[i2,2]+1] <- 0
    #'     }
    #' }
    #' rm(find_missing_interactions, get_env_coefs)
    
    
    # Build linear predictor for fixed part
    mu_fixed <- beta_intercept + beta_env[data_work$env_numeric]

    
    ##--------------------
    # Random effects component
    ##--------------------
    # Define priors for random effect standard deviations
    # Half-Cauchy priors are standard for variance components
    sigma_gen <- cauchy(0, 10, truncation = c(0, Inf))
    #sigma_rep_env  <- cauchy(0, 10, truncation = c(0, Inf))
    #sigma_gen_year <- cauchy(0, 10, truncation = c(0, Inf))
    #sigma_gen_loc  <- cauchy(0, 10, truncation = c(0, Inf))
    sigma_gen_env  <- cauchy(0, 10, truncation = c(0, Inf))
    
    # Define random effects
    # Each follows a normal distribution with mean 0 and respective sd
    u_gen <- normal(0, sigma_gen, dim = length(table(data_work$gen_numeric)))
    # u_rep_env <- normal(0, sigma_rep_env, dim = length(table(data_work$repenv_numeric)))
    # u_gen_year <- normal(0, sigma_gen_year, dim = length(table(data_work$genyear_numeric)))
    # u_gen_loc <- normal(0, sigma_gen_loc, dim = length(table(data_work$genloc_numeric)))
    u_gen_env <- normal(0, sigma_gen_env, dim = length(table(data_work$genenv_numeric)))
    
    ##--------------------
    # Outcome model
    ##--------------------
    # Calculate the linear predictor mu
    # Combines fixed and random effects using indices
    mu <- mu_fixed +
        u_gen[data_work$gen_numeric] +
        # u_rep_env[data_work$repenv_numeric] +
        # u_gen_year[data_work$genyear_numeric] +
        # u_gen_loc[data_work$genloc_numeric] +
        u_gen_env[data_work$genenv_numeric]

    # Half-Cauchy priors for environment-specific residual sd
    # Map residual sds to observations based on environment
    sigma_resid_env <- cauchy(0, 10, truncation = c(0, Inf), dim = length(unique(data_work$env_numeric)))
    sigma_resid_vec <- sigma_resid_env[data_work$env_numeric]

    # Observed outcome data to greta format
    y <- greta::as_data(data_work$yield)
    
    # Specify the likelihood: normal distribution with mean mu and sd sigma_resid_vec
    distribution(y) <- greta::normal(mu, sigma_resid_vec)
    
    
    ##--------------------
    ## Form model and sample...
    ##--------------------
    # Compile the greta model, tracking key parameters
    model_fit <- greta::model(beta_intercept,
                              beta_env, 
                              sigma_gen, 
                              #sigma_rep_env,
                              #sigma_gen_year, 
                              #sigma_gen_loc, 
                              sigma_gen_env,
                              sigma_resid_env,
                              u_gen,
                              u_gen_env)
                              
    
    # Run MCMC sampling for Bayesian inference (n_samples and warmup can be increased)
    m1 <- lm(yield ~ env, data = data_work)
    m1_greta <- greta::mcmc(model_fit, 
                            n_samples = 1000, 
                            warmup = 1000, 
                            initial_values = initials(beta_intercept = m1$coefficients[1]),
                            chains = 4, 
                            n_cores = 4)
    
    
    ##--------------------
    # Examine results -- using the posterior median and quantiles from greta as a basic check to start off with
    ##--------------------
    get_quantiles <- summary(m1_greta)$quantiles %>% 
        as.data.frame()
    rownames(get_quantiles) <- c("Intercept",
                                 paste0("env", levels(data_work$env)), 
                                 "sigma_variancecom_gen",
                                 "sigma_variancecom_gen_env",
                                 # paste0("loc", rep(levels(data_work$loc), n_yearf), 
                                 #        ":yearf", rep(levels(data_work$yearf), each = n_loc)),
                                 paste0("sigma_residual_env", levels(data_work$env)),
                                 paste0("gen_", levels(data_work$gen)),
                                 interaction(paste0("gen_", data_work$gen), paste0("env_", data_work$env), sep = ":") %>% factor %>% levels())
    
    
    ## Fixed effects
    data.frame(asreml = m1_asreml_fixed$estimate %>% setNames(m1_asreml_fixed$term),
               greta = get_quantiles[1:(1+length(levels(data_work$env))), "50%"])

    
    ## Variance components
    data.frame(asreml = m1_asreml_varcomp$estimate %>% setNames(m1_asreml_varcomp$term),
               greta = get_quantiles^2 %>% filter(str_detect(rownames(get_quantiles), "sigma_")) %>% pull(`50%`))
    
    
    ## Random effects 
    #' FKCH: Note this is slightly more finnicky to compare as the order of the levels in the interaction differ between asreml and greta
    data.frame(asreml = m1_asreml_random %>% 
                   filter(estimate != 0) %>% 
                   select(term, estimate) %>%
                   dplyr::slice(1:length(table(data_work$gen)))
               ) %>% 
        mutate(greta = get_quantiles %>% 
                   filter(str_detect(rownames(get_quantiles), "^gen_")) %>% 
                   dplyr::slice(1:length(table(data_work$gen))) %>% 
                   pull(`50%`))
    
    left_join(m1_asreml_random %>% 
                   filter(estimate != 0) %>% 
                   select(term, estimate) %>%
                   dplyr::slice(-(1:length(table(data_work$gen)))),
              get_quantiles %>% 
                  filter(str_detect(rownames(get_quantiles), "^gen_")) %>% 
                  dplyr::slice(-(1:length(table(data_work$gen)))) %>% 
                  rownames_to_column(var = "term") %>%
                  select(term, `50%`),
              by = "term") %>% 
        print(n = Inf)
    
    
    
    exec_time1 <- proc.time() - ptm1
    list_computational_times[[ii1]] <- exec_time1
    ii1 <- ii1+1
    
    print(exec_time1)
    
    # save the results
    file_names <- paste0(dir_out, "r_out_04a_greta_test_", i1, "_2000_",  years_vec[i1] ,".RData")
    #save(m1_greta, exec_time1, years_vec, seed_local, file = file_names)
    
    rm(model_fit, m1_greta, y, u_gen_year, u_gen_loc, u_gen_env)
    
    
# } # ends for (i1 in 2:s1)
