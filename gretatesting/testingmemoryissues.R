# source("r_code_04a_greta_test_v1.R")

##############
# Author: Max Moldovan
# max.moldovan@adelaide.edu.au
# version 0.0.1, 12 Sep 2025
##############


# The code is written with the purpose of testing how different Bayesian modelling packages
# (greta, INLA, brms, nimble and rjags) can handle the same (relatively complex) 
# linear mixed model over increasing sample sizes. The initial asreml model is as follows (as per R agridat::barrero.maize):

# the asreml model, see R agridat::barrero.maize
##################
# m1 <- asreml::asreml(yield ~ loc + yearf + env, data=dat,
#                random = ~ gen + rep:env +
#                  gen:yearf + gen:loc +
#                  gen:env,
#                residual = ~ dsum( ~ units|env),
#                workspace="500mb")
##################

rm(list = ls())
library(tidyverse)
library(data.table)
library(agridat)
library(greta)
library(Matrix)


# organise the dataset and create a directory for ouputs
##################
data(barrero.maize)
data_initial <- data.table(barrero.maize)

package_name <- "greta"

dir_out <- paste0("dir_packages_tesing/", package_name, "/")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)


# increase the memory limites to maximum
mem.maxVSize(vsize = Inf)
mem.maxNSize(nsize = Inf)

print(paste0("maximum sample size: ", nrow(data_initial)))

# 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010
# we will be taking the datasets 2000:2001, 2000:2002, ..., 2000:2010
years_vec <- sort(data_initial[,unique(year)])
list_computational_times <- list()
list_draws <- list()
ii1 <- 1


s1 <- length(years_vec)
set.seed(55)
rand_seeds_vec <- sample(10:1000, s1, replace = F)

# run through 2000:2001, 2000:2002, ..., 2000:2010
for (i1 in 2:s1) {
    years_needed <- years_vec[1:i1]
    data_work <- copy(data_initial)
    seed_local <- rand_seeds_vec[i1]
    
    print(years_needed)
    
    jj <- data_work[, !is.na(match(year, years_needed))]
    data_work <- data_work[jj,]
    
    # drop missing outcomes, as required for greta (to be investigated later)
    data_work <- data_work[!is.na(yield),]
    
    # set the factors
    data_work <- data_work[order(env),]
    data_work[, yearf := factor(year)]
    data_work[, env := factor(env)]
    data_work[, loc := factor(loc)]
    data_work[, gen := factor(gen)]
    data_work[, rep := factor(rep)]
    
    data_work$loc %>% table
    data_work$loc %>% table %>% length
    data_work$yearf %>% table
    data_work$yearf %>% table %>% length
    data_work$env %>% table
    data_work$env %>% table %>% length
    data_work$gen %>% table
    data_work$gen %>% unique %>% length
    data_work$rep %>% table
    data_work$rep %>% unique %>% length
    
    data_work <- as.data.frame(data_work) %>% 
        mutate(loc_numeric = loc %>% factor %>% as.numeric,
               yearf_numeric = yearf %>% factor %>% as.numeric,
               env_numeric = env %>% factor %>% as.numeric,
               gen_numeric = gen %>% factor %>% as.numeric,
               rep_numeric = rep %>% factor %>% as.numeric,
               repenv_numeric = interaction(rep, env) %>% factor %>% as.numeric,
               genyear_numeric = interaction(gen, yearf) %>% factor %>% as.numeric,
               genloc_numeric = interaction(gen, loc) %>% factor %>% as.numeric,
               genenv_numeric = interaction(gen, env) %>% factor %>% as.numeric) 
    
    gc()
    
    # computational time
    ptm1 <- proc.time() 
    
    # Create design matrix for fixed effects
    X <- model.matrix(~ loc + yearf + env, data = data_work)
    
    # Random effects
    # Z_rep_env  <- model.matrix(~ 0 + rep:env, data = data_work)
    # Z_gen_year <- model.matrix(~ 0 + gen:yearf, data = data_work)
    # Z_gen_loc  <- model.matrix(~ 0 + gen:loc, data = data_work)
    # Z_gen_env  <- model.matrix(~ 0 + gen:env, data = data_work)
    
    # Is there a way to make greta to work with sparse matrices
    # (see https://forum.greta-stats.org/t/sparse-matrix-format-for-count-and-binary-data-setting-dtype-of-tensors/124)
    # Sparse matrices CANNOT be passed to Tensorflow at this point (12 Sep 2025)
    # 
    # Z_gen      <- sparse.model.matrix(~ 0 + gen, data = data_work)
    # Z_rep_env  <- sparse.model.matrix(~ 0 + rep:env, data = data_work)
    # Z_gen_year <- sparse.model.matrix(~ 0 + gen:yearf, data = data_work)
    # Z_gen_loc  <- sparse.model.matrix(~ 0 + gen:loc, data = data_work)
    # Z_gen_env  <- sparse.model.matrix(~ 0 + gen:env, data = data_work)
    
    # Define priors for fixed effects coefficients
    beta_fixef <- normal(0, 100, dim = ncol(X))
    
    # Define priors for random effect standard deviations
    # Half-Cauchy priors are standard for variance components
    sigma_gen      <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_rep_env  <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_year <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_loc  <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_env  <- cauchy(0, 2, truncation = c(0, Inf))
    
    # Define random effects
    # Each follows a normal distribution with mean 0 and respective sd
    u_gen      <- normal(0, sigma_gen, dim = length(unique(data_work$gen_numeric)))
    u_rep_env  <- normal(0, sigma_rep_env, dim = length(unique(data_work$repenv_numeric)))
    u_gen_year <- normal(0, sigma_gen_year, dim = length(unique(data_work$genyear_numeric)))
    u_gen_loc  <- normal(0, sigma_gen_loc, dim = length(unique(data_work$genloc_numeric)))
    u_gen_env  <- normal(0, sigma_gen_env, dim = length(unique(data_work$genenv_numeric)))
    
    # Define residual standard deviations for each environment
    # Half-Cauchy priors for environment-specific residual sd
    env_levels <- unique(data_work$env)
    n_env <- length(env_levels)
    
    # Map residual sds to observations based on environment
    sigma_resid_env <- cauchy(0, 2, truncation = c(0, Inf), dim = n_env)
    env_index <- as.integer(data_work$env)
    sigma_resid_vec <- sigma_resid_env[env_index]
    
    
    # Calculate the linear predictor mu
    # Combines fixed and random effects
    mu <- X %*% beta_fixef +
        u_gen[data_work$gen_numeric] +
        u_rep_env[data_work$repenv_numeric] +
        u_gen_year[data_work$genyear_numeric] +
        u_gen_loc[data_work$genloc_numeric] +
        u_gen_env[data_work$genenv_numeric]
    
    # Observed outcome data to greta format
    y <- greta::as_data(data_work$yield)
    
    # Specify the likelihood: normal distribution with mean mu and sd sigma_resid_vec
    distribution(y) <- greta::normal(mu, sigma_resid_vec)
    
    # Compile the greta model, tracking key parameters
    model_fit <- greta::model(beta_fixef, sigma_gen, sigma_rep_env,
                              sigma_gen_year, sigma_gen_loc, sigma_gen_env,
                              sigma_resid_env)
    
    # Run MCMC sampling for Bayesian inference (n_samples and warmup can be incresed)
    set.seed(seed_local)
    m1_greta <- greta::mcmc(model_fit, n_samples = 1000, warmup = 500, chains = 4, n_cores = 6)
    
    exec_time1 <- proc.time() - ptm1
    list_computational_times[[ii1]] <- exec_time1
    ii1 <- ii1+1
    
    # save the results
    file_names <- paste0(dir_out, "r_out_04a_greta_test_", i1, "_2000_",  years_vec[i1] ,".RData")
    save(m1_greta, exec_time1, years_vec, seed_local, file = file_names)
    
    rm(model_fit, m1_greta, X, y, Z_gen, u_gen, Z_rep_env, u_rep_env, Z_gen_year, 
       u_gen_year, Z_gen_loc, u_gen_loc, Z_gen_env, u_gen_env)
    
    
} # ends for (i1 in 2:s1)
