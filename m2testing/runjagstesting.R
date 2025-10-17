# This is the code updated based on comments from Francis (email Friday, 19 September 2025), coding fixed effects in a different way
# The code is written with the purpose of testing how different Bayesian modelling packages
# (greta, INLA, brms, nimble and runjags) can handle the same (relatively complex) 
# linear mixed model over increasing sample sizes. The initial asreml model is as follows (as per R agridat::barrero.maize):

# the asreml model, see R agridat::barrero.maize
##################
# dat <- agridat::barrero.maize |> 
#     filter(year %in% 2000:2001) |> 
#     mutate(yearf = factor(year)) |> 
#     select(-c("yor", "daystoflower", "plantheight", "earheight", "population", "lodged", "moisture", "testweight")) |> 
#     na.omit()
# 
# m2 <- asreml(
#     yield ~ env,
#     data = dat,
#     random = ~ gen + gen:env,
#     residual = ~ dsum(~ units | env),
#     workspace = "500mb"
# )
##################
rm(list = ls())
library(tidyverse)
library(agridat)
library(runjags)
library(Matrix)
library(bayesplot)
library(coda)
here::i_am("m2testing/runjagstesting.R")
library(here)

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


# package_name <- "jags_index"
# 
# dir_out <- paste0("dir_packages_tesing/", package_name, "/")
#     data_initial_expended <- rbindlist(list_all_tabs)
#     data_initial <- data_initial_expended
#     
#     }


# package_name <- "jags_index"
# 
# dir_out <- paste0("dir_packages_tesing/", package_name, "/")
# dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
# 
# 
# # increase the memory limites to maximum
# mem.maxVSize(vsize = Inf)
# mem.maxNSize(nsize = Inf)
# 
# print(paste0("maximum sample size: ", nrow(data_initial)))

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
    
    # drop missing outcomes, although note jags can automatically handle this by either omitting or imputing if a distribution is assumed
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
    # m2_asreml <- asreml::asreml(yield ~ env, data = data_work,
    #                      random = ~ gen + gen:env,
    #                      residual = ~ dsum( ~ units|env),
    #                      workspace="500mb")
    
    m2_asreml_fixed <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-fixed.rds"))
    m2_asreml_random <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-random.rds"))
    m2_asreml_varcomp <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-vcomp.rds"))
    
    
    ##--------------------
    # Start runjags model
    ##--------------------
    # computational time
    ptm1 <- proc.time() 

    env_numeric <- data_work$env_numeric
    gen_numeric <- data_work$gen_numeric
    genenv_numeric <- data_work$genenv_numeric
    n_env <- table(env_numeric) %>% length
    n_gen <- table(gen_numeric) %>% length
    n_genenv <- table(genenv_numeric) %>% length
    N <- nrow(data_work)
    y <- data_work$yield
    
    model <- "model {
        # Fixed effects component
        beta_intercept ~ dnorm(0, 0.01)

        # Fixed effect for env (first level set as a reference)
        # The first level is fixed to 0 (reference level).
        beta_env[1] <- 0

        for (j in 2:n_env) { #data# n_env
              beta_env[j] ~ dnorm(0, 0.01)
              }        
    
        # Random effects component
        # Define priors for random effect standard deviations
        # Half-Cauchy priors are standard for variance components
        sigma_gen ~ dt(0, 0.01, 1)T(0,)
        #sigma_rep_env ~ dt(0, 0.01, 1)T(0,)
        #sigma_gen_year ~ dt(0, 0.01, 1)T(0,)
        #sigma_gen_loc ~ dt(0, 0.01, 1)T(0,)
        sigma_gen_env ~ dt(0, 0.01, 1)T(0,)
    
    
        # Define random effects
        for (i in 1:n_gen) { #data# n_gen
            u_gen[i] ~ dnorm(0, pow(sigma_gen, -2))
            }
        # u_rep_env ~ normal(0, sigma_rep_env, dim = length(table(data_work$repenv_numeric)))
        # u_gen_year ~ normal(0, sigma_gen_year, dim = length(table(data_work$genyear_numeric)))
        # u_gen_loc ~ normal(0, sigma_gen_loc, dim = length(table(data_work$genloc_numeric)))
        for (i in 1:n_genenv) { #data# n_genenv
            u_gen_env[i] ~ dnorm(0, pow(sigma_gen_env, -2))
            }
        
        
        # Outcome model
        # Calculate the linear predictor mu
        # Combines fixed and random effects using indices
        for(i in 1:N) { #data# N
            mu[i] <- beta_intercept + beta_env[env_numeric[i]] + u_gen[gen_numeric[i]] + u_gen_env[genenv_numeric[i]] 
            #data# env_numeric, gen_numeric, genenv_numeric
            }  
        # u_rep_env[data_work$repenv_numeric] +
        # u_gen_year[data_work$genyear_numeric] +
        # u_gen_loc[data_work$genloc_numeric] +

    
        # Half-Cauchy priors for environment-specific residual sd
        # Map residual sds to observations based on environment
        for(j in 1:n_env) { 
            sigma_resid_env[j] ~ dt(0, 0.01, 1)T(0,)
            tau_resid_env[j] <- pow(sigma_resid_env[j], -2)
            }
        for(i in 1:N) {
            y[i] ~ dnorm(mu[i], tau_resid_env[env_numeric[i]]) #data# y
            }

    
    #monitor# beta_intercept, beta_env, sigma_gen, sigma_gen_env, sigma_resid_env, u_gen, u_gen_env
    }"
    

    ##--------------------
    ## Form model and sample...
    ##--------------------
    # Run MCMC sampling for Bayesian inference (burnin, sample, and adapt are set to defaults; these are not necessarily compatible with what brms, greta and other MCMC samples use!)
    m2_jags <- run.jags(model = model, 
                        n.chains = 4, 
                        method = "parallel")
    
    
    ##--------------------
    # Examine results -- using the posterior median and quantiles from jags as a basic check to start off with
    ##--------------------
    s <- summary(m2_jags) #' Note this command takes a while because runjags does many of summaries. There is probably a way to produce this bad or use CODA instead, but this is a future enhancement issue
    get_quantiles <- s %>% 
        as.data.frame() %>% 
        select(Lower95:Upper95) 
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
    data.frame(asreml = m2_asreml_fixed$estimate %>% setNames(m2_asreml_fixed$term),
               jags = get_quantiles[1:(1+length(levels(data_work$env))),] %>% pull("Median"))

    
    ## Variance components
    data.frame(asreml = m2_asreml_varcomp$estimate %>% setNames(m2_asreml_varcomp$term),
               jags = get_quantiles^2 %>% filter(str_detect(rownames(get_quantiles), "sigma_")) %>% pull("Median"))
    
    
    ## Random effects 
    #' FKCH: Note this is slightly more finnicky to compare as the order of the levels in the interaction differ between asreml and jags
    data.frame(asreml = m2_asreml_random %>% 
                   filter(estimate != 0) %>% 
                   select(term, estimate) %>%
                   dplyr::slice(1:length(table(data_work$gen)))
               ) %>% 
        mutate(jags = get_quantiles %>% 
                   filter(str_detect(rownames(get_quantiles), "^gen_")) %>% 
                   dplyr::slice(1:length(table(data_work$gen))) %>% 
                   pull("Median"))
    
    left_join(m2_asreml_random %>% 
                   filter(estimate != 0) %>% 
                   select(term, estimate) %>%
                   dplyr::slice(-(1:length(table(data_work$gen)))),
              get_quantiles %>% 
                  filter(str_detect(rownames(get_quantiles), "^gen_")) %>% 
                  dplyr::slice(-(1:length(table(data_work$gen)))) %>% 
                  rownames_to_column(var = "term") %>%
                  select(term, `Median`),
              by = "term") %>% 
        print(n = Inf)
    
    
    
    exec_time1 <- proc.time() - ptm1
    list_computational_times[[ii1]] <- exec_time1
    ii1 <- ii1+1
    
    print(exec_time1)
    
    # save the results
    file_names <- paste0(dir_out, "r_out_04a_jags_test_", i1, "_2000_",  years_vec[i1] ,".RData")
    #save(m1_jags, exec_time1, years_vec, seed_local, file = file_names)
    
    rm(model_fit, m1_jags, y, u_gen_year, u_gen_loc, u_gen_env)
    
    
# } # ends for (i1 in 2:s1)

    

    