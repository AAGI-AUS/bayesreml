# This is the code updated based on comments from Francis (email Friday, 19 September 2025), coding fixed effects in a different way
# The code is written with the purpose of testing how different Bayesian modelling packages (greta, INLA, brms, nimble and runjags) can handle the same (relatively complex) linear mixed model over increasing sample sizes. The initial asreml model is as follows (as per R agridat::barrero.maize):

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

# 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010
# we will be taking the datasets 2000:2001, 2000:2002, ..., 2000:2010
years_vec <- sort(data_initial$year %>% unique)
list_computational_times <- list()
sample_sizes_vec <- NULL
do_summaries <- FALSE


# run through 2000:2001, 2000:2002, ..., 2000:length(years_vec)
for (i1 in 2:length(years_vec)) {
    years_needed <- years_vec[1:i1]
    print(paste0("number of years to be processed: ", length(years_needed)))
    
    # drop missing outcomes, although note jags can automatically handle this by either omitting or imputing if a distribution is assumed
    data_work <- data_initial %>%
        filter(year %in% years_needed) %>% 
        na.omit()
    
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
    
    # str(data_work)
    
    ##--------------------
    # Model in asreml (benchmark)
    #' FKCH: Since I do not won asreml, then I asked Emi to run the (simpler) model and save the results on GitHub. Thanks to Emi!
    ##--------------------
    # m2_asreml <- asreml::asreml(yield ~ env, data = data_work,
    #                      random = ~ gen + gen:env,
    #                      residual = ~ dsum( ~ units|env),
    #                      workspace="500mb")
    # m2_asreml_fixed <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-fixed.rds"))
    # m2_asreml_random <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-random.rds"))
    # m2_asreml_varcomp <- readRDS(file = here("outputs", "01-asreml-barrero-maize-m2-vcomp.rds"))
    
    
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
    # Note default settings of burnin = 4000, sample = 10000, adapt = 1000 MCMC samples are used
    m2_jags <- run.jags(model = model, 
                        n.chains = 4, 
                        method = "parallel")
    
    
    ##--------------------
    # Examine results -- using the posterior median and quantiles from jags as a basic check to start off with
    ##--------------------
    if(do_summaries){
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
        }
    
    ##--------------------
    # Final collection of computation times
    ##--------------------
    exec_time1 <- proc.time() - ptm1
    sample_sizes_vec <- c(sample_sizes_vec, nrow(data_work))
    list_computational_times[[i1-1]] <- exec_time1

    print(exec_time1)
    rm(m2_jags, y, n_env, n_gen, n_genenv, N, env_numeric, gen_numeric, genenv_numeric, data_work)
    gc()

    } # ends for (i1 in 2:s1)

    

tab_times_summary <- data.frame(sample_size = sample_sizes_vec, 
                                times_sec = sapply(list_computational_times, function(x) x[3]))
writexl::write_xlsx(tab_times_summary, path = "m2testing/xlsx_tab_runjags_times_m2.xlsx")


##---------------------
sessioninfo::session_info()
##---------------------
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.3 (2024-02-29)
# os       Linux Mint 22.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2025-10-23
# rstudio  2025.05.1+513 Mariposa Orchid (desktop)
# pandoc   3.1.3 @ /usr/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version  date (UTC) lib source
# agridat     * 1.24     2024-10-27 [1] CRAN (R 4.3.3)
# bayesplot   * 1.11.1   2024-02-15 [3] CRAN (R 4.3.3)
# cli           3.6.5    2025-04-23 [1] CRAN (R 4.3.3)
# coda        * 0.19-4.1 2024-01-31 [3] CRAN (R 4.3.2)
# colorspace    2.1-1    2024-07-26 [1] CRAN (R 4.3.3)
# dplyr       * 1.1.4    2023-11-17 [1] CRAN (R 4.3.3)
# forcats     * 1.0.0    2023-01-29 [3] CRAN (R 4.2.2)
# generics      0.1.4    2025-05-09 [1] CRAN (R 4.3.3)
# ggplot2     * 3.5.1    2024-04-23 [1] CRAN (R 4.3.3)
# glue          1.8.0    2024-09-30 [1] CRAN (R 4.3.3)
# gtable        0.3.6    2024-10-25 [1] CRAN (R 4.3.3)
# here        * 1.0.2    2025-09-15 [1] CRAN (R 4.3.3)
# hms           1.1.3    2023-03-21 [3] CRAN (R 4.3.1)
# lattice       0.22-5   2023-10-24 [4] CRAN (R 4.3.3)
# lifecycle     1.0.4    2023-11-07 [1] CRAN (R 4.3.3)
# lubridate   * 1.9.3    2023-09-27 [3] CRAN (R 4.3.1)
# magrittr      2.0.4    2025-09-12 [1] CRAN (R 4.3.3)
# Matrix      * 1.6-5    2024-01-11 [4] CRAN (R 4.3.2)
# munsell       0.5.1    2024-04-01 [1] CRAN (R 4.3.3)
# pillar        1.11.1   2025-09-17 [1] CRAN (R 4.3.3)
# pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 4.3.3)
# purrr       * 1.0.4    2025-02-05 [1] CRAN (R 4.3.3)
# R6            2.6.1    2025-02-15 [1] CRAN (R 4.3.3)
# readr       * 2.1.5    2024-01-10 [3] CRAN (R 4.3.2)
# RhpcBLASctl   0.23-42  2023-02-11 [3] CRAN (R 4.3.2)
# rlang         1.1.6    2025-04-11 [1] CRAN (R 4.3.3)
# rprojroot     2.1.1    2025-08-26 [1] CRAN (R 4.3.3)
# rstudioapi    0.17.1   2024-10-22 [1] CRAN (R 4.3.3)
# runjags     * 2.2.2-5  2025-04-09 [1] CRAN (R 4.3.3)
# scales        1.3.0    2023-11-28 [1] CRAN (R 4.3.3)
# sessioninfo   1.2.2    2021-12-06 [3] CRAN (R 4.1.2)
# stringi       1.8.4    2024-05-06 [1] CRAN (R 4.3.3)
# stringr     * 1.5.1    2023-11-14 [1] CRAN (R 4.3.3)
# tibble      * 3.3.0    2025-06-08 [1] CRAN (R 4.3.3)
# tidyr       * 1.3.1    2024-01-24 [1] CRAN (R 4.3.3)
# tidyselect    1.2.1    2024-03-11 [1] CRAN (R 4.3.3)
# tidyverse   * 2.0.0    2023-02-22 [3] CRAN (R 4.3.1)
# timechange    0.3.0    2024-01-18 [3] CRAN (R 4.3.2)
# tzdb          0.4.0    2023-05-12 [3] CRAN (R 4.3.1)
# vctrs         0.6.5    2023-12-01 [1] CRAN (R 4.3.3)
# withr         3.0.2    2024-10-28 [1] CRAN (R 4.3.3)
# writexl       1.5.0    2024-02-09 [3] CRAN (R 4.3.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.3
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    