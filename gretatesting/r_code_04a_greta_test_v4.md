This is the code updated based on comments from Francis (email Friday, 19 September 2025), coding fixed effects in adifferent way.

One more minor improvement made is replacing `as.numeric` with `as.integer` for factors.

From eyeballing outputs, `summary(m1_greta)` estimates still differ from the initial "model.matrix" approach.

```
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
# m1 <- asreml::asreml(yield ~ loc + yearf + env, data=dat,
#                random = ~ gen + rep:env +
#                  gen:yearf + gen:loc +
#                  gen:env,
#                residual = ~ dsum( ~ units|env),
#                workspace="500mb")
##################

rm(list = ls())
#library(tidyverse)
library(data.table)
library(agridat)
library(greta)
library(Matrix)


# organise the dataset and create a directory for ouputs
##################
data(barrero.maize)
data_initial <- data.table(barrero.maize)
n_initial <- nrow(data_initial)

# only retain necessary columns
cols_to_drop <- c("yor", "daystoflower", "plantheight", "earheight", "population", "lodged", "moisture", "testweight")
data_initial[, (cols_to_drop) := NULL]

# if needed, expend the initial table x times
###############################
dataset_expend <- FALSE

if (dataset_expend) {

	x <- 100

	list_all_tabs <- list()
	list_all_tabs[[1]] <- data_initial
	ii0 <- 1
	set.seed(55)
	for (i0 in 2:x) {

		tab_new <- copy(data_initial)
		tab_new[, year := year + 11*ii0]
		tab_new[, env := paste0(yield, loc)]
		tab_new[, yield := sample(yield, n_initial, replace = T)]
		list_all_tabs[[i0]] <- tab_new
		ii0 <- ii0 + 1

	}

	data_initial_expended <- rbindlist(list_all_tabs)
	data_initial <- data_initial_expended

}
###############################


package_name <- "greta_index"

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
sample_sizes_vec <- 1:s1*0

# run through 2000:2001, 2000:2002, ..., 2000:s1
for (i1 in 2:s1) {
#for (i1 in s1) {
    years_needed <- years_vec[1:i1]
    data_work <- copy(data_initial)
    seed_local <- rand_seeds_vec[i1]
    
    print(paste0("number of years to be processed: ", length(years_needed)))
    
    jj <- data_work[, !is.na(match(year, years_needed))]
    data_work <- data_work[jj,]
    
    # drop missing outcomes, as required for greta (to be investigated later)
    data_work <- data_work[!is.na(yield),]
    sample_sizes_vec[i1] <- nrow(data_work)
    
    # set the factors and their integer representations
    data_work <- data_work[order(env),]
    data_work[, yearf := factor(year)][, yearf_numeric := as.integer(yearf)]
    data_work[, env := factor(env)][, env_numeric := as.integer(env)]
    data_work[, loc := factor(loc)][, loc_numeric := as.integer(loc)]
    data_work[, gen := factor(gen)][, gen_numeric := as.integer(gen)]
    data_work[, rep := factor(rep)][, rep_numeric := as.integer(rep)]
    data_work[, repenv_numeric := as.integer(factor(interaction(rep, env)))]
    data_work[, genyear_numeric := as.integer(factor(interaction(gen, yearf)))]
    data_work[, genloc_numeric := as.integer(factor(interaction(gen, loc)))]
    data_work[, genenv_numeric := as.integer(factor(interaction(gen, env)))]
    
       
    gc()
    
    # computational time
    ptm1 <- proc.time() 
    
    # Define priors for fixed effects coefficients using index approach
	n_loc <- max(data_work$loc_numeric)
	n_yearf <- max(data_work$yearf_numeric)
	n_env <- max(data_work$env_numeric)
		
	# Intercept
	beta_intercept <- normal(0, 100)

	# Fixed effect for loc (first level set as a reference)
	beta_loc_rest <- normal(0, 100, dim = n_loc - 1)
	beta_loc <- greta_array(dim = n_loc)
	beta_loc[1] <- 0
	beta_loc[2:n_loc] <- beta_loc_rest

	# Fixed effect for yearf (first level set as a reference)
	beta_yearf_rest <- normal(0, 100, dim = n_yearf - 1)
	beta_yearf <- greta_array(dim = n_yearf)
	beta_yearf[1] <- 0
	beta_yearf[2:n_yearf] <- beta_yearf_rest

	# Fixed effect for env (first level set as a reference)
	beta_env_rest <- normal(0, 100, dim = n_env - 1)
	beta_env <- greta_array(dim = n_env)
	beta_env[1] <- 0
	beta_env[2:n_env] <- beta_env_rest

	# Build linear predictor for fixed part
	mu_fixed <- beta_intercept +
				beta_loc[data_work$loc_numeric] +
				beta_yearf[data_work$yearf_numeric] +
				beta_env[data_work$env_numeric]
	
	
	# Impose sum-to-zero constraints (??)
	# beta_loc_raw <- normal(0, 100, dim = n_loc)
# 	beta_loc <- beta_loc_raw - mean(beta_loc_raw)
# 	
# 	beta_yearf_raw <- normal(0, 100, dim = n_yearf)
# 	beta_yearf <- beta_yearf_raw - mean(beta_yearf_raw)
# 	
# 	beta_env_raw <- normal(0, 100, dim = n_env)
# 	beta_env <- beta_env_raw - mean(beta_env_raw)

    # Define priors for random effect standard deviations
    # Half-Cauchy priors are standard for variance components
    sigma_gen      <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_rep_env  <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_year <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_loc  <- cauchy(0, 2, truncation = c(0, Inf))
    sigma_gen_env  <- cauchy(0, 2, truncation = c(0, Inf))
    
    # Define random effects using index approach
	n_gen <- max(data_work$gen_numeric)
	n_repenv <- max(data_work$repenv_numeric)
	n_genyear <- max(data_work$genyear_numeric)
	n_genloc <- max(data_work$genloc_numeric)
	n_genenv <- max(data_work$genenv_numeric)
	u_gen <- normal(0, sigma_gen, dim = n_gen)
	u_rep_env <- normal(0, sigma_rep_env, dim = n_repenv)
	u_gen_year <- normal(0, sigma_gen_year, dim = n_genyear)
	u_gen_loc <- normal(0, sigma_gen_loc, dim = n_genloc)
	u_gen_env <- normal(0, sigma_gen_env, dim = n_genenv)

	# Define residual standard deviations for each environment
	# Half-Cauchy priors for environment-specific residual sd
	sigma_resid_env <- cauchy(0, 2, truncation = c(0, Inf), dim = n_env)
	sigma_resid_vec <- sigma_resid_env[data_work$env_numeric]
    
        
    # Calculate the linear predictor mu
	# Combines fixed and random effects using indices
	mu <- mu_fixed +
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
    # model_fit <- greta::model(beta_fixef, sigma_gen, sigma_rep_env,
#                               sigma_gen_year, sigma_gen_loc, sigma_gen_env,
#                               sigma_resid_env)
                              
    model_fit <- greta::model(beta_loc, beta_yearf, beta_env, sigma_gen, sigma_rep_env,
                          sigma_gen_year, sigma_gen_loc, sigma_gen_env,
                          sigma_resid_env)
    
    # Run MCMC sampling for Bayesian inference (n_samples and warmup can be incresed)
    set.seed(seed_local)
    m1_greta <- greta::mcmc(model_fit, n_samples = 1000, warmup = 500, chains = 4, n_cores = 10)
    
    exec_time1 <- proc.time() - ptm1
    list_computational_times[[ii1]] <- exec_time1
    ii1 <- ii1+1
    
    print(exec_time1)
    
    # save the results
    file_names <- paste0(dir_out, "r_out_04a_greta_test_", i1, "_2000_",  years_vec[i1] ,".RData")
    #save(m1_greta, exec_time1, years_vec, seed_local, file = file_names)
    
    rm(model_fit, m1_greta, y, u_gen_year, u_gen_loc, u_gen_env)
    
    
} # ends for (i1 in 2:s1)

```
