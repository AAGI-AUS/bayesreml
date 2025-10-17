# source("r_code_04d_nimble_test_v2.R")

##############
# Author: Max Moldovan
# max.moldovan@adelaide.edu.au
# version 0.0.1, 12 Sep 2025
# version 0.0.2, 16 Oct 2025 reduced asreml m2 model, see the same date emails
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



require(data.table)
require(agridat)
require(nimble)
require(fst)


rm(list = ls())


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
dataset_expend <- F

if (dataset_expend) {

	x <- 25

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

package_name <- "nimblev2"

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
#for (i1 in 2:s1) {
for (i1 in 2) {

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
	
	gc()
	
	# computational time
	ptm1 <- proc.time() 

	# Create unique identifiers for the random effects
	data_work[, rep_env := interaction(rep, env, drop = TRUE)]
	data_work[, gen_yearf := interaction(gen, yearf, drop = TRUE)]
	#data_work[, gen_loc := interaction(gen, loc, drop = TRUE)]
	data_work[, gen_env := interaction(gen, env, drop = TRUE)]

	# Prepare data for NIMBLE
	nimble_data <- list(yield = data_work$yield)

	# Constants for the model
	constants <- list(
	  n_obs = nrow(data_work),
	  env = as.numeric(data_work$env),
	  gen = as.numeric(data_work$gen),
	  gen_env = as.numeric(factor(data_work$gen_env)),
	  n_env = nlevels(data_work$env),
	  n_gen = nlevels(data_work$gen),
	  n_gen_env = nlevels(factor(data_work$gen_env))
	)

	# Check if any constants are zero
	if(any(unlist(constants) == 0)) {
	  warning("Some constants are zero. Adjusting model structure...")
  
	  # Remove components with zero levels
	  constants <- constants[!sapply(constants, function(x) length(x) == 1 && x == 0)]
	}

	# NIMBLE code with intercept
	model_code <- nimbleCode({
	  # Fixed effects with intercept
	  for (i in 1:n_obs) {
		mu[i] <- beta_0 + beta_env[env[i]] +
				 u_gen[gen[i]] + u_gen_env[gen_env[i]]
		yield[i] ~ dnorm(mu[i], tau_env[env[i]])
	  }
	  # Priors for fixed effects
	  beta_0 ~ dnorm(0, sd = 100)  # Intercept prior (broad, as in prior implementations)
	  beta_env[1] <- 0  # Fix first level for identifiability
	  for (i in 2:n_env) {
		beta_env[i] ~ dnorm(0, sd = 100)  # Relaxed prior (1e-6 precision was too tight)
	  }
	  
	  # Random effects priors
	  for (i in 1:n_gen) {
		u_gen[i] ~ dnorm(0, tau_gen)
	  }
	  for (i in 1:n_gen_env) {
		u_gen_env[i] ~ dnorm(0, tau_gen_env)
	  }
	  # Variance components (precision parameters)
	  tau_gen ~ dgamma(0.001, 0.001)
	  tau_gen_env ~ dgamma(0.001, 0.001)
	  
	  # Environment-specific residual variances
	  for (i in 1:n_env) {
		tau_env[i] ~ dgamma(0.001, 0.001)
		sigma_env[i] <- 1 / sqrt(tau_env[i])
	  }
	  
	  # Convert precision to standard deviations
	  sigma_gen <- 1 / sqrt(tau_gen)
	  sigma_gen_env <- 1 / sqrt(tau_gen_env)
	})

	# Initial values
	inits <- list(
	  beta_0 = 0,
	  beta_env = c(0, if (constants$n_env > 1) rnorm(constants$n_env - 1, 0, 1) else numeric(0)),
	  u_gen = if (constants$n_gen > 0) rnorm(constants$n_gen, 0, 1) else numeric(0),
	  tau_gen = 1,
	  tau_gen_env = 1,
	  tau_env = rep(1, constants$n_env)
	)

	# Build and compile the model
	model <- nimbleModel(code = model_code, constants = constants, data = nimble_data, inits = inits)
	mcmc_config <- configureMCMC(model)
	mcmc_config$setMonitors(c("beta_0", "beta_env", "sigma_gen", "sigma_gen_env", "sigma_env"))
	mcmc <- buildMCMC(mcmc_config)
	compiled <- compileNimble(model, mcmc)

	# Run MCMC
	samples <- runMCMC(compiled$mcmc, niter = 10000, nburnin = 2000, thin = 2, nchains = 3, 
			samplesAsCodaMCMC = TRUE)

	
	exec_time1 <- proc.time() - ptm1
	list_computational_times[[ii1]] <- exec_time1
	ii1 <- ii1+1
	
	# save the results
	file_names <- paste0(dir_out, "r_out_04a_", package_name,"_test_", i1, "_2000_",  years_vec[i1] ,".RData")
	#save(samples, compiled, mcmc, mcmc_config, model, inits, nimble_data, constants, model_code, 
	# 			exec_time1, years_vec, seed_local,
# 				file = file_names)
	
	print(exec_time1)
	
	#rm(data_work, samples, compiled, mcmc, mcmc_config, model, inits, nimble_data, constants, model_code)


} # ends for (i1 in 2:s1)


# post processing the outputs
summary_out <- summary(samples)

# Post-process summary to match factor names in data_work
# Extract statistics and quantiles
stats_mat <- summary_out$statistics
quant_mat <- summary_out$quantiles
param_names <- rownames(stats_mat)

# Define mappings to factor levels in data_work
env_levels <- levels(data_work$env)
beta_map <- c("fixed_env_ref", paste("env", env_levels[-1], sep = ""))  # beta_env[1] is fixed at 0
sigma_map <- paste("sigma_", env_levels, sep = "")
scalar_map <- c("sigma_gen" = "sd_gen", "sigma_gen_env" = "sd_gen:env")

# Map indices to descriptive names
explicit_names <- param_names
beta_idx <- grep("^beta_env\\[\\d+\\]$", param_names)
if (length(beta_idx) > 0) {
  explicit_names[beta_idx] <- beta_map[seq_along(beta_idx)]
}
sigma_idx <- grep("^sigma_env\\[\\d+\\]$", param_names)
if (length(sigma_idx) > 0) {
  explicit_names[sigma_idx] <- sigma_map[seq_along(sigma_idx)]
}
for (orig in names(scalar_map)) {
  scalar_idx <- which(param_names == orig)
  if (length(scalar_idx) > 0) {
    explicit_names[scalar_idx] <- scalar_map[orig]
  }
}

# Reassign rownames
rownames(stats_mat) <- explicit_names
rownames(quant_mat) <- explicit_names

# Combine into a single data.table
tab_estimates <- data.table(
  rn = explicit_names,
  as.data.table(stats_mat),
  as.data.table(quant_mat)
)



# collect the estimates for asreml m2 reduced model, i1 = 2, to compare agains the rest of implementations
#write_fst(tab_estimates, path = "dir_mehtods_estimates_compare/fst_esimates_nimble_i1_2.fst")

