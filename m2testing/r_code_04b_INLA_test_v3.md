This is the R INLA implementation of the reduced asreml m2 model. It appears it can only run to certan sample sizes without producing an error, smaller than greta.

For i1 = 2, the fixed effects `m1_inla$summary.fixed` appear to be close to asreml and greta (just visual examination).

Random effects `m1_inla$summary.random` need to be further organised since they do not contain variable names.

```
# source("r_code_04b_INLA_test_v3.R")

##############
# Author: Max Moldovan
# max.moldovan@adelaide.edu.au
# version 0.0.1, 12 Sep 2025
# version 0.0.3, 09 Oct 2025 reduced asreml m2 model, see the same date emails
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
require(INLA)
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
dataset_expend <- TRUE

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


package_name <- "INLAv2"

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

times_collect_vec <- c()
n_collect_vec <- c()
i55 <- 1


s1 <- length(years_vec)
set.seed(55)
rand_seeds_vec <- sample(10:1000, s1, replace = F)

years_to_check <- c(5, 10, 30, 100, 150, 200, 275)

# run through 2000:2001, 2000:2002, ..., 2000:2010
for (i1 in 2) {
#for (i1 in years_to_check) {

	years_needed <- years_vec[1:i1]
	data_work <- copy(data_initial)
	seed_local <- rand_seeds_vec[i1]
	
	print(paste0("number of years to be processed: ", length(years_needed)))
	
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

	# Create indices for random effects
	data_work[, gen_id := as.numeric(gen)]
	data_work[, rep_env_id := as.numeric(interaction(rep, env))]
	data_work[, gen_yearf_id := as.numeric(interaction(gen, yearf))]
	data_work[, gen_loc_id := as.numeric(interaction(gen, loc))]
	data_work[, gen_env_id := as.numeric(interaction(gen, env))]

	# INLA formula
	formula_inla <- yield ~ env +
	  f(gen_id, model = "iid") +
	  #f(rep_env_id, model = "iid") +
	  #f(gen_yearf_id, model = "iid") +
	  #f(gen_loc_id, model = "iid") +
	  f(gen_env_id, model = "iid")

	# Note: INLA doesn't directly support environment-specific residual variances
	# We can approximate by using replicate structure
	set.seed(seed_local)
	m1_inla <- inla(
	  formula_inla,
	  data = data_work,
	  control.compute = list(dic = TRUE, waic = TRUE),
	  control.predictor = list(compute = TRUE)
	)
	
	
	exec_time1 <- proc.time() - ptm1
	list_computational_times[[ii1]] <- exec_time1
	ii1 <- ii1+1
	
	times_collect_vec[i55] <- exec_time1[[3]]
	n_collect_vec[i55] <- nrow(data_work)
	i55 <- i55 + 1
	
	# save the results
	file_names <- paste0(dir_out, "r_out_04a_", package_name,"_test_", i1, "_2000_",  years_vec[i1] ,".RData")
	#save(m1_inla, exec_time1, years_vec, seed_local, file = file_names)
	
	#rm(m1_inla, formula_inla, data_work)


} # ends for (i1 in 2:s1)

#tab_times_summary <- data.table(sample_size = n_collect_vec, times_sec = times_collect_vec)
#write_fst(tab_times_summary, path = "fst_tab_INLA_times_m2.fst")

# plot(tab_times_summary[, sample_size], tab_times_summary[, times_sec],
#      	xlab = "Sample sizes",  # X-axis label (customize text as needed)
#      	ylab = "exacution times, sec",   # Y-axis label (optional, for symmetry)
#      	main = "INLA")  # Main title (optional)

```
