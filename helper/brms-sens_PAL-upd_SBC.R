library(ggplot2)         # plots
library(brms)            # Bayesian lmms
library(tidyverse)       # tibble stuff
library(SBC)
library(bayesplot)

# set options for SBC
use_cmdstanr = getOption("SBC.vignettes_cmdstanr", TRUE) # Set to false to use rstan instead
options(brms.backend = "cmdstanr")

# using parallel processing
library(future)
plan(multisession)

# set working directory
setwd("..")

pr.descriptions = c("chosen",
                    "sdx2",    "sdx4",   "sdx8",
                    "sdx0.5", "sdx0.25", "sdx0.125"
)


# Setup caching of results
cache_dir = "./_brms_SBC_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
brms_dir = "./_brms_models"
if(!dir.exists(brms_dir)) {
  dir.create(brms_dir)
}

# load the function to perform the sensitivity analysis
source('helper/fun_bf-sens.R')

# set the directory in which to save results
sense_dir = file.path(getwd(), "_brms_sens_cache")
main.code = "hgf_alpha"

# rerun the model with more iterations for bridgesampling
set.seed(5544)
m.orig = brm(f.alpha, family = lognormal,
             df.upd, prior = priors,
             iter = 40000, warmup = 10000,
             backend = "cmdstanr", threads = threading(8),
             file = file.path(brms_dir, "m_hgf_alpha_bf"), silent = 2,
             save_pars = save_pars(all = TRUE)
)

# loop through the priors that have not been used before
for (pr.desc in pr.descriptions) {
  tryCatch({
    # use function to compute BF with the described priors
    bf_sens_3int(m.orig, "diagnosis", "level", "change",
                 pr.desc,
                 main.code, # prefix for all models and MLL
                 file.path("logfiles", "log_PAL_bf.txt"), # log file
                 sense_dir # where to save the models and MLL
    )
  },
  error = function(err) {
    message(sprintf("Error for %s: %s", pr.desc, err))
  }
  )
}