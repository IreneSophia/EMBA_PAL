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

# setup caching of results
cache_dir = "./_brms_SBC_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

# number of simulations
nsim = 500

# set number of iterations and warmup for models
iter = 8000
warm = 2000

# set the maximum treedepth
maxdepth = 15

# get the data
dt.path = 'data/CV_pup_sum'
df = read_csv(file.path(dt.path, "CV_data_C(expected, Sum)[S.expected].csv"),
              show_col_types = F) %>%
  drop_na() %>%
  merge(., read_csv(file.path('data', 'PAL_data.csv'), show_col_types = F) %>%
          select(subID, diagnosis) %>%
          distinct()) %>%
  mutate_if(is.character, as.factor)

# print the contrasts
contrasts(df$diagnosis) = contr.sum(4)
contrasts(df$diagnosis)
contrasts(df$expected)  = contr.sum(2)
contrasts(df$expected)

# code for filenames
code = "PAL-pup"

# model formula
f.pup = brms::bf(rel_pupil ~ diagnosis * expected + rts +
                   (1 | subID)
)

# set weakly informative priors
priors = c(
  prior(normal(0,    0.05),  class = Intercept),
  prior(normal(0.15, 0.10),  class = sigma),
  prior(normal(0.10, 0.10),  class = sd),
  prior(normal(0,    0.03),  class = b)
)

# change Intercept based on mean in the data > only thing adjusted between models
priors = priors %>%
  mutate(
    prior = if_else(
      class == "Intercept", 
      gsub("\\(.*,", paste0("(", round(mean(df$rel_pupil),3), ", "), prior), prior)
  )

# load the data and prepare SBC
gen = SBC_generator_brms(f.pup, data = df, prior = priors, 
                         thin = 50, warmup = 10000, refresh = 2000,  
                         control = list(max_treedepth = maxdepth),
                         generate_lp = TRUE, init = 0.1)

if (!file.exists(file.path(cache_dir, paste0("dat_", code, ".rds")))){
  dat = generate_datasets(gen, nsim)
  saveRDS(dat, file = file.path(cache_dir, paste0("dat_", code, ".rds")))
} else {
  dat = readRDS(file = file.path(cache_dir, paste0("dat_", code, ".rds")))
}

bck = SBC_backend_brms_from_generator(gen, chains = 4, thin = 1,
                                      warmup = warm, iter = iter)

# get the last number
ls.files = list.files(path = cache_dir, pattern = sprintf("res_%s_([0-9]+).rds", code))
if (is_empty(ls.files)) {
  i = 1
} else {
  i = max(as.numeric(gsub(sprintf("res_%s_(.+).rds", code), "\\1", ls.files))) + 1
}
set.seed(2468+i)
m = 25

# perform the SBC 
write(sprintf("%s: Starting %i simulations (total: %d).", Sys.time(), m, warm+iter), file.path("logfiles", "log_PAL.txt"), append = TRUE)
dat_part = SBC_datasets(dat$variables[((i-1)*m + 1):(i*m),],
                        dat$generated[((i-1)*m + 1):(i*m)])
res = compute_SBC(dat_part, bck,
                  cache_mode = "results",
                  cache_location = file.path(cache_dir, sprintf("res_%s_%02d", code, i)))
write(sprintf("%s: DONE.", Sys.time()), file.path("logfiles/", "log_PAL.txt"), append = TRUE)

