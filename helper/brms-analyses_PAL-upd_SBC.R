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

nsim = 500

# get data 
df.hgf = read_csv(file.path("HGF_results/main", "eHGF-C_results.csv")) %>%
  mutate_if(is.character, as.factor)

# get alpha updates in long format
df.upd = df.hgf %>%
  select(subID, diagnosis, 
         alpha2_pre2vol, alpha2_vol2post, 
         alpha3_pre2vol, alpha3_vol2post) %>%
  pivot_longer(
    cols = c(alpha2_pre2vol, alpha2_vol2post, 
             alpha3_pre2vol, alpha3_vol2post)) %>%
  mutate(
    # take the absolute value because we are interested in the size of the change
    value = abs(value)
  ) %>% 
  separate(name, into = c("level", "change")) %>%
  mutate_if(is.character, as.factor)

# set and print the contrasts
contrasts(df.upd$diagnosis) = contr.sum(4)
contrasts(df.upd$diagnosis)
contrasts(df.upd$change) = contr.sum(2)
contrasts(df.upd$change)
contrasts(df.upd$level) = contr.sum(2)
contrasts(df.upd$level)

# code for filenames
code = "PAL_alpha"

# model formula
f.alpha = brms::bf( value ~ diagnosis * level * change + (level + change | subID) )

# set weakly informative priors taking Lawson 2017 into consideration
priors = c(
  prior(normal(-5, 2),    class = Intercept),
  prior(normal(0.5, 0.5), class = sigma),
  prior(normal(0.5, 0.5), class = sd),
  prior(lkj(2),           class = cor),
  prior(normal(0, 0.5),   class = b)
)

# load the data and prepare SBC
gen = SBC_generator_brms(f.alpha, data = df.upd, prior = priors, family = lognormal,
                         thin = 50, warmup = 10000, refresh = 2000,
                         generate_lp = TRUE, init = 0.1)

if (!file.exists(file.path(cache_dir, paste0("dat_", code, ".rds")))){
  dat = generate_datasets(gen, nsim)
  saveRDS(dat, file = file.path(cache_dir, paste0("dat_", code, ".rds")))
} else {
  dat = readRDS(file = file.path(cache_dir, paste0("dat_", code, ".rds")))
}
warm = 1000
iter = 3000

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
m = 100

# perform the SBC 
write(sprintf("%d: Starting %i simulations (total: %d).", i, m, warm+iter), file.path("logfiles", "log_PAL-HGF.txt"), append = TRUE)
dat_part = SBC_datasets(dat$variables[((i-1)*m + 1):(i*m),],
                        dat$generated[((i-1)*m + 1):(i*m)])
res = compute_SBC(dat_part, bck,
                  cache_mode = "results",
                  cache_location = file.path(cache_dir, sprintf("res_%s_%02d", code, i)))
write("DONE.", file.path("logfiles/", "log_PAL.txt"), append = TRUE)

