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
setwd("/home/emba/Documents/EMBA/EMBA_PAL_scripts")

# setup caching of results
cache_dir = "./_brms_SBC_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

# load the trial order
df.trl = read_csv('data/PAL_scheme.csv', show_col_types = F) %>%
  select(trl, emo, tone, difficulty, phase, expected)

# create the data
df.pal = data.frame()
groups = c("A", "B")
nos    = 22 # number of subjects per group
for (g in groups) {
  for (i in 1:nos) {
    df.pal = rbind(df.pal, 
                   df.trl %>% mutate(subID = sprintf("%s%02d", g, i),
                                     diagnosis = g))
  }
}

# add bogus rts
df.pal = df.pal %>% mutate(rt.cor = 1)

# drop the neutral condition for the analysis
df.pal = df.pal %>%
  filter(expected != "neutral") %>% droplevels() %>%
  mutate_if(is.character, as.factor) %>%
  ungroup()

# set and print the contrasts
contrasts(df.pal$diagnosis) = contr.sum(2)
contrasts(df.pal$diagnosis)
contrasts(df.pal$expected) = contr.sum(2)
contrasts(df.pal$expected)
contrasts(df.pal$phase) = contr.sum(3)
contrasts(df.pal$phase)
contrasts(df.pal$difficulty) = contr.sum(3)
contrasts(df.pal$difficulty)

nsim = 250

code = "PAL-rt"

# set the formula
f.pal = brms::bf(rt.cor ~ diagnosis * expected * phase * difficulty +
                   (expected * phase * difficulty | subID) + 
                   (diagnosis | trl))

# informative priors based Lawson et al. and Schad, Betancourt & Vasishth (2019)
priors = c(
  prior(normal(6.0,   0.3),   class = Intercept),
  prior(normal(0.0,   0.5),   class = sigma),
  prior(normal(0,     0.1),   class = sd),
  prior(lkj(2),              class = cor),
  prior(normal(100, 100.0),   class = ndt),
  prior(normal(0.00,  0.04),  class = b)
)

# load the data and prepare SBC
gen = SBC_generator_brms(f.pal, data = df.pal, prior = priors, 
                         thin = 50, warmup = 10000, refresh = 2000,
                         generate_lp = TRUE, family = shifted_lognormal, init = 0.1)
if (!file.exists(file.path(cache_dir, paste0("dat_", code, ".rds")))){
  dat = generate_datasets(gen, nsim)
  saveRDS(dat, file = file.path(cache_dir, paste0("dat_", code, ".rds")))
} else {
  dat = readRDS(file = file.path(cache_dir, paste0("dat_", code, ".rds")))
}
warm = 1500
iter = 6000

bck = SBC_backend_brms_from_generator(gen, chains = 4, thin = 1,
                                      warmup = warm, iter = iter,
                                      init = 0.1)

# get the last number
ls.files = list.files(path = cache_dir, pattern = sprintf("res_%s_([0-9]+).rds", code))
if (is_empty(ls.files)) {
  i = 1
} else {
  i = max(as.numeric(gsub(sprintf("res_%s_(.+).rds", code), "\\1", ls.files))) + 1
}
set.seed(2468+i)
m = 50

# perform the SBC 
write(sprintf("%d: Starting %i simulations (total: %d).", i, m, warm+iter), 
      file.path("/home/emba/pCloudDrive/LOGFILES", "log_PAL.txt"), append = TRUE)
dat_part = SBC_datasets(dat$variables[((i-1)*m + 1):(i*m),],
                        dat$generated[((i-1)*m + 1):(i*m)])
res = compute_SBC(dat_part, bck, keep_fits = F,
                  cache_mode = "results",
                  cache_location = file.path(cache_dir, sprintf("res_%s_%02d", code, i)))
write("DONE.", 
      file.path("/home/emba/pCloudDrive/LOGFILES", "log_PAL.txt"), append = TRUE)

