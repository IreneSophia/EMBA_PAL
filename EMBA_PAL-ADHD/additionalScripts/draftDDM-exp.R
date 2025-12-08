library(ggdmc)
library(tidyverse)

# options for summarise
options(dplyr.summarise.inform = FALSE)

# set cores
options(mc.cores = parallel::detectCores())

# set options for brms
options(brms.backend = "cmdstanr")

setwd("/Users/vilya/Downloads/ddm-scripts")

# Setup caching of results
ddm_dir = "./_ddm_models-exp"
if(!dir.exists(ddm_dir)) {
  dir.create(ddm_dir)
}



# build the model
model = BuildModel(
  p.map     = list(a = "1", v = "E", z = "1", d = "1", sz = "1", sv = "1",
                   t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2"), E = c("exp", "un")),
  constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
  responses = c("r1", "r2"),
  type      = "rd")
npar = length(GetPNames(model))
nop = 2   # number of levels

# determine the priors
pop.mean  = c(a = 2.0, v.exp = 2.5, v.un = 2.0, z = 0.5, t0 = 0.3)
pop.scale = c(a = 0.5, v.exp = 0.5, v.un = 0.5, z = 0.1, t0 = 0.05)
pop.lower = c(0, rep(-5, nop), rep(0, npar-1-nop))
pop.upper = c(5, rep( 7, nop), rep(1, npar-1-nop))

p.prior = BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale*5,
  lower = pop.lower,
  upper = pop.upper)
mu.prior = p.prior
plot(p.prior)

sigma.prior = BuildPrior(
  dists = rep("beta", npar),
  p1    = rep(1, npar),
  p2    = rep(1, npar),
  upper = rep(2, npar))

names(sigma.prior) = GetPNames(model)

# collect the priors for the hierarchical model
priors = list(pprior = p.prior, location = mu.prior, scale = sigma.prior)


# load in the data
df = read_csv("../data/PAL-ADHD_data.csv", show_col_types = F) %>%
  # filter out super long and super fast reactions
  filter(rt > 100 & rt < 1500) %>%
  # add all the trial information
  merge(., read_csv("../data/PAL_scheme.csv"), show_col_types = F) %>%
  mutate(
    R = if_else((acc & emo == "positive") | (!acc & emo == "negative"), 
                "r1", "r2"),
    P = case_match(phase, "prevolatile" ~ "pre", "volatile" ~ "vol", "postvolatile" ~ "post"),
    E = case_match(expected, "expected" ~ "exp", "unexpected" ~ "un"),
    S = if_else(emo == "positive", "s1", "s2"),
    s = as.factor(subID),
    RT = rt/1000
  ) %>% select(diagnosis, s, S, E, R, RT) %>% drop_na() %>%
  mutate_if(is.character, as.factor)

m       = list()
df.part = data.frame()

groups = unique(df$diagnosis)
groups = "COMP"

for (group in groups) {
  
  if (!file.exists(file.path(ddm_dir, sprintf("m_%s.rds", group)))) {
    # select part of the data
    df.grp = df %>% filter(diagnosis == group) %>% 
      select(-diagnosis) %>% droplevels()
    
    # create the data model instance 
    dmi = BuildDMI(df.grp, model)
    
    # do this until rhats are okay
    rhat  = 5
    count = 0
    
    while (rhat >= 1.1) {
      
      # increase counter and thinning
      count = count + 1
      print(sprintf("%d: %s", count, group))
      tictoc::tic()
      
      # start sampling
      fit0 = StartNewsamples(data = dmi, prior = priors)
      
      # run more to get it to be stable
      fit = run(fit0, thin = 32, nmc = 1000) # ~33min
      
      # check the rhats
      rhats = hgelman(fit, verbose = TRUE)
      rhat  = max(rhats)# if rhats are fine, then save the fit
      
      tictoc::toc()
      
    }
    
  } else {
    # if it already exists, just load the fit
    fit = readRDS(file.path(ddm_dir, sprintf("m_%s.rds", group)))
  }
  
  # if rhats are fine, then save the "ultimate" fit
  saveRDS(fit, file.path(ddm_dir, sprintf("m_%s.rds", group)))
  
}