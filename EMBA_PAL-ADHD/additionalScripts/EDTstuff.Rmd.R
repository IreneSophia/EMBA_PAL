# Association between EDT and HGF parameters

## Model setup and posterior predictive checks

```{r edt, fig.height=6}

# set some priors
priors = c(
  prior(normal(0, 0.50), class = Intercept),
  prior(normal(0, 1.00), class = b),
  prior(normal(0, 0.50), class = sigma)
)

# fit the model
m.edt = brm(EDT ~ sbe1 + sbe2 + sbe3 + sze + som2 + som3 , 
            seed = 8282,
            df.hgf, prior = priors,
            iter = iter, warmup = warm,
            backend = "cmdstanr", threads = threading(t),
            file = file.path(brms_dir, "m_hgf_edt"),
            save_pars = save_pars(all = TRUE)
)
rstan::check_hmc_diagnostics(m.edt$fit)

# check that rhats are below 1.01
sum(brms::rhat(m.edt) >= 1.01, na.rm = T)

# check the trace plots
post.draws = as_draws_df(m.edt)
mcmc_trace(post.draws, regex_pars = "^b_",
           facet_args = list(ncol = 3)) +
  scale_x_continuous(breaks=scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3))

```
```{r postpc_edt, fig.height=2}

# check the fit of the predicted data compared to the real data
p = pp_check(m.edt, ndraws = nsim) + 
  theme_bw() + theme(legend.position = "none") + xlim(0, 1)

annotate_figure(p, top = text_grob("Posterior predictive checks", 
                                   face = "bold", size = 14))

```

This model fits the data well. 

## Inferences

Now that we are convinced that we can trust our model, we have a look at its estimate and use the hypothesis function to assess our hypotheses and perform explorative tests. 

```{r inf_edt, fig.height=6}

# print a summary
summary(m.edt)

# get the estimates and compute group comparisons
df.m.edt = post.draws %>% 
  select(starts_with("b_"))

# plot the posterior distributions
df.m.edt %>%
  select(starts_with("b_")) %>%
  pivot_longer(cols = starts_with("b_"), names_to = "coef", values_to = "estimate") %>%
  subset(!startsWith(coef, "b_Int")) %>%
  mutate(
    coef = substr(coef, 3, nchar(coef)),
    coef = fct_reorder(coef, desc(estimate))
  ) %>% 
  group_by(coef) %>%
  mutate(
    cred = case_when(
      (mean(estimate) < 0 & quantile(estimate, probs = 0.95) < 0) |
        (mean(estimate) > 0 & quantile(estimate, probs = 0.05) > 0) ~ "strong",
      T ~ "x"
    )
  ) %>% ungroup() %>%
  ggplot(aes(x = estimate, y = coef, fill = cred)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ggdist::stat_halfeye(alpha = 0.7) + ylab(NULL) + theme_bw() +
  scale_fill_manual(values = c(c_dark, c_light)) + theme(legend.position = "none")

# assess the parameters
e.om2 = hypothesis(m.edt, "0 > som2")
e.om2
e.ze  = hypothesis(m.edt, "0 < sze")
e.ze

# calculate effect sizes
df.effect = post.draws %>%
  mutate(across(starts_with("sd")|starts_with("sigma"), ~.^2)) %>%
  mutate(
    sumvar = sqrt(rowSums(select(., starts_with("sd")|starts_with("sigma"))))
  ) %>% select(-b_Intercept) %>% pivot_longer(cols = starts_with("b_")) %>%
  mutate(
    value = value/sumvar
  ) %>% pivot_wider(names_prefix = "e")

kable(df.effect %>% select(starts_with("e")) %>%
        pivot_longer(cols = everything(), values_to = "estimate") %>%
        group_by(name) %>%
        summarise(
          ci.lo = lower_ci(estimate),
          mean  = mean(estimate),
          ci.hi = upper_ci(estimate),
          interpret = interpret_cohens_d(mean)
        ), digits = 3
)


```

# Plots for HGF parameters

```{r plot_hgf, fig.height=6}

p = df.hgf %>%
  mutate(diagnosis = if_else(diagnosis == "BOTH", "ADHD+ASD", diagnosis)) %>%
  select(subID, diagnosis, be1, be2, be3, ze, om2, om3) %>% #
  pivot_longer(cols = c(be1, be2, be3, ze, om2, om3), 
               names_to = "parameter") %>%
  mutate(
    parameter = factor(case_match(parameter,
                                  "be1" ~ "stimulus surprise",
                                  "be2" ~ "precision-weighted PE",
                                  "be3" ~ "phasic volatility",
                                  "ze"  ~ "Sigma (decision noise)",
                                  "om2" ~ "2nd tonic volatility",
                                  "om3" ~ "3rd tonic volatility"
    ), levels = c("2nd tonic volatility", 
                  "3rd tonic volatility", 
                  "stimulus surprise", 
                  "precision-weighted PE", 
                  "phasic volatility", 
                  "Sigma (decision noise)"))
  ) %>%
  ggplot(aes(x = 1, y = value, fill = diagnosis, colour = diagnosis)) + #
  geom_rain(rain.side = 'r',
            boxplot.args = list(color = "black", outlier.shape = NA, show.legend = FALSE, alpha = .8),
            violin.args = list(color = "black", outlier.shape = NA, alpha = .8),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = 0, width = 0.3), width = 0.3
            ),
            point.args = list(show_guide = FALSE, alpha = .5),
            violin.args.pos = list(
              width = 0.6, position = position_nudge(x = 0.16)),
            point.args.pos = list(position = ggpp::position_dodgenudge(x = -0.25, width = 0.1))) +
  scale_fill_manual(values = col.grp) +
  scale_color_manual(values = col.grp) +
  facet_wrap(. ~ parameter, scales = "free", ncol = 3) +
  theme_bw() + 
  theme(legend.position = "bottom", plot.title = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        text = element_text(size = 13), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,0,0,0))

# Alternative plotting
plist = list()
cols  = c("om2", "om3", "be1", "be2", "be3")
names = c("tonic volatility (cue-outcome)", "tonic volatility (environment)", 
          "stimulus surprise", "precision-weighted PE", "phasic volatility")
for (i in 1:length(cols)) {
  df.hgf$value     = df.hgf[, cols[i]]
  df.hgf$parameter = names[i]
  if (cols[i] == "be2") {
    plist[[i]] = ggplot(df.hgf, aes(x = 1, y = value, 
                                    fill = diagnosis, colour = diagnosis)) + 
      geom_rain(rain.side = 'r',
                boxplot.args = list(color = "black", outlier.shape = NA, show.legend = FALSE, alpha = .8),
                violin.args = list(color = "black", outlier.shape = NA, alpha = .8),
                boxplot.args.pos = list(
                  position = ggpp::position_dodgenudge(x = 0, width = 0.3), width = 0.3
                ),
                point.args = list(show_guide = FALSE, alpha = .5),
                violin.args.pos = list(
                  width = 0.6, position = position_nudge(x = 0.16)),
                point.args.pos = list(position = ggpp::position_dodgenudge(x = -0.25, width = 0.1))) +
      scale_fill_manual(values = col.grp) +
      scale_color_manual(values = col.grp) +
      facet_wrap(. ~ parameter) +
      ylim(-0.25, 0.005) + 
      theme_bw() + 
      theme(legend.position = "inside", legend.position.inside = c(0.5, 0.1),
            legend.direction = "horizontal", legend.title = element_blank(), 
            plot.title = element_blank(), 
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            text = element_text(size = 13), axis.text.x=element_blank(), 
            axis.ticks.x=element_blank())
    
  } else {
    plist[[i]] = ggplot(df.hgf, aes(x = 1, y = value, 
                                    fill = diagnosis, colour = diagnosis)) + 
      geom_rain(rain.side = 'r',
                boxplot.args = list(color = "black", outlier.shape = NA, show.legend = FALSE, alpha = .8),
                violin.args = list(color = "black", outlier.shape = NA, alpha = .8),
                boxplot.args.pos = list(
                  position = ggpp::position_dodgenudge(x = 0, width = 0.3), width = 0.3
                ),
                point.args = list(show_guide = FALSE, alpha = .5),
                violin.args.pos = list(
                  width = 0.6, position = position_nudge(x = 0.16)),
                point.args.pos = list(position = ggpp::position_dodgenudge(x = -0.25, width = 0.1))) +
      scale_fill_manual(values = col.grp) +
      scale_color_manual(values = col.grp) +
      facet_wrap(. ~ parameter) +
      theme_bw() + 
      theme(legend.position = "none", plot.title = element_blank(), 
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            text = element_text(size = 13), axis.text.x=element_blank(), 
            axis.ticks.x=element_blank())
  }
}

# add the correlation
df.hgf$parameter = "association with emotion recognition"
plist[[i+1]] = df.hgf %>%
  ggplot(., aes(x = EDT, y = om2)) + 
  geom_point(aes(fill = diagnosis, colour = diagnosis), alpha = 0.8) + 
  geom_smooth(method = "lm", colour = "#D81B60", fill = "#D81B60", alpha = 0.25) +
  scale_fill_manual(values = col.grp) +
  scale_color_manual(values = col.grp) +
  xlab("emotion discrimination threshold") + ylab("tonic volatility (cue-outcome)") +
  facet_wrap(. ~ parameter) +
  theme_bw() + 
  theme(legend.position = "none", plot.title = element_blank(), 
        text = element_text(size = 13), axis.title = element_text(size = 11))

p = ggarrange(plotlist = plist, ncol = 3, nrow = 2)
p.a = annotate_figure(p, top = text_grob("Participant-specific HGF parameters", 
                                         face = "bold", size = 14))

ggsave("plots/FigHGF.svg", plot = p.a, units = "cm", width = 27, height = 16)

# include medication
df.hgf %>%
  mutate(diagnosis = if_else(diagnosis == "BOTH", "ADHD+ASD", diagnosis),
         adhd.meds.bin = case_when(adhd.meds.bin == "TRUE" ~ "medicated", 
                                   T ~ ""),
         group = paste0(diagnosis, adhd.meds.bin)) %>%
  select(subID, diagnosis, group, be1, be2, be3, ze, om2, om3) %>% #
  pivot_longer(cols = c(be1, be2, be3, ze, om2, om3), 
               names_to = "parameter") %>%
  mutate(
    parameter = factor(case_match(parameter,
                                  "be1" ~ "stimulus surprise",
                                  "be2" ~ "precision-weighted PE",
                                  "be3" ~ "phasic volatility",
                                  "ze"  ~ "Sigma (decision noise)",
                                  "om2" ~ "2nd tonic volatility",
                                  "om3" ~ "3rd tonic volatility"
    ), levels = c("2nd tonic volatility", 
                  "3rd tonic volatility", 
                  "stimulus surprise", 
                  "precision-weighted PE", 
                  "phasic volatility", 
                  "Sigma (decision noise)"))
  ) %>%
  ggplot(aes(x = diagnosis, y = value, fill = group, colour = group)) + #
  geom_rain(rain.side = 'r',
            boxplot.args = list(color = "black", outlier.shape = NA, show.legend = FALSE, alpha = .8),
            violin.args = list(color = "black", outlier.shape = NA, alpha = .8),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = 0, width = 0.3), width = 0.3
            ),
            point.args = list(show_guide = FALSE, alpha = .5),
            violin.args.pos = list(
              width = 0.6, position = position_nudge(x = 0.16)),
            point.args.pos = list(position = ggpp::position_dodgenudge(x = -0.25, width = 0.1))) +
  #scale_fill_manual(values = col.grp) +
  #scale_color_manual(values = col.grp) +
  facet_wrap(. ~ parameter, scales = "free", ncol = 3) +
  theme_bw() + 
  theme(legend.position = "bottom", plot.title = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        text = element_text(size = 13), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,0,0,0))


```

