library(tidyverse)
library(ggpubr)

col = c("#40B0A6", "#E1BE6A")
col.grp = c("#FFC107", "#004D40", "#1E88E5")
sz  = 3

df.sim = read_csv(file.path("HGF_results", "main", "ppc_data_eHGF-L21.csv")) %>% 
  mutate(
    simulated = exp(yhat)
  )

load("../data/PAL-ADHD_data.RData")

df = merge(df.tsk, df.sim) %>%
  mutate(difficulty = if_else(difficulty == "difficult", "hard", difficulty),
         difficulty = factor(difficulty, levels = c("easy", "medium", "hard")),
         diagnosis  = if_else(diagnosis == "BOTH", "ADHD+ASD", diagnosis)) %>%
  filter(expected != "neutral") %>% droplevels() %>% 
  group_by(subID, diagnosis, phase, expected, difficulty, trl) %>%
  summarise(
    rt.cor    = mean(rt.cor),
    simulated = mean(exp(yhat))
  )

# expectancy effect as line plot
p1 = df %>% 
  group_by(subID, diagnosis, phase, expected) %>%
  summarise(rt.sim = median(simulated, na.rm = T),
            rt.cor = median(rt.cor, na.rm = T)) %>%
  group_by(subID, diagnosis, phase) %>%
  summarise(simulated = diff(rt.sim), 
            observed  = diff(rt.cor)) %>%
  pivot_longer(cols = c(simulated, observed), values_to = "rt", names_to = "type") %>%
  group_by(diagnosis, phase, type) %>%
  summarise(
    rt.mn = mean(rt),
    rt.se = sd(rt)
  ) %>%
  ggplot(aes(y = rt.mn, x = phase, 
             group = diagnosis, colour = diagnosis)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = 5) + 
  labs (x = "phase", y = "rt difference (ms)") + 
  ggtitle("Expectancy effect") + 
  facet_wrap(. ~ type, ncol = 2, scales = "free_y") + 
  scale_fill_manual(values = col.grp) +
  scale_color_manual(values = col.grp) +
  geom_hline(yintercept = 0) +
  theme_bw()  + 
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.78),
        axis.title.x = element_blank(), legend.title = element_blank(),
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5), text = element_text(size = 13))

# do the same for difficulty effects
p2 = df %>% filter(difficulty != "medium") %>%
  group_by(subID, diagnosis, phase, difficulty) %>%
  summarise(rt.sim = median(simulated, na.rm = T),
            rt.cor = median(rt.cor, na.rm = T)) %>%
  group_by(subID, diagnosis, phase) %>%
  summarise(simulated = diff(rt.sim), 
            observed  = diff(rt.cor)) %>%
  pivot_longer(cols = c(simulated, observed), values_to = "rt", names_to = "type") %>%
  group_by(diagnosis, phase, type) %>%
  summarise(
    rt.mn = mean(rt),
    rt.se = sd(rt)
  ) %>%
  ggplot(aes(y = rt.mn, x = phase, 
             group = diagnosis, colour = diagnosis)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = 5) + 
  labs (x = "phase", y = "rt difference (ms)") + 
  ggtitle("Difficulty effect") + 
  facet_wrap(. ~ type, ncol = 2, scales = "free_y") + 
  scale_fill_manual(values = col.grp) +
  scale_color_manual(values = col.grp) +
  theme_bw()  + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        text = element_text(size = 13), plot.title = element_text(hjust = 0.5))

p3 = df %>% rename("observed" = "rt.cor") %>% 
  pivot_longer(cols = c(observed, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, expected, type, phase) %>%
  summarise(value = sd(value, na.rm = T)) %>%
  group_by(type, diagnosis, phase) %>%
  summarise(
    rt.mn = mean(value),
    rt.se = sd(value)/sqrt(n())
  ) %>%
  ggplot(aes(y = rt.mn, x = phase, 
             group = diagnosis, colour = diagnosis)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = sz) + 
  labs(x = "", y = "SD rt (ms)") +
  scale_fill_manual(values = col.grp) +
  scale_color_manual(values = col.grp) +
  facet_wrap(. ~ type, ncol = 2, scales = "free_y") + 
  theme_bw() + 
  ggtitle("Reaction time variance") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        text = element_text(size = 13), plot.title = element_text(hjust = 0.5))

p = ggarrange(p1, p2, p3, nrow = 3, ncol = 1)
annotate_figure(p, top = text_grob("Reaction time and reaction time variance: observed and simulated data", 
                                   face = "bold", size = 14))
ggsave("plots/FigPPC_new.svg", plot = p, units = "cm", width = 27, height = 19)

# First version -----------------------------------------------------------

# expectancy
p1 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, expected) %>%
  summarise(value = median(value, na.rm = T)) %>%
  group_by(type, diagnosis, expected) %>%
  summarise(
    rt.mn = mean(value),
    rt.se = sd(value)/sqrt(n())
  ) %>%
  ggplot(aes(y = rt.mn, x = expected, 
             group = type, colour = type)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = sz) + 
  labs(x = "", y = "mean rt (ms)") +
  ylim(515, 685) + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  facet_wrap(. ~ diagnosis, ncol = 3) + 
  theme_bw() + 
  theme(legend.position.inside = c(0.2, 0.15), 
        legend.position = "inside", legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# difficulty
p2 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, difficulty) %>%
  summarise(value = median(value, na.rm = T)) %>%
  group_by(type, diagnosis, difficulty) %>%
  summarise(
    rt.mn = mean(value),
    rt.se = sd(value)/sqrt(n())
  ) %>%
  ggplot(aes(y = rt.mn, x = difficulty, 
             group = type, colour = type)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = sz) + 
  ylim(515, 685) + 
  labs(x = "", y = "mean rt (ms)") + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  facet_wrap(. ~ diagnosis, ncol = 3) + 
  theme_bw() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# expectancy
p4 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, phase) %>%
  summarise(value = sd(value, na.rm = T)) %>%
  group_by(type, diagnosis, phase) %>%
  summarise(
    rt.mn = mean(value),
    rt.se = sd(value)/sqrt(n())
  ) %>%
  ggplot(aes(y = rt.mn, x = phase, 
             group = type, colour = type)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  geom_errorbar(aes(ymin = rt.mn - rt.se, 
                    ymax = rt.mn + rt.se), linewidth = 1, width = 0.5, 
                position = position_dodge(0.4)) + 
  geom_point(position = position_dodge(0.4), size = sz) + 
  labs(x = "", y = "SD rt (ms)") +
  ylim(0, 160) + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  facet_wrap(. ~ diagnosis, ncol = 3) + 
  theme_bw() + 
  theme(legend.position = "none", legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


p = ggarrange(p1, p2, p4,
              ncol = 1, nrow = 3)
annotate_figure(p, top = text_grob("Posterior predictive checks: winning HGF model", 
                                   face = "bold", size = 14))
ggsave("plots/FigPPC.svg", units = "cm", width = 27, height = 18)

