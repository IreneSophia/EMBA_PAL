library(tidyverse)
library(ggpubr)

col = c("#40B0A6", "#E1BE6A")
sz  = 3

df.sim = read_csv(file.path("HGF_results", "main", "ppc_data_HGF-L17.csv")) %>% 
  group_by(subID, trl) %>%
  summarise(
    simulated = mean(exp(yhat))
  )

load("../data/PAL-ASD_data.RData")

df = merge(df.tsk, df.sim) %>%
  mutate(difficulty = if_else(difficulty == "difficult", "hard", difficulty),
         difficulty = factor(difficulty, levels = c("easy", "medium", "hard"))) %>%
  filter(expected != "neutral") %>% droplevels()

# expectancy
p1 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, expected) %>%
  summarise(value = median(value, na.rm = T)) %>%
  group_by(diagnosis, type, expected) %>%
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
  labs (x = "", y = "rt (ms)") +
  ylim(515, 665) + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  facet_wrap(. ~ diagnosis, nrow = 2) + 
  theme_bw() + 
  theme(legend.position.inside = c(0.35, 0.6), 
        legend.position = "inside", legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# difficulty
p2 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, difficulty) %>%
  summarise(value = median(value, na.rm = T)) %>%
  group_by(diagnosis, type, difficulty) %>%
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
  ylim(515, 665) + 
  labs (x = "", y = "") + 
  facet_wrap(. ~ diagnosis, nrow = 2) + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_bw() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


# phase
p3 = df %>% rename("data" = "rt.cor") %>% 
  pivot_longer(cols = c(data, simulated), names_to = "type") %>%
  group_by(subID, diagnosis, type, phase) %>%
  summarise(value = median(value, na.rm = T)) %>%
  group_by(diagnosis, type, phase) %>%
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
  labs (x = "", y = "") + 
  ylim(515, 665) + 
  facet_wrap(. ~ diagnosis, nrow = 2) + 
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_bw() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


p = ggarrange(p1, p2, p3, 
              nrow = 1, ncol = 3)
annotate_figure(p, top = text_grob("Posterior predictive checks: HGF", 
                                   face = "bold", size = 14))
ggsave("plots/FigPPC.svg", units = "cm", width = 24, height = 12)

