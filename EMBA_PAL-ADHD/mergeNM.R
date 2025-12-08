library(tidyverse)

# get HGF parameters  
df = read_csv(file.path("HGF_results/main", "eHGF-L21_traj.csv")) %>%
  select(subID, diagnosis, trl, alpha2, alpha3) %>% ungroup() %>%
  mutate(
    # code the phases > only take the beginning and end of volatile
    phase = case_when(
      trl < 73  ~ "pre",
      trl > 264 ~ "post",
      trl < 145 ~ "vol1",
      trl > 192 ~ "vol2"
    )
  ) %>%
  drop_na() %>%
  group_by(subID, diagnosis, phase) %>%
  summarise(
    alpha2 = median(alpha2),
    alpha3 = median(alpha3)
  ) %>%
  pivot_wider(names_from = phase, id_cols = c(subID, diagnosis), values_from = starts_with("alpha")) %>%
  group_by(subID, diagnosis) %>%
  summarise(
    alpha2_pre2vol  = abs(alpha2_pre  - alpha2_vol1),
    alpha2_vol2post = abs(alpha2_post - alpha2_vol2),
    alpha3_pre2vol  = abs(alpha3_pre  - alpha3_vol1),
    alpha3_vol2post = abs(alpha3_post - alpha3_vol2)
  ) %>%
  merge(., read_csv(file.path("HGF_results/main", "eHGF-L21_results.csv"))) %>%
  select(subID, diagnosis, starts_with("alpha"), be1, be2, be3,
         ze, om2, om3) %>%
  mutate(
    ADHD = if_else(diagnosis == "COMP", 0, 1)
  )

load(file.path("./_ddm_models", "fit_hddm.RData"))

df = merge(df, df.part %>%
             rename("subID" = "s")) %>%
  relocate(subID, diagnosis, ADHD) %>%
  arrange(ADHD)

write_csv(df, "NM_data.csv")
