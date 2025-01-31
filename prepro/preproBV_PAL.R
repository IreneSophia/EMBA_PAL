library(tidyverse)       # tibble stuff

dt.path = paste('/home/emba/Documents/EMBA', 'BVET', sep = "/")
#dt.path = paste('/home/emba/Documents/EMBA', 'BVET-explo', sep = "/")

# load the relevant data in long format
df.tsk = list.files(path = dt.path, pattern = "PAL-BV*", full.names = T) %>%
  map_df(~read_csv(., show_col_types = F, col_types = "ccddddccddcd")) %>% 
  filter(is.na(key)|(key != "tstart" & key != "tend")) %>%
  mutate(
    # sometimes participants used other rows instead of the middle one
    # or moved right with their ring finger to plus instead of 6
    key = case_when(
        key == "+" | key == "9" | key == "3" | key == "6" ~ "right",
        key == "7" | key == "1" | key == "4" ~ "left"
      ),
    key_pos = case_when(
      key_pos == "rught" ~ "right",
      T ~ key_pos
    )
    ) %>%
  select(subID, key_pos, trl, key, rt)

# load the information on the trials from the experiment file and merge with data
df.tsk = merge(df.tsk, 
               read_csv(file.path(dt.path, 'PAL_scheme-pic.csv'), show_col_types = F)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(
    choice = case_when(
           key == "right" & key_pos == "right" ~ "positive",
           key == "left"  & key_pos == "right" ~ "negative",
           key == "right" & key_pos == "left"  ~ "negative",
           key == "left"  & key_pos == "left"  ~ "positive",
           ),
    acc = case_when(emo == choice ~ T,
                    T ~ F),
    # change factor levels for phase and difficulty
    phase      = fct_relevel(phase, c("prevolatile", "volatile", "postvolatile")),
    difficulty = fct_relevel(difficulty, c("easy", "medium", "difficult"))
  ) %>%
  group_by(subID) %>%
  mutate(
         # replace very short and very long RTs with NAs
         rt.upper = quantile(rt, 0.75, na.rm = T) + 1.5 * IQR(rt, na.rm = T),
         rt.lower = quantile(rt, 0.25, na.rm = T) - 1.5 * IQR(rt, na.rm = T),
         use      = if_else(rt > rt.lower & rt < rt.upper, T, F),
         rt.cor   = if_else(use == T & acc == T, rt, NA),
         rt.use   = if_else(use == T & !is.na(key), rt, NA)
         ) %>%
  select(subID, phase, trl, expected, emo, tone, ut, difficulty, 
         rt.cor, rt.use, acc) 

# does anyone have to be excluded?
exc = df.tsk %>% group_by(subID) %>% summarise(acc = mean(!is.na(rt.cor))) %>% filter(acc < 2/3)
exc = as.character(exc$subID)

# load pilot participants
pilot = read_csv(paste0(dt.path, "/pilot-subIDs.csv"), show_col_types = F)
# save the excluded subjects minus the pilot participants
write(setdiff(exc, pilot$subID), file.path(dt.path, "PAL_exc.txt"))
# add pilot participants to the list
exc   = c(exc, pilot$subID)

# exclude these participants
df.tsk = df.tsk %>% filter(!(subID %in% exc)) %>%
  arrange(subID, trl)

# save data frames
saveRDS(df.tsk, file = file.path(dt.path, "df_tsk.rds"))

# how many are still left
length(unique(df.tsk$subID))
