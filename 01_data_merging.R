# Data prep for the climate-fitness manuscript

library(tidyverse)

# Load data ---------------------------------------------------------------

specimens_compiled = read_csv("rainfall-fitness/data/compiled_specimen_data.csv") %>% 
  mutate(species = str_c(genus, "_", species)) %>% 
  filter(!(species %in% c("s_barbatus", "s_cordatus", "s_morrisonii", "s_tortuosus"))) %>% 
  rename(collection_year = year, collection_day = day, collection_month = month) %>% 
  # mutate(collection_day_assumed = if_else(is.na(collection_day), 15, collection_day)) %>% 
  mutate(collection_date = make_date(year = collection_year, day = collection_day, month = collection_month)) %>% 
  mutate(sept_1_date = make_date(year = collection_year-1, month = 9, day = 1),
         collection_DOWY = as.numeric(collection_date - sept_1_date))
# Created with script "compile_all_data.R"
# 4862


re_levels = read_csv("rainfall-fitness/data/re_levels.csv") %>% 
  select(-coordUncertaintyAdjusted) %>% 
  filter(!(species %in% c("s_barbatus", "s_cordatus", "s_morrisonii", "s_tortuosus"))) %>% 
  distinct()
# created with script "generate_spatial_random_effects.R"
# 1639

flint_climate_summaries = read_csv("rainfall-fitness/data/climate_summaries.csv") %>% 
  filter(!(species %in% c("s_barbatus", "s_cordatus", "s_morrisonii", "s_tortuosus")))  
# Created with script "climate_summaries.R"

prism_data = read_csv("rainfall-fitness/data/climate_from_daily_data.csv")
# created with script "climate summaries from prism.R"

all_data = left_join(specimens_compiled, re_levels) %>% 
  left_join(., flint_climate_summaries) %>%
  mutate(species = str_sub(species, 3)) %>% 
  left_join(., prism_data)

write_csv(all_data, "rainfall-fitness/data/all_data_compiled.csv")


