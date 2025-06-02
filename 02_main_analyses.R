library(tidyverse)
library(brms)
library(phytools)
library(cowplot)
library(ggtree)
library(glue)
library(tidybayes)
library(modelr)
library(posterior)
library(tidybayes)
library(ggeffects)
library(ggrepel)
library(ggforce)
library(ggridges)
library(nlme)
library(ape)
library(caper)

# Load data ---------------------------------------------

dat = read_csv("rainfall-fitness/data/all_data_compiled.csv") 

mean_lats = dat %>% 
  group_by(species) %>% 
  summarize(mean_lat = mean(decimalLatitude))

species_list = unique(dat$species)

species_names = dat %>% 
  mutate(species_nice = str_c(toupper(genus), ". ", species)) %>% 
  dplyr::select(species, species_nice) %>% 
  distinct() %>% 
  mutate(species_phy = str_replace(species_nice, "C. ", "Caulanthus_"),
         species_phy = str_replace(species_phy, "S. ", "Streptanthus_")) %>% 
  mutate(species_phy_ordered = factor(species_phy, ordered = TRUE, 
                                      levels = rev(c("Streptanthus_diversifolius", 
                                                     "Streptanthus_polygaloides", 
                                                     "Caulanthus_hallii", 
                                                     "Caulanthus_amplexicaulis", 
                                                     "Streptanthus_drepanoides", 
                                                     "Streptanthus_breweri",  
                                                     "Caulanthus_simulans", 
                                                     "Streptanthus_insignis", 
                                                     "Streptanthus_glandulosus", 
                                                     "Caulanthus_anceps", 
                                                     "Caulanthus_coulteri", 
                                                     "Caulanthus_inflatus", 
                                                     "Caulanthus_cooperi", 
                                                     "Caulanthus_heterophyllus"))))


phy = read.newick("rainfall-fitness/data/tree_pruned.new")

pruned.tree = drop.tip(phy, phy$tip.label[-match(species_names$species_phy, phy$tip.label)])


# REPRODUCTION ----
# Q1A) Which components of climate influenced specimen fitness? ----------------

# Focus on a few potentially important predictors of reproductive output:
# Temperature Sept-May
# CWD Sept-May
# PPT Sept-May

## Prep main dataset ----

dat1 = dat %>% 
  filter(coordinateUncertaintyInMeters <= 5000) %>% 
  filter(entire_plant_confidence != "unsure") %>% 
  filter(!is.na(ppt_mm_rain_avg), 
         !is.na(tave_rain_avg), 
         !is.na(cwd_rain_avg), 
         !is.na(collection_DOWY), 
         !is.na(start_month_rain_avg)) %>% 
  # Comment out from here for sensitivity analysis #
  mutate(total_repro = bud + infBud + flwr + infFlwr + immFrt + fillFrt + infFrt + unfFrt + infImmFrt + unkRepStr,
         total_buds = bud + infBud,
         total_flowers = flwr + infFlwr,
         total_immFrts = immFrt + infImmFrt,
         total_matFrts = fillFrt + infFrt + unfFrt,
         total_frts = immFrt + infImmFrt + fillFrt + infFrt + unfFrt,
         weighted_fruits = round(0.33*total_flowers + 0.67*total_immFrts + 1*total_matFrts, 0),
         pheno_stage = weighted_fruits/total_repro) %>%
  # For sensitivity analysis with uncertain structures omitted and plants with >=25% uncertain structures omitted #
  # mutate(total_repro = bud + infBud + flwr + infFlwr + immFrt + fillFrt + infFrt + unfFrt + infImmFrt + unkRepStr) %>%
  # mutate(pct_uncertain  = (infBud + infFlwr + infImmFrt + infFrt + unkRepStr)/total_repro*100) %>%
  # filter(pct_uncertain < 25) %>%
  # mutate(total_repro = bud + flwr + immFrt + fillFrt + unfFrt,
  #        total_buds = bud,
  #        total_flowers = flwr,
  #        total_immFrts = immFrt,
  #        total_matFrts = fillFrt + unfFrt,
  #        total_frts = immFrt + fillFrt + unfFrt,
  #        weighted_fruits = round(0.33*total_flowers + 0.67*total_immFrts + 1*total_matFrts, 0),
  #        pheno_stage = weighted_fruits/total_repro) %>% 
  # End sensitivity section #
  mutate(start_month_rain_WY = if_else(start_month_rain_avg >= 9 , 
                                       start_month_rain_avg - 8, 
                                       start_month_rain_avg + 4)) %>% 
  dplyr::select(specimen, species, group, 
         total_repro, pheno_stage, weighted_fruits, 
         total_frts, total_matFrts,
         collection_DOWY, start_month_rain_WY, 
         cwd_rain_avg, tave_rain_avg, ppt_mm_rain_avg,
         first_rain_DOWY, 
         decimalLatitude, MAP_norm_sep_may, CWD_norm_sep_may # needed? maybe delete
         ) %>% 
  mutate(
    log_ppt_mm_rain_avg = log(ppt_mm_rain_avg),
    log_ppt_mm_rain_avg_scaled = (log_ppt_mm_rain_avg - mean(log_ppt_mm_rain_avg))/sd(log_ppt_mm_rain_avg),
    cwd_rain_avg_scaled = (cwd_rain_avg - mean(cwd_rain_avg))/sd(cwd_rain_avg),
    tave_rain_avg_scaled = (tave_rain_avg - mean(tave_rain_avg))/sd(tave_rain_avg),
    first_rain_month_scaled = (start_month_rain_WY - mean(start_month_rain_WY))/sd(start_month_rain_WY),
    collection_DOWY_scaled = (collection_DOWY - mean(collection_DOWY))/sd(collection_DOWY)) %>% 
  mutate(re_group = as.integer(as.factor(group))) %>% 
  left_join(., species_names)
       
dat1 %>% filter(total_frts > 0) # length 3283

dat1 %>% dplyr::select(specimen, species) %>% distinct %>% group_by(species) %>% summarize(n= n()) %>% arrange(n) %>% ungroup %>% summarize(mean(n))

dat1 %>% dplyr::select(specimen) %>% group_by(specimen) %>% summarize(n= n()) %>% ungroup %>% summarize(mean(n), min(n), max(n))

length(unique(dat1$specimen))
### Generate correlation plots for each species -----

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}


for (i in 1:length(species_list)) {

  species_choose = species_list[i]

  print(i)
  print(species_choose)

  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    dplyr::select(cwd_rain_avg, collection_DOWY,
           tave_rain_avg, log_ppt_mm_rain_avg,
           start_month_rain_WY)

  png(file = glue("rainfall-fitness/figs/predictor_correlations/{species_choose}.png"), height = 800, width = 800)
  pairs(dat_slim, upper.panel = panel.cor)
  dev.off()

}

## Model 1: temperature and precipitation on repro ----

mod_list1 = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept
  
  mod = brm(total_repro ~ collection_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + (1|re_group),
           family = negbinomial(),
           prior = priors,
           cores = 2,
           iter = 4000,
           init = 0,
           data = dat_slim, 
           refresh = 0,
           file = glue("rainfall-fitness/model_fits/repro_model1_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list1[[i]] = mod
  
}

names(mod_list1) = species_list

# look at default priors etc.
stancode(mod_list1[[1]])

### Make table of parameters and rhats ----

m1_all = map_df(mod_list1, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  mutate(model = "model1") %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`) 

# Check Rhats for convergence 
# sum should be 0
sum(m1_all$rhat >= 1.01)
# All converge with 4000 iterations

# Make a slimmed down table for point-range plots
param_list_m1 = m1_all %>% 
  filter(variable %in% c("b_tave_rain_avg_scaled", 
                         "b_log_ppt_mm_rain_avg_scaled",
                         "b_collection_DOWY_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "Average temperature\n(first rainy month until species'\naverage collection month)" =
                                      "b_tave_rain_avg_scaled",
                                    "Total precipitation\n(first rainy month until species'\naverage collection month)" =
                                      "b_log_ppt_mm_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"
                                    ),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "repro", model = "model1", model_eq = "total_repro ~ collection_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + (1|re_group)")


## Model 2: cwd on repro ----

mod_list2 = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept
  
  mod = brm(total_repro ~ collection_DOWY_scaled + cwd_rain_avg_scaled + (1|re_group),
          family = negbinomial(),
          prior = priors,
          cores = 2,
          iter = 4000,
          init = 0,
          data = dat_slim, 
          refresh = 0,
          file = glue("rainfall-fitness/model_fits/repro_model2_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list2[[i]] = mod
  
}

names(mod_list2) = species_list

### Make table of parameters and rhats ----

m2_all = map_df(mod_list2, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  mutate(model = "model2") %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`)

# Check Rhats for convergence 
# sums should all be 0
sum(m2_all$rhat >= 1.01)

param_list_m2 = m2_all %>% 
  filter(variable %in% c("b_cwd_rain_avg_scaled",
                         "b_collection_DOWY_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "Total climatic water deficit\n(first rainy month until species'\naverage collection month)" = "b_cwd_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "repro", model = "model2", model_eq = "total_repro ~ collection_DOWY_scaled + cwd_rain_avg_scaled + (1|re_group)") 


# Figure 3: repro vs. climate multipanel ----

plot_list1p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list1p[[i]] = mod_list1[[i]] %>%
    ggemmeans(., terms = c("log_ppt_mm_rain_avg_scaled")) %>%
    # ggemmeans(., terms = c("log_ppt_mm_rain_avg_scaled"), bias_correction = TRUE) %>%
    as_tibble() %>%
    mutate(species = species_choose)
  
  print(species_choose)
  
}

repro_ppt_panel_plot = plot_list1p %>%
  bind_rows() %>%
  mutate(ppt_mm_rain_avg = exp(x * sd(dat1$log_ppt_mm_rain_avg) + mean(dat1$log_ppt_mm_rain_avg))) %>%
  mutate(variable = "b_log_ppt_mm_rain_avg_scaled") %>%
  left_join(., param_list_m1) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = ppt_mm_rain_avg, y = total_repro), alpha = 0.2) +
  geom_line(aes(x = ppt_mm_rain_avg, y = predicted, color = sig)) +
  geom_ribbon(aes(x = ppt_mm_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Total reproductive structures") +
  scale_x_log10(name = "Precipitation (mm), log scale") +
  guides(color = "none", fill = "none")+
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)); repro_ppt_panel_plot



p = ggtree(pruned.tree) + 
  # geom_tiplab(hjust = 0, offset = 0.0001) +
  # xlim(0, 0.0027) +
  xlim(0, 0.00155) +
  theme_tree(); p

ests1 = ggplot(data = filter(param_list_m1, variable != "b_collection_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_x_continuous(limits = c(-0.9, 0.9)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nreproduction") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests1

ests2 = ggplot(data = filter(param_list_m2, variable != "b_collection_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_x_continuous(limits = c(-0.9, 0.9)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nreproduction") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()); ests2

p_centered = plot_grid(NULL, p, NULL, rel_heights = c(0.1, 1, 0.1), ncol = 1)

p_all = plot_grid(p_centered, ests1, NULL, ests2, rel_widths = c(0.3, 1, 0.05, 0.38), rel_heights = c(0.8, 1, 1), nrow = 1,
                  labels = c("", "b.", "c."), label_x = c(0, 0.25, 0)); p_all

plot_grid(repro_ppt_panel_plot, p_all, rel_heights = c(0.8, 1, 1), nrow = 2,
          labels = c("a.", "", ""), label_x = c(0, 0.25, 0))

ggsave("rainfall-fitness/figs/Fig3_repro_mod1_multipanel.pdf", height = 10, width = 9)
ggsave("rainfall-fitness/figs/Fig3_repro_mod1_multipanel.png", height = 10, width = 9)



# Figure S10: temp and cwd on repro panel plots ----

# Create list to put model predictions in with ggemmeans
plot_list1t = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list1t[[i]] = mod_list1[[i]] %>%
    ggemmeans(., terms = c("tave_rain_avg_scaled")) %>%
    # ggemmeans(., terms = c("tave_rain_avg_scaled"), bias_correction = TRUE) %>%
    as_tibble() %>%
    mutate(species = species_choose)
  
  print(species_choose)
  
}

repro_temp_panel_plot = plot_list1t %>%
  bind_rows() %>%
  mutate(tave_rain_avg = x * sd(dat1$tave_rain_avg) + mean(dat1$tave_rain_avg)) %>%
  mutate(variable = "b_tave_rain_avg_scaled") %>%
  left_join(., param_list_m1) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = tave_rain_avg, y = total_repro), alpha = 0.2) +
  geom_line(aes(x = tave_rain_avg, y = predicted, color = sig), size = 1) +
  geom_ribbon(aes(x = tave_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  guides(color = "none", fill = "none")+
  labs(y = "Total reproductive structures", x = "Average temperature (C)") +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text(hjust = 0, vjust = 0.5)); repro_temp_panel_plot


plot_list2 = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {

  species_choose = species_list[i]

  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))

  plot_list2[[i]] = mod_list2[[i]] %>%
    ggemmeans(., terms = c("cwd_rain_avg_scaled")) %>%
    # ggemmeans(., terms = c("cwd_rain_avg_scaled"), bias_correction = TRUE) %>%
    as_tibble() %>%
    mutate(species = species_choose)

  print(species_choose)
  
}

repro_cwd_panel_plot = plot_list2 %>%
  bind_rows() %>%
  mutate(cwd_rain_avg = x * sd(dat1$cwd_rain_avg) + mean(dat1$cwd_rain_avg)) %>%
  mutate(variable = "b_cwd_rain_avg_scaled") %>%
  filter(cwd_rain_avg > 0) %>%
  left_join(., param_list_m2) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = cwd_rain_avg, y = total_repro), alpha = 0.2) +
  geom_line(aes(x = cwd_rain_avg, y = predicted, color = sig)) +
  geom_ribbon(aes(x = cwd_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Total reproductive structures", color = "Species", fill = "Species") +
  scale_x_continuous(name = "Climatic water deficit") +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  guides(color = "none", fill = "none")+
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), 
        axis.text.x = element_text( hjust = 0, vjust = 0.5)); repro_cwd_panel_plot

plot_grid(repro_temp_panel_plot, repro_cwd_panel_plot, rel_heights = c(1,1), nrow = 2,
          labels = c("a.", "b."))

ggsave("rainfall-fitness/figs/FigS10_repro_mod12_cwd_temp_panel.pdf", height = 10, width = 9)
ggsave("rainfall-fitness/figs/FigS10_repro_mod12_cwd_temp_panel.png", height = 10, width = 9)



# Q2A) Does later onset of seasonal rainfall, predicted to become more frequent under climate change, reduce fitness? ----

## Prep daily dataset -------------------

dat3 = dat %>% 
  filter(!is.na(first_rain_DOWY), !is.na(tave_sep_may)) %>% 
  filter(coordinateUncertaintyInMeters <= 5000) %>% 
  filter(entire_plant_confidence != "unsure") %>% 
  filter(first_rain_DOWY < 182) %>% # filter out eight specimens with the first rain event after March 1
  filter(first_rain_DOWY < collection_DOWY) %>% # filter out one specimen with negative growth window
  mutate(start_month_rain_WY = if_else(start_month_rain_avg >= 9 , start_month_rain_avg - 8, start_month_rain_avg + 4)) %>% 
  mutate(total_repro = bud + infBud + flwr + infFlwr + immFrt + fillFrt + infFrt + unfFrt + infImmFrt + unkRepStr,
         total_buds = bud + infBud,
         total_flowers = flwr + infFlwr,
         total_immFrts = immFrt + infImmFrt,
         total_matFrts = fillFrt + infFrt + unfFrt,
         weighted_fruits = round(0.33*total_flowers + 0.67*total_immFrts + 1*total_matFrts, 0),
         pheno_stage = weighted_fruits/total_repro) %>% 
  dplyr::select(specimen, species, first_rain_DOWY, collection_DOWY, late_herbivory, stem_diam_cm, group, total_repro,
         ppt_mm_sep_may, cwd_sep_may, tave_sep_may, total_rain_events, total_ppt_mm_prism, GDD, tmean_first_rain_to_coll, start_month_rain_WY, weighted_fruits,
         ppt_mm_rain_avg, cwd_rain_avg, tave_rain_avg, pheno_stage) %>% 
  mutate(first_rain_DOWY_scaled = (first_rain_DOWY - mean(first_rain_DOWY))/sd(first_rain_DOWY),
         collection_DOWY_scaled = (collection_DOWY - mean(collection_DOWY))/sd(collection_DOWY),
         log_ppt_mm_sep_may = log(ppt_mm_sep_may),
         ppt_mm_sep_may_scaled = (ppt_mm_sep_may - mean(ppt_mm_sep_may))/sd(ppt_mm_sep_may),
         log_ppt_mm_sep_may_scaled = (log_ppt_mm_sep_may - mean(log_ppt_mm_sep_may))/sd(log_ppt_mm_sep_may),
         cwd_sep_may_scaled = (cwd_sep_may - mean(cwd_sep_may))/sd(cwd_sep_may),
         tave_sep_may_scaled = (tave_sep_may - mean(tave_sep_may))/sd(tave_sep_may),
         log_ppt_mm_rain_avg = log(ppt_mm_rain_avg),
         ppt_mm_rain_avg_scaled = (ppt_mm_rain_avg - mean(ppt_mm_rain_avg))/sd(ppt_mm_rain_avg),
         log_ppt_mm_rain_avg_scaled = (log_ppt_mm_rain_avg - mean(log_ppt_mm_rain_avg))/sd(log_ppt_mm_rain_avg),
         cwd_rain_avg_scaled = (cwd_rain_avg - mean(cwd_rain_avg))/sd(cwd_rain_avg),
         tave_rain_avg_scaled = (tave_rain_avg - mean(tave_rain_avg))/sd(tave_rain_avg),
         total_rain_events_scaled = (total_rain_events - mean(total_rain_events))/sd(total_rain_events),
         log_total_ppt_mm_prism = log(total_ppt_mm_prism),
         first_rain_month_scaled = (start_month_rain_WY - mean(start_month_rain_WY))/sd(start_month_rain_WY),
         log_total_ppt_mm_prism_scaled = (log_total_ppt_mm_prism - mean(log_total_ppt_mm_prism))/sd(log_total_ppt_mm_prism),
         tmean_first_rain_to_coll_scaled = (tmean_first_rain_to_coll - mean(tmean_first_rain_to_coll))/sd(tmean_first_rain_to_coll),
         collection_DOWY_scaled = (collection_DOWY - mean(collection_DOWY))/sd(collection_DOWY),
         GDD_scaled = (GDD-mean(GDD))/sd(GDD)) %>% 
  mutate(re_group = as.integer(as.factor(group))) %>% 
  left_join(., species_names) 

length(unique(dat3$specimen))

### Generate correlation plots for each species -----

for (i in 1:length(species_list)) {

  species_choose = species_list[i]

  print(i)
  print(species_choose)

  dat_slim = dat3 %>%
    filter(species == species_choose) %>%
    dplyr::select(collection_DOWY, first_rain_DOWY, cwd_rain_avg,
           log_ppt_mm_rain_avg,
           tave_rain_avg,
           start_month_rain_WY)

  png(file = glue("rainfall-fitness/figs/predictor_correlations/{species_choose}_daily_dataset.png"), height = 800, width = 800)
  pairs(dat_slim, upper.panel = panel.cor)
  dev.off()

}


## Model 3: precipitation timing -----

mod_list3 = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat3 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept

  mod = brm(total_repro ~ collection_DOWY_scaled + first_rain_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + (1|re_group),
           family = negbinomial(),
           prior = priors,
           cores = 2,
           iter = 4000,
           init = 0,
           data = dat_slim, 
  refresh = 0,
  file = glue("rainfall-fitness/model_fits/repro_model3_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list3[[i]] = mod
  
}

names(mod_list3) = species_list


### Make table of parameters and rhats ----

m3_all = map_df(mod_list3, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  mutate(model = "model3") %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`)

# Check Rhats for convergence 
# sums should all be 0
sum(m3_all$rhat >= 1.01)

param_list_m3 = m3_all %>% 
  filter(variable %in% c("b_first_rain_DOWY_scaled", 
                         "b_tave_rain_avg_scaled", 
                         "b_log_ppt_mm_rain_avg_scaled", 
                         "b_collection_DOWY_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "First rain event DOWY" = "b_first_rain_DOWY_scaled",
                                    "Average temperature\n(first rainy month until species'\naverage collection month)" = "b_tave_rain_avg_scaled",
                                    "Total precipitation\n(first rainy month until species'\naverage collection month)" = "b_log_ppt_mm_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "repro", model = "model3", model_eq = "total_repro ~ collection_DOWY_scaled + first_rain_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + (1|re_group)")



# Figure S11: point ranges model 3 on repro ----

ests3 = ggplot(data = param_list_m3) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  # scale_x_continuous(limits = c(-0.85, 0.85)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nreproduction") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests3

ggsave("rainfall-fitness/figs/FigS11_repro_mod3_estimates.pdf", height = 6, width = 9)
ggsave("rainfall-fitness/figs/FigS11_repro_mod3_estimates.png", height = 6, width = 9)

# Figure 4: first rain vs. repro -----

plot_list3 = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat3 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list3[[i]] = mod_list3[[i]] %>%
    ggemmeans(., terms = c("first_rain_DOWY_scaled")) %>% 
    # ggemmeans(., terms = c("first_rain_DOWY_scaled"), bias_correction = TRUE) %>% 
    as_tibble() %>%
    mutate(species = species_choose)
  
}

# sept 30 = 30
# oct 31 = 61
# nov 30  = 91
# dec 31 = 122 
# jan 31 = 153
# feb 28 = 181
# march 31 = 212
# april 30 = 242
# may 31 = 273
# june 30 = 303
# July 31 = 334
# aug 31 = 365
# sept 30 = 395
# oct 31 = 426

repro_germ_panel_plot_tall = plot_list3 %>%
  bind_rows() %>%
  mutate(first_rain_DOWY = x * sd(dat3$first_rain_DOWY) + mean(dat3$first_rain_DOWY)) %>%
  mutate(variable = "b_first_rain_DOWY_scaled") %>%
  left_join(., species_names) %>%
  left_join(., param_list_m3) %>% # add on significance to color lines
  ggplot() +
  geom_point(data = dat3, aes(x = first_rain_DOWY, y = total_repro), alpha = 0.2) +
  geom_line(aes(x = first_rain_DOWY, y = predicted, color = sig)) +
  geom_ribbon(aes(x = first_rain_DOWY, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Total reproductive structures") +
  scale_x_continuous(name = "Timing of first rain event", breaks = c(0, 31, 62, 92, 123, 154, 182), labels = c("Sep 1", "Oct 1", "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1")) +
  scale_color_manual(values = c("grey60", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "dodgerblue")) +
  guides(color = "none", fill = "none")+
  facet_wrap(. ~ species_nice, nrow = 5, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text( hjust = 0, vjust = 0.5, angle = 90)); repro_germ_panel_plot_tall

ests3b = ggplot(data = filter(param_list_m3, variable == "b_first_rain_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "dodgerblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nreproduction") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests3b

p1 = plot_grid(p_centered, ests3b, rel_widths = c(0.3, 1), rel_heights = c(0.8, 1), nrow = 1, labels = c("", "a."), label_x = c(0, 0.25)); p1

p2 = plot_grid(p1, repro_germ_panel_plot_tall, rel_widths = c(0.8, 1), rel_heights = c(0.8, 1), nrow = 1, labels = c("", "b.")); p2

ggsave("rainfall-fitness/figs/Fig4_repro_mod3_estimates_multipanel.pdf", height = 6, width = 12)
ggsave("rainfall-fitness/figs/Fig4_repro_mod3_estimates_multipanel.png", height = 6, width = 12)


# PHENOLOGY ----
# Q1B) Which components of climate influenced specimen phenology? ----------------

## Model 1: temperature and precipitation on pheno ----

mod_list1_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept
  
  mod = brm(weighted_fruits|trials(total_repro) ~ 
              collection_DOWY_scaled + 
              tave_rain_avg_scaled + 
              log_ppt_mm_rain_avg_scaled + 
              (1|re_group), 
            family = binomial(),
            prior = priors,
            cores = 2,
            iter = 4000,
            init = 0,
            data = dat_slim,
            refresh = 0,
            file = glue("rainfall-fitness/model_fits/pheno_model1_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list1_p[[i]] = mod
  
}

names(mod_list1_p) = species_list

### Make table of parameters and rhats ----

m1_all_p = map_df(mod_list1_p, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`) 

# Check Rhats for convergence 
# sum should be 0
sum(m1_all_p$rhat >= 1.01)
# All converge with 4000 iterations

# Make a slimmed down table for point-range plots

param_list_m1_p = m1_all_p %>% 
  dplyr::select(-species_nice, -species_phy) %>%
  filter(variable %in% c( "b_tave_rain_avg_scaled", "b_log_ppt_mm_rain_avg_scaled", "b_collection_DOWY_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "Average temperature\n(first rainy month until species'\naverage collection month)" =
                                      "b_tave_rain_avg_scaled",
                                    "Total precipitation\n(first rainy month until species'\naverage collection month)" =
                                      "b_log_ppt_mm_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "pheno", model = "model1", model_eq = "weighted_fruits|trials(total_repro) ~ collection_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + (1|re_group)")


## Model 2: CWD ----

mod_list2_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept
  
  mod = brm(weighted_fruits|trials(total_repro) ~ 
              collection_DOWY_scaled + 
              cwd_rain_avg_scaled + 
              (1|re_group), 
            family = binomial(),
            prior = priors,
            cores = 2,
            iter = 4000,
            init = 0,
            data = dat_slim,
            refresh = 0,
            file = glue("rainfall-fitness/model_fits/pheno_model2_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list2_p[[i]] = mod
  
}

names(mod_list2_p) = species_list


### Make table of parameters and rhats ----

m2_all_p = map_df(mod_list2_p, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`)

# Check Rhats for convergence 
# sums should all be 0
sum(m2_all_p$rhat >= 1.01)
# all good but anceps ess are low

param_list_m2_p = m2_all_p %>% 
  dplyr::select(-species_nice, -species_phy) %>%
  filter(variable %in% c("b_cwd_rain_avg_scaled", "b_collection_DOWY_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "Total climatic water deficit\n(first rainy month until species'\naverage collection month)" = "b_cwd_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "pheno", model = "model2", model_eq = "weighted_fruits|trials(total_repro) ~ collection_DOWY_scaled + cwd_rain_avg_scaled + (1|re_group)")


# Figure 1: pheno vs. climate multipanel ----


plot_list1p_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list1p_p[[i]] = mod_list1_p[[i]] %>%
    ggemmeans(., terms = c("log_ppt_mm_rain_avg_scaled [all]")) %>%
    # ggemmeans(., terms = c("log_ppt_mm_rain_avg_scaled [all]"), bias_correction = TRUE) %>% 
    as_tibble() %>%
    mutate(species = species_choose)
  
}

pheno_ppt_panel_plot = plot_list1p_p %>%
  bind_rows() %>%
  mutate(ppt_mm_rain_avg = exp(x * sd(dat1$log_ppt_mm_rain_avg) + mean(dat1$log_ppt_mm_rain_avg))) %>%
  mutate(variable = "b_log_ppt_mm_rain_avg_scaled") %>%
  left_join(., param_list_m1_p) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = ppt_mm_rain_avg, y = pheno_stage), alpha = 0.2) +
  geom_line(aes(x = ppt_mm_rain_avg, y = predicted, color = sig)) +
  geom_ribbon(aes(x = ppt_mm_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Phenological stage") +
  scale_x_log10(name = "Precipitation (mm), log scale") +
  guides(color = "none", fill = "none")+
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)); pheno_ppt_panel_plot


ests1_p = ggplot(data = filter(param_list_m1_p, variable != "b_collection_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nphenology") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests1_p

ests2_p = ggplot(data = filter(param_list_m2_p, variable != "b_collection_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nphenology") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()); ests2_p

p_centered = plot_grid(NULL, p, NULL, rel_heights = c(0.1, 1, 0.1), ncol = 1)

p_all_p = plot_grid(p_centered, ests1_p, NULL, ests2_p, rel_widths = c(0.3, 1, 0.05, 0.38), rel_heights = c(0.8, 1, 1), nrow = 1,
                    labels = c("", "b.", "c."), label_x = c(0, 0.25, 0)); p_all_p

plot_grid(pheno_ppt_panel_plot, p_all_p, rel_heights = c(0.8, 1, 1), nrow = 2,
          labels = c("a.", "", ""), label_x = c(0, 0.25, 0))

ggsave("rainfall-fitness/figs/Fig1_pheno_mod1_multipanel.pdf", height = 10, width = 9)
ggsave("rainfall-fitness/figs/Fig1_pheno_mod1_multipanel.png", height = 10, width = 9)



# Figure S8: temp and cwd on pheno panel plot ----

plot_list1t_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list1t_p[[i]] = mod_list1_p[[i]] %>%
    ggemmeans(., terms = c("tave_rain_avg_scaled [all]")) %>%
    # ggemmeans(., terms = c("tave_rain_avg_scaled [all]", bias_correction = TRUE)) %>% 
    as_tibble() %>%
    mutate(species = species_choose)
  
}

pheno_temp_panel_plot = plot_list1t_p %>%
  bind_rows() %>%
  mutate(tave_rain_avg = x * sd(dat1$tave_rain_avg) + mean(dat1$tave_rain_avg)) %>%
  mutate(variable = "b_tave_rain_avg_scaled") %>%
  left_join(., param_list_m1_p) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = tave_rain_avg, y = pheno_stage), alpha = 0.2) +
  # geom_smooth(data = dat1, aes(x = tave_rain_avg, y = total_repro), method = "lm", alpha = 0.2) +
  geom_line(aes(x = tave_rain_avg, y = predicted, color = sig), size = 1) +
  geom_ribbon(aes(x = tave_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  guides(color = "none", fill = "none")+
  labs(y = "Phenological stage", x = "Average temperature (C)") +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text(hjust = 0, vjust = 0.5)); pheno_temp_panel_plot


plot_list2_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {

  species_choose = species_list[i]

  dat_slim = dat1 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))

  plot_list2_p[[i]] = mod_list2_p[[i]] %>%
    ggemmeans(., terms = c("cwd_rain_avg_scaled [all]")) %>% 
    # ggemmeans(., terms = c("cwd_rain_avg_scaled [all]"), bias_correction = TRUE) %>% 
    as_tibble() %>%
    mutate(species = species_choose)

}

pheno_cwd_panel_plot = plot_list2_p %>%
  bind_rows() %>%
  mutate(cwd_rain_avg = x * sd(dat1$cwd_rain_avg) + mean(dat1$cwd_rain_avg)) %>%
  mutate(variable = "b_cwd_rain_avg_scaled") %>%
  filter(cwd_rain_avg > 0) %>%
  left_join(., param_list_m2_p) %>% # add on significance to color lines
  left_join(., species_names) %>%
  ggplot() +
  geom_point(data = dat1, aes(x = cwd_rain_avg, y = pheno_stage), alpha = 0.2) +
  geom_line(aes(x = cwd_rain_avg, y = predicted, color = sig)) +
  geom_ribbon(aes(x = cwd_rain_avg, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Phenological stage", color = "Species", fill = "Species") +
  scale_x_continuous(name = "Climatic water deficit") +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  guides(color = "none", fill = "none")+
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text( hjust = 0, vjust = 0.5)); pheno_cwd_panel_plot


plot_grid(pheno_temp_panel_plot, pheno_cwd_panel_plot, rel_heights = c(1,1), nrow = 2,
          labels = c("a.", "b."))

ggsave("rainfall-fitness/figs/FigS8_pheno_mod12_cwd_temp_panel.pdf", height = 10, width = 9)
ggsave("rainfall-fitness/figs/FigS8_pheno_mod12_cwd_temp_panel.png", height = 10, width = 9)



# Q2B) Does later onset of seasonal rainfall, predicted to become more frequent under climate change, change phenology? ----

## Model 3: rainfall timing with recent dataset ----

mod_list3_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  print(i)
  print(species_choose) 
  
  dat_slim = dat3 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  priors <- c(prior(normal(0, 1), class = b))
  # Use default student t prior on intercept
  
  
  mod = brm(weighted_fruits|trials(total_repro) ~ 
              collection_DOWY_scaled + 
              tave_rain_avg_scaled + 
              log_ppt_mm_rain_avg_scaled + 
              first_rain_DOWY_scaled + 
              (1|re_group),
            family = binomial(),
            prior = priors,
            cores = 2,
            iter = 4000,
            init = 0,
            data = dat_slim,
            refresh = 0,
            file = glue("rainfall-fitness/model_fits/pheno_model3_{species_choose}.rds"))
  
  print(summary(mod))
  
  mod_list3_p[[i]] = mod
  
}

names(mod_list3_p) = species_list


### Make table of parameters and rhats ----

m3_all_p = map_df(mod_list3_p, summarise_draws, .id="species", ~quantile(.x, probs = c(0.025, 0.975)), mean, median, sd, rhat, ess_bulk, ess_tail) %>% 
  left_join(species_names) %>% 
  rename(lower = `2.5%`, upper = `97.5%`)

# Check Rhats for convergence 
# sums should all be 0
sum(m3_all_p$rhat >= 1.01)

param_list_m3_p = m3_all_p %>% 
  dplyr::select(-species_nice, -species_phy) %>%
  filter(variable %in% c("b_collection_DOWY_scaled", "b_first_rain_DOWY_scaled", "b_tave_rain_avg_scaled", "b_log_ppt_mm_rain_avg_scaled")) %>% 
  mutate(variable_nice = fct_recode(variable,
                                    "First rain event DOWY" = "b_first_rain_DOWY_scaled",
                                    "Average temperature\n(first rainy month until species'\naverage collection month)" = "b_tave_rain_avg_scaled",
                                    "Total precipitation\n(first rainy month until species'\naverage collection month)" = "b_log_ppt_mm_rain_avg_scaled",
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%
  mutate(response = "pheno", model = "model3", model_eq = "weighted_fruits|trials(total_repro) ~ collection_DOWY_scaled + tave_rain_avg_scaled + log_ppt_mm_rain_avg_scaled + first_rain_DOWY_scaled + (1|re_group)")


# Figure S9: point ranges model 3 on pheno ----

ests3_p = ggplot(data = param_list_m3_p) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  # scale_x_continuous(limits = c(-0.85, 0.85)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nphenology") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests3_p

ggsave("rainfall-fitness/figs/FigS9_pheno_mod3_estimates.pdf", height = 6, width = 9)
ggsave("rainfall-fitness/figs/FigS9_pheno_mod3_estimates.png", height = 6, width = 9)


# Figure 2: first rain vs. pheno ----


plot_list3_p = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat3 %>%
    filter(species == species_choose) %>%
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_list3_p[[i]] = mod_list3_p[[i]] %>%
    ggemmeans(., terms = c("first_rain_DOWY_scaled [all]")) %>% 
    # ggemmeans(., terms = c("first_rain_DOWY_scaled [all]"), bias_correction = TRUE) %>% 
    as_tibble() %>%
    mutate(species = species_choose)
  
}

pheno_germ_panel_plot = plot_list3_p %>%
  bind_rows() %>%
  mutate(first_rain_DOWY = x * sd(dat3$first_rain_DOWY) + mean(dat3$first_rain_DOWY)) %>%
  mutate(variable = "b_first_rain_DOWY_scaled") %>%
  left_join(., species_names) %>%
  left_join(., param_list_m3_p) %>% # add on significance to color lines
  ggplot() +
  geom_point(data = dat3, aes(x = first_rain_DOWY, y = pheno_stage), alpha = 0.2) +
  geom_line(aes(x = first_rain_DOWY, y = predicted, color = sig)) +
  geom_ribbon(aes(x = first_rain_DOWY, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Phenological stage") +
  scale_x_continuous(name = "Timing of first rain event", breaks = c(0, 31, 62, 92, 123, 154, 182), labels = c("Sep 1", "Oct 1", "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1")) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  guides(color = "none", fill = "none")+
  facet_wrap(. ~ species_nice, nrow = 3, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text( hjust = 0, vjust = 0.5, angle = 90)); pheno_germ_panel_plot

# sept 30 = 30
# oct 31 = 61
# nov 30  = 91
# dec 31 = 122 
# jan 31 = 153
# feb 28 = 181
# march 31 = 212
# april 30 = 242
# may 31 = 273
# june 30 = 303
# July 31 = 334
# aug 31 = 365
# sept 30 = 395
# oct 31 = 426

ests3b_p = ggplot(data = filter(param_list_m3_p, variable == "b_first_rain_DOWY_scaled")) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  facet_grid(.~ variable_nice) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  # scale_x_continuous(limits = c(-0.85, 0.85)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect on\nphenology") +
  theme(axis.text.y = element_text(face = "italic", hjust = 0.5), axis.title.y = element_blank()); ests3b_p

p_centered = plot_grid(NULL, p, NULL, rel_heights = c(0.05, 1, 0.1), ncol = 1)

p1_p = plot_grid(p_centered, ests3b_p, rel_widths = c(0.3, 1), rel_heights = c(0.8, 1), nrow = 1, labels = c("", "a."), label_x = c(0, 0.25)); p1_p

pheno_germ_panel_plot_tall = plot_list3_p %>%
  bind_rows() %>%
  mutate(first_rain_DOWY = x * sd(dat3$first_rain_DOWY) + mean(dat3$first_rain_DOWY)) %>%
  mutate(variable = "b_first_rain_DOWY_scaled") %>%
  left_join(., species_names) %>%
  left_join(., param_list_m3_p) %>% # add on significance to color lines
  ggplot() +
  geom_point(data = dat3, aes(x = first_rain_DOWY, y = pheno_stage), alpha = 0.2) +
  geom_line(aes(x = first_rain_DOWY, y = predicted, color = sig)) +
  geom_ribbon(aes(x = first_rain_DOWY, ymin = conf.low, ymax = conf.high, fill = sig), alpha = 0.6) +
  labs(y = "Phenological stage") +
  scale_x_continuous(name = "Timing of first rain event", breaks = c(0, 31, 62, 92, 123, 154, 182), labels = c("Sep 1", "Oct 1", "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1")) +
  scale_color_manual(values = c("grey60", "tomato", "dodgerblue")) +
  scale_fill_manual(values = c("grey60", "tomato", "dodgerblue")) +
  guides(color = "none", fill = "none")+
  facet_wrap(. ~ species_nice, nrow = 5, scales = "free_y") +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text = element_text(face = "italic"), axis.text.x = element_text( hjust = 0, vjust = 0.5, angle = 90)); pheno_germ_panel_plot_tall

p2_p = plot_grid(p1_p, pheno_germ_panel_plot_tall, rel_widths = c(0.8, 1), rel_heights = c(0.8, 1), nrow = 1, labels = c("", "b.")); p2_p

ggsave("rainfall-fitness/figs/Fig2_pheno_mod3_estimates_multipanel.pdf", height = 6, width = 12)
ggsave("rainfall-fitness/figs/Fig2_pheno_mod3_estimates_multipanel.png", height = 6, width = 12)


# Appendix 1 ----

all_params = bind_rows(param_list_m1, param_list_m2, param_list_m3, param_list_m1_p, param_list_m2_p, param_list_m3_p)

table(all_params$model, all_params$variable)
write_csv(all_params, "rainfall-fitness/results/Appendix1_all_parameter_estimates.csv")


## Q3) Across species, are phenological and estimated reproduction responses to climate variables correlated? ----

all_slopes = bind_rows(param_list_m1, param_list_m1_p, param_list_m2, param_list_m2_p, param_list_m3_p, param_list_m3) %>% 
  pivot_wider(id_cols = c(species, variable, model), names_from = response, values_from = mean)

slopes_m12 = all_slopes %>% filter(model %in% c("model1", "model2"))
slopes_m3 = all_slopes %>% filter(model %in% c("model3"))

pic_results = list()
# phy indep contrasts:
# http://lukejharmon.github.io/ilhabela/instruction/2015/07/02/phylogenetic-independent-contrasts/
bm = corBrownian(1, pruned.tree)
bm

# Temp
temp_slopes = filter(slopes_m12, variable == "b_tave_rain_avg_scaled") %>% 
  left_join(., species_names) %>% 
  as.data.frame()
cor.test(temp_slopes$repro, temp_slopes$pheno)

ix = pic(temp_slopes$pheno, pruned.tree)
iy = pic(temp_slopes$repro, pruned.tree)
fit = lm(iy~ix-1) ## we have to fit the model without an intercept term
summary(fit)
pic_results[[1]] = as_tibble(summary(fit)$coefficients) %>% 
  mutate(variable = "temp")

pic_temp = ggplot(temp_slopes, aes(y = repro, x = pheno)) +
  # ggtitle("Temperature") +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Slope of phenological\nresponse to temperature", y = "Slope of reproductive\nresponse to temperature") +
  geom_text_repel(aes(y = repro, x = pheno, label = species), size = 3); pic_temp


# PPT
ppt_slopes = filter(slopes_m12, variable == "b_log_ppt_mm_rain_avg_scaled") %>% 
  left_join(., species_names) %>% 
  as.data.frame()
cor.test(ppt_slopes$repro, ppt_slopes$pheno)

ix = pic(ppt_slopes$pheno, pruned.tree)
iy = pic(ppt_slopes$repro, pruned.tree)
fit = lm(iy~ix-1) ## we have to fit the model without an intercept term
summary(fit)
pic_results[[2]] = as_tibble(summary(fit)$coefficients) %>% 
  mutate(variable = "ppt")

pic_ppt = ggplot(ppt_slopes, aes(y = repro, x = pheno)) +
  # ggtitle("Precipitation") +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Slope of phenological\nresponse to precipitation", y = "Slope of reproductive\nresponse to precipitation") +
  geom_text_repel(aes(y = repro, x = pheno, label = species), size = 3); pic_ppt

# CWD
cwd_slopes = filter(slopes_m12, variable == "b_cwd_rain_avg_scaled") %>% 
  left_join(., species_names) %>% 
  as.data.frame()
cor.test(cwd_slopes$repro, cwd_slopes$pheno)

ix = pic(cwd_slopes$pheno, pruned.tree)
iy = pic(cwd_slopes$repro, pruned.tree)
fit<-lm(iy~ix-1) ## we have to fit the model without an intercept term
summary(fit)
pic_results[[3]] = as_tibble(summary(fit)$coefficients) %>% 
  mutate(variable = "CWD")

pic_cwd = ggplot(cwd_slopes, aes(y = repro, x = pheno)) +
  # ggtitle("Climatic water deficit") +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Slope of phenological\nresponse to CWD", y = "Slope of reproductive\nresponse to CWD") +
  geom_text_repel(aes(y = repro, x = pheno, label = species), size = 3); pic_cwd

# Germ timing
germ_slopes = filter(slopes_m3, variable == "b_first_rain_DOWY_scaled") %>% 
  left_join(., species_names) %>% 
  as.data.frame()
cor.test(germ_slopes$repro, germ_slopes$pheno)

ix = pic(germ_slopes$pheno, pruned.tree)
iy = pic(germ_slopes$repro, pruned.tree)
fit<-lm(iy~ix-1) ## we have to fit the model without an intercept term
summary(fit)

pic_results[[4]] = as_tibble(summary(fit)$coefficients) %>% 
  mutate(variable = "first_rain_DOY")

pic_germ = ggplot(germ_slopes, aes(y = repro, x = pheno)) +
  # ggtitle("First rain DOY") +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Slope of phenological\nresponse to first rain DOWY", y = "Slope of reproductive\nresponse to first rain DOWY") +
  geom_text_repel(aes(y = repro, x = pheno, label = species), size = 3); pic_germ

pic_df = pic_results %>% bind_rows()

plot_grid(pic_temp, pic_ppt, pic_cwd, pic_germ, 
          ncol = 2,
          labels = c("a.", "b.", "c.", "d."),
          align = "hv")
          # labels = c("a. Temperature", "b. Precipitation", "c. Climatic water deficit", "d. First rain event DOWY"))

ggsave("rainfall-fitness/figs/FigS12_pheno_v_repro.pdf", width = 8, height = 7.5)
ggsave("rainfall-fitness/figs/FigS12_pheno_v_repro.png", width = 8, height = 7.5)

write_csv(pic_df, "rainfall-fitness/results/TableS2_pics.csv")

## Q4) Are responses to climate evolutionarily conserved across the clade? ----

### Write out model slopes and ses for phylogenetic analyses -----

mods = ls(pattern  = "mod_list*")

param_list = vector(mode = "list", length = length(mods))
param_list[[1]] = vector(mode = "list", length = length(species_list))
param_list[[2]] = vector(mode = "list", length = length(species_list))
param_list[[3]] = vector(mode = "list", length = length(species_list))
param_list[[4]] = vector(mode = "list", length = length(species_list))
param_list[[5]] = vector(mode = "list", length = length(species_list))
param_list[[6]] = vector(mode = "list", length = length(species_list))


for (j in 1:length(mods)) {
  
  mod = mget(mods[j])[[1]]
  name = mods[j]
  
  for (i in 1:length(species_list)){
    
    fe = summary(mod[[i]])$fixed %>% 
      as_tibble(rownames = "parameter") %>% 
      mutate(species = species_list[i],
             model = name)
    re = summary(mod[[i]])$random$re_group %>% 
      as_tibble(rownames = "parameter") %>% 
      mutate(species = species_list[i],
             model = name)
    sp = summary(mod[[i]])$spec_pars %>% 
      as_tibble(rownames = "parameter") %>% 
      mutate(species = species_list[i],
             model = name)
    
    param_list[[j]][[i]] = bind_rows(fe, re, sp)
    
  }
}

param_list_all = bind_rows(param_list) %>% 
  rename(upper = `u-95% CI`, 
         lower = `l-95% CI`) %>% 
  left_join(species_names) %>% 
  unite(parameter_model, parameter, model, remove = FALSE)

param_list_for_phy_se = param_list_all %>% 
  dplyr::select(species, species_nice, species_phy, Est.Error, parameter_model) %>% 
  pivot_wider(names_from = parameter_model, values_from = Est.Error) %>%
  dplyr::select(-contains("Intercept"), -contains("sigma"), -contains("shape"))

param_list_for_phy = param_list_all %>% 
  dplyr::select(species, species_nice, species_phy, Estimate, parameter_model) %>% 
  pivot_wider(names_from = parameter_model, values_from = Estimate) %>% 
  dplyr::select(-contains("Intercept"), -contains("sigma"), -contains("shape"))

names(param_list_for_phy_se) == names(param_list_for_phy)

param_list_for_phy_with_se = param_list_all %>% 
  dplyr::select(species, species_nice, species_phy, Est.Error, parameter_model) %>% 
  pivot_wider(names_from = parameter_model, values_from = Est.Error, names_glue = "{parameter_model}_se") %>%
  left_join(param_list_for_phy) %>% 
  dplyr::select(-contains("Intercept"), -contains("sigma"), -contains("shape")) 

write_csv(param_list_for_phy_with_se, "rainfall-fitness/results/slopes_se.csv")

# Then see script phylo.signal.R




# ADDITIONAL FIGURES -----

## Slimmed data frame (one row = one sheet) ----

specimens_slim = read_csv("rainfall-fitness/data/all_data_compiled.csv")  %>% 
  filter(coordinateUncertaintyInMeters <= 5000) %>% 
  filter(entire_plant_confidence != "unsure") %>% 
  filter(!is.na(ppt_mm_rain_avg), 
         !is.na(tave_rain_avg), 
         !is.na(cwd_rain_avg), 
         !is.na(collection_DOWY), 
         !is.na(start_month_rain_avg)) %>%
  dplyr::select(genus, species, specimen, decimalLatitude, decimalLongitude, elev_m, collection_date, collection_DOWY, collection_year, start_month_rain_avg, first_rain_DOWY) %>% 
  distinct() %>% 
  mutate(species_nice = str_c(toupper(genus), ". ", species)) %>% 
  mutate(start_month_rain_WY = if_else(start_month_rain_avg >= 9 , start_month_rain_avg - 8, start_month_rain_avg + 4)) %>% 
mutate(days_since_germination = collection_DOWY-first_rain_DOWY) 



# Figure S1 ----------------

world = map_data("world")
states = map_data("state")

find_hull <- function(df) df[chull(df$decimalLatitude, df$decimalLongitude), ]
hulls <- plyr::ddply(specimens_slim, "species_nice", find_hull)

map = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey40", size = 1) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey40", size = 0.5) +
  coord_quickmap(xlim = c(-124, -114), ylim = c(33, 43)) +  
  geom_point(data = specimens_slim, aes(x = decimalLongitude, y = decimalLatitude, color = species_nice), alpha = 0.4) +
  geom_shape(data = hulls, aes(x = decimalLongitude, y = decimalLatitude, color = species_nice, fill = species_nice), alpha = 0.1, linewidth = 1, expand = 0.01, radius = 0.01) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic")) +
  scale_fill_discrete(name = "Species") +
  scale_color_discrete(name = "Species") +
  NULL; map

elevations  = ggplot(specimens_slim, aes(x = elev_m, color = species_nice, fill = species_nice)) +
  geom_density_ridges(aes(y = reorder(species_nice, desc(species_nice))), alpha = 0.3, show.legend = FALSE, bandwidth = 25, panel_scaling = TRUE) +
  scale_x_continuous(name = "Elevation (m)") +
  theme_ridges() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5), axis.text.y = element_text(face = "italic")); elevations

plot_grid(NULL, elevations, map, ncol = 3, rel_widths = c(0.05, 1, 1), labels = c("","a.", "b."), label_x= c(0, -0.05, 0))
ggsave("rainfall-fitness/figs/FigS1_map_elevations.pdf", width = 12, height = 5)
ggsave("rainfall-fitness/figs/FigS1_map_elevations.png", width = 12, height = 5)


sum(!is.na(specimens_slim$elev_m))






# Figure S2: collection years ----

summary(specimens_slim$collection_year)

coll_year = ggplot(specimens_slim, aes(x = collection_year, fill = species_nice)) +
  geom_histogram(binwidth = 10) +
  # geom_histogram(aes(y = reorder(species_nice, desc(species_nice))), alpha = 0.3, show.legend = FALSE, bandwidth = 2, panel_scaling = TRUE) +
  scale_x_continuous(name = "Collection year") +
  # theme_ridges() +
  labs(y = "Count of specimens") +
  facet_grid(species_nice~.) +
  guides(fill = "none") +
  theme(axis.title.x = element_text(hjust = 0.5), strip.text.y = element_text(angle = 0, face = "italic")); coll_year
ggsave("rainfall-fitness/figs/FigS2_specimens_years.pdf", height = 7, width = 4)
ggsave("rainfall-fitness/figs/FigS2_specimens_years.png", height = 7, width = 4)



# Figure S5 -------------------

specimens = read_csv("rainfall-fitness/data/all_data_compiled.csv") %>% 
  mutate(total_repro = bud + infBud + flwr + infFlwr + immFrt + fillFrt + infFrt + unfFrt + infImmFrt + unkRepStr,
         total_buds = bud + infBud,
         total_flowers = flwr + infFlwr,
         total_immFrts = immFrt + infImmFrt,
         total_matFrts = fillFrt + infFrt + unfFrt,
         weighted_fruits = round(0.33*total_flowers + 0.67*total_immFrts + 1*total_matFrts, 0),
         pheno_stage = weighted_fruits/total_repro,
         prop_fruits = (total_matFrts)/total_repro) %>% 
  dplyr::select(genus, species, specimen, decimalLatitude, decimalLongitude, elev_m, collection_date, collection_DOWY, pheno_stage, total_repro, weighted_fruits, prop_fruits) %>% 
  mutate(species_nice = str_c(toupper(genus), ". ", species)) 

plot(specimens$prop_fruits, specimens$pheno_stage)
cor(specimens$prop_fruits, specimens$pheno_stage)

phenostages =
  ggplot(dat1, aes(y = fct_rev(species_nice), x = pheno_stage, fill = species_nice)) +
  geom_boxplot(alpha = 0.7) +
  labs(y = "Species", x = "Phenological stage") +
  guides(fill = "none") +
  theme_cowplot() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5), axis.text.y = element_text(face = "italic")) ; phenostages

coll_date = ggplot(specimens_slim, aes(x = collection_DOWY, color = species_nice, fill = species_nice)) +
  geom_density_ridges(aes(y = reorder(species_nice, desc(species_nice))), alpha = 0.3, show.legend = FALSE, bandwidth = 2, panel_scaling = TRUE) +
  scale_x_continuous(name = "Collection day-of-year", breaks = c(154, 182, 213, 243, 274, 304, 335, 366, 396), labels = c("Feb 1", "Mar 1", "Apr 1", "May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1", "Oct 1")) +
  theme_ridges() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)); coll_date


phenostage_subpanel = plot_grid(NULL, phenostages, NULL, ncol = 1, rel_heights = c(0.025, 1, 0.095)) 
plot_grid(NULL, phenostage_subpanel, NULL, coll_date, ncol = 4, rel_widths = c(0.05,1,0.05,1),labels = c("a.","", "b.", ""))

ggsave("rainfall-fitness/figs/FigS5_collection_dates.pdf", width = 12, height = 5)
ggsave("rainfall-fitness/figs/FigS5_collection_dates.png", width = 12, height = 5)




# Figure S6: phenological stage and reproductive output vs. collection DOWY -----

param_listP = m1_all_p %>% 
  dplyr::select(species, variable, mean, lower, upper, species_nice, species_phy_ordered) %>% 
  filter(variable %in% c("b_collection_DOWY_scaled")) %>%
  mutate(variable_nice = fct_recode(variable,
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>%  
  mutate(response = "pheno", model = "model1")


plot_listP = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_listP[[i]] = mod_list1_p[[i]] %>%
    ggemmeans(., terms = c("collection_DOWY_scaled")) %>%
    # ggemmeans(., terms = c("collection_DOWY_scaled"), bias_correction = TRUE) %>%
    as_tibble() %>%
    mutate(species = species_choose)
  
}

pheno_DOY_plot = plot_listP %>% 
  bind_rows() %>% 
  mutate(collection_DOWY = x * sd(dat1$collection_DOWY) + mean(dat1$collection_DOWY)) %>% 
  left_join(., species_names) %>% 
  ggplot() +
  geom_point(data = dat1, aes(x = collection_DOWY, y = pheno_stage, color = species_nice), alpha = 0.2) +
  geom_line(aes(x = collection_DOWY, y = predicted, color = species_nice)) +
  geom_ribbon(aes(x = collection_DOWY, ymin = conf.low, ymax = conf.high, fill = species_nice), alpha = 0.3) +
  labs(y = "Phenological score", color = "Species", fill = "Species") +
  guides(fill = "none", color = "none") +
  scale_x_continuous(name = "Collection DOWY", breaks = c(154, 182, 213, 243, 274, 304, 335, 366, 396), labels = c("Feb 1", "Mar 1", "Apr 1", "May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1", "Oct 1")) +
  # facet_wrap(. ~ species_nice, nrow = 3) +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), legend.text = element_text(face = "italic"), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)); pheno_DOY_plot

param_listR = m1_all %>% 
  dplyr::select(species, variable, mean, lower, upper, species_nice, species_phy_ordered) %>% 
  filter(variable %in% c("b_collection_DOWY_scaled")) %>%
  mutate(variable_nice = fct_recode(variable,
                                    "Collection DOWY" = "b_collection_DOWY_scaled"),
         sig = if_else((upper > 0 & lower > 0), "sig_pos",
                       if_else((upper <0 & lower < 0), "sig_neg", "ns"))) %>% 
  mutate(response = "repro", model = "model1")

plot_listR = vector(mode = "list", length = length(species_list))

for (i in 1:length(species_list)) {
  
  species_choose = species_list[i]
  
  dat_slim = dat1 %>% 
    filter(species == species_choose) %>% 
    mutate(re_group = as.integer(as.factor(group)))
  
  plot_listR[[i]] = mod_list1[[i]] %>% 
    ggemmeans(., terms = c("collection_DOWY_scaled")) %>%
    as_tibble() %>%
    mutate(species = species_choose)
  
  print(species_choose)
}

repro_DOY_plot = plot_listR %>% 
  bind_rows() %>% 
  mutate(collection_DOWY = x * sd(dat1$collection_DOWY) + mean(dat1$collection_DOWY)) %>% 
  left_join(., species_names) %>% 
  ggplot() +
  geom_point(data = dat1, aes(x = collection_DOWY, y = total_repro, color = species_nice), alpha = 0.2) +
  geom_line(aes(x = collection_DOWY, y = predicted, color = species_nice)) +
  geom_ribbon(aes(x = collection_DOWY, ymin = conf.low, ymax = conf.high, fill = species_nice), alpha = 0.3) +
  labs(y = "Total reproductive structures", color = "Species", fill = "Species") +
  scale_x_continuous(name = "Collection DOWY", breaks = c(154, 182, 213, 243, 274, 304, 335, 366, 396), labels = c("Feb 1", "Mar 1", "Apr 1", "May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1", "Oct 1")) +
  # facet_wrap(. ~ species_nice, nrow = 3) +
  scale_y_log10() +
  theme_cowplot() +
  theme(axis.title.x = element_text(hjust = 0.5), legend.text = element_text(face = "italic"), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)); repro_DOY_plot

plot_grid(pheno_DOY_plot, repro_DOY_plot, labels = c("a.", "b."), rel_widths = c(1, 1.3))
ggsave("rainfall-fitness/figs/FigS6_repro_pheno_vs_DOY.pdf", width = 10, height = 4)
ggsave("rainfall-fitness/figs/FigS6_repro_pheno_vs_DOY.png", width = 10, height = 4)




# Figure S7 ----

estsR = ggplot(data = param_listR) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  scale_color_manual(values = c("grey60", "dodgerblue")) +
  # scale_x_continuous(limits = c(-0.85, 0.85)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect of collection\nDOWY on reproduction") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()); estsR

estsP = ggplot(data = param_listP) +
  geom_point(aes(x = mean, y = species_phy_ordered, color = sig), size = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = species_phy_ordered, color = sig), height = 0, size = 0.7) +
  scale_color_manual(values = c("grey60", "dodgerblue")) +
  # scale_x_continuous(limits = c(-0.85, 0.85)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = "none") +
  labs(y = "Species", x = "Estimated effect of collection\nDOWY on phenology") +
  theme(axis.text.y = element_text(face = "italic", hjust = 1), axis.title.y = element_blank()); estsP

plot_grid(estsP, estsR, rel_widths = c(1, 0.4), rel_heights = c(1, 1), nrow = 1, labels = c("a.", "b."), label_x = c(0.25, -0.06))

ggsave("rainfall-fitness/figs/FigS7_repro_pheno_CDOY_pointranges.pdf", height = 6, width = 7.5)
ggsave("rainfall-fitness/figs/FigS7_repro_pheno_CDOY_pointranges.png", height = 6, width = 7.5)

