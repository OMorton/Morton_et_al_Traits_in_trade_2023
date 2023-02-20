####################################
##---Traded, Listed and traits ---##
####################################


#### Packages and presets ####
options(scipen = 999) ## remove scientific notation, namely "e"  to show very small and large numbers in full.
options(na.action = "na.pass")
library(brms)
library(tidyverse)
library(tidybayes)
library(bayestestR)

#### Data sources ####

## data for 10254 bird species
Trade_dat <- read.csv("Data/CITES/Traded_Traits_Data") %>% 
  rename("phyloname" = "phyloname.x")

#### Exploration ####
Trade_dat <- Trade_dat %>% mutate(Listing_bin = ifelse(Listing == "Listed", 1, 0))
Trade_dat %>% summarise(Sp2000 = sum(Trade_since_2000), Sp2015 = sum(Trade_since_2015), SpList = sum(Listing_bin))

## Standardize the continuous trait data.
## For bodymass we standardize log bodymass.
unique(Trade_dat$pet_prod)


Trade_dat <- Trade_dat %>% 
  mutate(Age_at_first_breeding_z = (log2(Age_at_first_breeding) - mean(log2(Age_at_first_breeding)))/
           sd(log2(Age_at_first_breeding)),
         Bodymass_z = (Bodymass - mean(Bodymass))/sd(Bodymass),
         Bodymass_logz = (log2(Bodymass) - mean(log2(Bodymass)))/sd(log2(Bodymass)),
         Max_longevity_z = (log2(Max_longevity) - mean(log2(Max_longevity)))/sd(log2(Max_longevity)))


#### Multilevel Model ####

library(ape)
library(picante)
## Read in 1000 avian bird trees.
## Available from https://data.vertlife.org/
#Bird_trees <- read.tree("Data/Vertlife_supertrees/AllBirdsHackett1.tre")

## Get the MCC tree. This takes some time. 
#MCC_tree <- maxCladeCred(Bird_trees)

## To save time write out the MCC tree
#write.nexus(MCC_tree, file = "Data/Vertlife_supertrees/AllBirdsMCC.nex")
MCC_tree <- read.nexus("Data/Vertlife_supertrees/AllBirdsMCC.nex")

## Prune the tree to our list of species.
## make names in phylogeny row names and add underscores.
Phylo_names <- Trade_dat %>% group_by(phyloname) %>% slice_head() %>%
  mutate(phyloname = str_replace(phyloname, " ", "_")) %>%
  column_to_rownames(var = 'phyloname')

Trade_dat <- Trade_dat %>%
  mutate(phyloname = str_replace(phyloname, " ", "_"))

## pruned to now have the 9839 species
phy_bird <- match.phylo.data(MCC_tree, Phylo_names)$phy

VCV_mat <- vcv(phy_bird)

sum(Trade_dat$Trade_since_2000)
sum(Trade_dat$Listing_bin)

## For traded species as per the Scheffers data.
TR1_sp <- brm(Traded ~ Age_at_first_breeding_z + Bodymass_logz + Max_longevity_z + (1|gr(phyloname, cov = A)),
              family = bernoulli(),
              sample_prior = TRUE,
              prior = c(
                prior(normal(0,1), "b"),
                prior(normal(0,1), "Intercept"),
                prior(normal(0,.5), "sd")),
              #control = list(adapt_delta = 0.9),
              data = Trade_dat,
              data2 = list(A = VCV_mat),
              backend = "cmdstanr",
              chains = 4, iter = 1000, thin = 1, cores = 4, warmup = 500)


## For listed species
LI1_sp <- brm(Listing_bin ~ Age_at_first_breeding_z + Bodymass_logz + Max_longevity_z + (1|gr(phyloname, cov = A)),
              family = bernoulli(),
              sample_prior = TRUE,
              prior = c(
                prior(normal(0,1), "b"),
                prior(normal(0,1), "Intercept"),
                prior(normal(0,.5), "sd")),
              #control = list(adapt_delta = 0.9),
              data = Trade_dat,
              data2 = list(A = VCV_mat),
              backend = "cmdstanr",
              chains = 4, iter = 1000, thin = 1, cores = 4, warmup = 500)

## For CITES traded species
CITES20001_sp <- brm(Trade_since_2000 ~ Age_at_first_breeding_z + Bodymass_logz + Max_longevity_z + (1|gr(phyloname, cov = A)),
                     family = bernoulli(),
                     sample_prior = TRUE,
                     prior = c(
                       prior(normal(0,1), "b"),
                       prior(normal(0,1), "Intercept"),
                       prior(normal(0,.5), "sd")),
                     #control = list(adapt_delta = 0.9),
                     data = Trade_dat,
                     data2 = list(A = VCV_mat),
                     backend = "cmdstanr",
                     chains = 4, iter = 1000, thin = 1, cores = 4, warmup = 500)

TR_PD <- p_direction(TR1_sp) %>% as.data.frame()
LI_PD <- p_direction(LI1_sp) %>% as.data.frame()
CT_PD <- p_direction(CITES20001_sp) %>% as.data.frame()

#### Plotting ####
## Create simulated data for examing each trend conditionally
Age_dat <- Trade_dat %>% summarise(Age_at_first_breeding_z = seq(from = min(Age_at_first_breeding_z), to = max(Age_at_first_breeding_z), length.out = 100),
                                   Max_longevity_z = 0, Bodymass_logz = 0)
Long_dat <- Trade_dat %>% summarise(Max_longevity_z = seq(from = min(Max_longevity_z), to = max(Max_longevity_z), length.out = 100),
                                    Age_at_first_breeding_z = 0, Bodymass_logz = 0)

Body_dat <- Trade_dat %>% summarise(Bodymass_logz = seq(from = min(Bodymass_logz), to = max(Bodymass_logz), length.out = 100),
                                    Max_longevity_z = 0, Age_at_first_breeding_z = 0)



## Extract fitted posterior population level estimates
Age_fit_TR1 <-  add_epred_draws(TR1_sp, newdata =  Age_dat, re_formula = NA) %>% group_by(Age_at_first_breeding_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Traded")
Long_fit_TR1 <-  add_epred_draws(TR1_sp, newdata =  Long_dat, re_formula = NA) %>% group_by(Max_longevity_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Traded")
Body_fit_TR1 <-  add_epred_draws(TR1_sp, newdata =  Body_dat, re_formula = NA) %>% group_by(Bodymass_logz) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Traded")

Age_fit_LI1 <-  add_epred_draws(LI1_sp, newdata =  Age_dat, re_formula = NA) %>% group_by(Age_at_first_breeding_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Listed")
Long_fit_LI1 <-  add_epred_draws(LI1_sp, newdata =  Long_dat, re_formula = NA) %>% group_by(Max_longevity_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Listed")
Body_fit_LI1 <-  add_epred_draws(LI1_sp, newdata =  Body_dat, re_formula = NA) %>% group_by(Bodymass_logz) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "Listed")

Age_fit_CT20001 <-  add_epred_draws(CITES20001_sp, newdata =  Age_dat, re_formula = NA) %>% group_by(Age_at_first_breeding_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "CITES2000")
Long_fit_CT20001 <-  add_epred_draws(CITES20001_sp, newdata =  Long_dat, re_formula = NA) %>% group_by(Max_longevity_z) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "CITES2000")
Body_fit_CT20001 <-  add_epred_draws(CITES20001_sp, newdata =  Body_dat, re_formula = NA) %>% group_by(Bodymass_logz) %>% median_hdci(.epred, .width = .9) %>% mutate(Model = "CITES2000")

## Bind the results of all models
Age_all_models <- rbind(Age_fit_TR1, Age_fit_LI1, Age_fit_CT20001)
Long_all_models <- rbind(Long_fit_TR1, Long_fit_LI1, Long_fit_CT20001)
Body_all_models <- rbind(Body_fit_TR1, Body_fit_LI1, Body_fit_CT20001)



## Age at first breeding

Age_all_models <- Age_all_models %>% mutate(Age_at_first_breeding = 
                                              2^(Age_at_first_breeding_z*sd(log2(Trade_dat$Age_at_first_breeding)) + 
                                                   mean(log2(Trade_dat$Age_at_first_breeding))))

## Seperate the models to plot
TR_age <- Age_all_models %>% filter(Model == "Traded")
LI_age <- Age_all_models %>% filter(Model == "Listed")
CT2000_age <- Age_all_models %>% filter(Model == "CITES2000")

## All effects but trade neither certain positive or negative.
(Age_fit_plot <- ggplot(, aes(Age_at_first_breeding, .epred, colour = Model, group = Model)) + 
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .25, linetype = "dashed", size = 1) +
  scale_colour_manual(values = c("black", "darkblue", "skyblue3") ,limits=c("Traded", "Listed", "CITES2000")) +
  geom_ribbon(data = TR_age, aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA, fill = "black") +
  geom_line(data = TR_age, size = 1.5, linetype = "solid") +
  geom_line(data = LI_age, size = 1.5, linetype = "dashed") +
  geom_line(data = CT2000_age, size = 1.5, linetype = "dashed") +
  ylab("Probability of trade/listing") +
  xlab("First reproduction (Years)") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), legend.position = "none"))


## Bodymass
##  plot together as all certain direction
Body_all_models <- Body_all_models %>% mutate(Bodymass_log = 
                                                2^(Bodymass_logz*sd(log2(Trade_dat$Bodymass)) + 
                                                     mean(log2(Trade_dat$Bodymass))))

(Body_fit_plot <- ggplot(Body_all_models, aes(Bodymass_log, .epred, colour = Model, fill = Model, group = Model)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA) +
  geom_line(size = 1.5) +
  scale_colour_manual(values = c("black", "darkblue", "skyblue3") ,limits=c("Traded", "Listed", "CITES2000")) +
  scale_fill_manual(values = c("black", "darkblue", "skyblue3") ,limits=c("Traded", "Listed", "CITES2000")) +
  ylab("Probability of trade/listing") +
  xlab("Body mass (g)") +
  coord_cartesian(ylim = c(0,1), xlim = c(1, 111001)) +
  scale_x_log10(breaks = c(1, 100, 10000)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), legend.position = "none"))


## seperate to plot
Long_all_models <- Long_all_models %>% mutate(Max_longevity = 
                                                2^(Max_longevity_z*sd(log2(Trade_dat$Max_longevity)) + 
                                                     mean(log2(Trade_dat$Max_longevity))))

TR_long <- Long_all_models %>% filter(Model == "Traded")
LI_long <- Long_all_models %>% filter(Model == "Listed")
CT2000_long <- Long_all_models %>% filter(Model == "CITES2000")

(Long_fit_plot <- ggplot(, aes(Max_longevity, .epred, colour = Model, group = Model)) + 
  geom_ribbon(data = TR_long, aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA, fill = "black") +
  #geom_ribbon(data = CT2000_long, aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA, fill = "skyblue3") +
  geom_line(data = TR_long, size = 1.5) +
  geom_line(data = LI_long, size = 1.5, linetype = "dashed") +
  geom_line(data = CT2000_long, size = 1.5, linetype = "dashed") +
  scale_colour_manual(values = c("black", "darkblue", "skyblue3") ,limits=c("Traded", "Listed", "CITES2000")) +
  ylab("Probability of trade/listing") +
  xlab("Maximum longevity (Years)") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), legend.position = "none"))

#### Fixef ####
library(tidybayes)

TR_fixef <- fixef(TR1_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par") %>%
  group_by(Par) %>% median_hdci(Est, .width = .9)
TR_fixef_raw <- fixef(TR1_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par")

LI_fixef <- fixef(LI1_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par") %>%
  group_by(Par) %>% median_hdci(Est, .width = .9)
LI_fixef_raw <- fixef(LI1_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par")

CT_fixef <- fixef(CITES20001_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par") %>%
  group_by(Par) %>% median_hdci(Est, .width = .9)
CT_fixef_raw <- fixef(CITES20001_sp, summary = FALSE) %>% as.data.frame() %>% select(-Intercept) %>% 
  pivot_longer(everything(), values_to = "Est", names_to = "Par")

library(ggridges)
## Plot
(TR_fixef_plt <- ggplot(TR_fixef, aes(Est, Par)) +
  geom_density_ridges(data = TR_fixef_raw, aes(Est), alpha=0.5, 
                      stat="binline", bins=60, draw_baseline = FALSE, scale = .9, fill = "black", colour = NA) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0, size = .75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("First reproduction", "Body mass", "Longevity")) +
  xlab("Possible parameter values") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), panel.grid = element_blank()))

(LI_fixef_plt <- ggplot(LI_fixef, aes(Est, Par)) +
  geom_density_ridges(data = LI_fixef_raw, aes(Est), alpha=0.5, 
                      stat="binline", bins=60, draw_baseline = FALSE, scale = .9, fill = "darkblue", colour = NA) +
  geom_point(size = 2, colour = "darkblue") +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0, colour = "darkblue", size = .75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("First reproduction", "Body mass", "Longevity")) +
  xlab("Possible parameter values") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), panel.grid = element_blank()))

(CT_fixef_plt <- ggplot(CT_fixef, aes(Est, Par)) +
  geom_density_ridges(data = CT_fixef_raw, aes(Est), alpha=0.5, 
                      stat="binline", bins=60, draw_baseline = FALSE, scale = .9, fill = "skyblue3", colour = NA) +
  geom_point(size = 2, colour = "skyblue3") +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0, colour = "skyblue3", size = .75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c( "First reproduction", "Body mass", "Longevity")) +
  xlab("Possible parameter values") +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank(), panel.grid = element_blank()))

#### Arrangement ####

library(ggpubr)

## Legend
plt_leg <- ggplot(, aes(Max_longevity, .epred, colour = Model, group = Model)) + 
  geom_ribbon(data = TR_long, aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA, fill = "black") +
  geom_ribbon(data = CT2000_long, aes(ymin = .lower, ymax = .upper), alpha = .35, colour = NA, fill = "skyblue3") +
  geom_line(data = TR_long, size = 1.5) +
  geom_line(data = LI_long, size = 1.5, linetype = "dashed") +
  geom_line(data = CT2000_long, size = 1.5) +
  scale_colour_manual(values = c("black", "darkblue", "skyblue3") ,limits=c("Traded", "Listed", "CITES2000")) +
  ylab("Probability of trade/listing") +
  xlab("Maximum longevity (Years)") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.title = element_blank())

Legend <- get_legend(plt_leg)

empty <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "white", colour = NA))
Coef_plot <- ggarrange(empty, empty, empty, 
          TR_fixef_plt, LI_fixef_plt, CT_fixef_plt,
          heights = c(0.1,1), hjust = -.1,
          ncol = 3, nrow = 2, labels = 
            c("A. Probability of being broadly traded", 
              "B. Probability of being CITES Listed", 
              "C. Probability of being CITES traded"))

Cond_lines <- ggarrange(Body_fit_plot, Age_fit_plot, Long_fit_plot, nrow = 1, labels = c("D.", "E.", "F."))

Fig2 <- ggarrange(Coef_plot, empty, Cond_lines, legend.grob = Legend, 
                  legend = "bottom", nrow = 3, heights = c(1, .05, .75))

#### Outputs ####

TR_PD_Sum <- left_join(TR_fixef, mutate(TR_PD, Parameter = gsub("b_", "", Parameter)), by = c("Par" = "Parameter"))
LI_PD_Sum <- left_join(LI_fixef, mutate(LI_PD, Parameter = gsub("b_", "", Parameter)), by = c("Par" = "Parameter"))
CT_PD_Sum <- left_join(CT_fixef, mutate(CT_PD, Parameter = gsub("b_", "", Parameter)), by = c("Par" = "Parameter"))

write.csv(TR_PD_Sum, "D:/CITES_Trade/1_Life_History/Outputs/Figure2/TR_PD_Sum.csv")
write.csv(LI_PD_Sum, "D:/CITES_Trade/1_Life_History/Outputs/Figure2/LI_PD_Sum.csv")
write.csv(CT_PD_Sum, "D:/CITES_Trade/1_Life_History/Outputs/Figure2/CT_PD_Sum.csv")

TR_fixef_raw %>% group_by(Par) %>% mutate(OR = exp(Est)) %>% median_hdci(OR, .width = .9)
LI_fixef_raw %>% group_by(Par) %>% mutate(OR = exp(Est)) %>% median_hdci(OR, .width = .9)
CT_fixef_raw %>% group_by(Par) %>% mutate(OR = exp(Est)) %>% median_hdci(OR, .width = .9)
