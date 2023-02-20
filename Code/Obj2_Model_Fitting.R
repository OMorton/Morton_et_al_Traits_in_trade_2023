
#######################
#### Model Fitting ####
#######################

## Purpose
## Script to do the final processing of the CITES data and fit the model.
## The model is written to run on the Sharc HPC, running on a standard computer is not feasible.
## Setting up the cmdstan toolchain is required for the within chain threading to make the model tractable.

## global options settings
options(scipen = 999) ## remove scientific notation, namely "e"  to show very small and large numbers in full.
options(na.action = "na.pass")

library(brms)
library(tidyverse)
library(ape)
library(picante)
library(cmdstanr)

## 758 sp to start
LH_Final <- read.csv("Data/Life_History/CITES_Life_History/CITES_IUCN_LH_Match_2_22_NoOcc.csv") %>% select(-X)
n_distinct(LH_Final$Name_for_CITESdb) ## 795

#### Tidying ####
## Standardise the predictors and tidy data.
## We keep series where all entries are 0, for example a species was traded only captive sourced therefore we keep is wild sourced series as 0's
## This is to reflect instances like above where that species has a probs approaching 0 of being in wild sourced trade.
## The lack of any >0 values for some series could cause issues - when we add the species effects this might help.
SLH_Final <- LH_Final %>%
  filter(Name_for_LH != "Rhodonessa caryophyllacea") %>%## Remove - missing bodymass also poss extinct
  ## Round up traded volumes to whole digits, 1.5 would be 2 as you cant trade .5 of an animal.
  ## Convert source and threat to a composite category
  mutate(n = ceiling(n),
         WildSourceB = if_else(WildSource == "Yes", 1, 0),
         WildSourceB = as.factor(WildSourceB),
         WildThreat = case_when((WildSource == "Yes" & Threat_code == "Non-threatened")~ "WNT",
                                (WildSource == "Yes" & Threat_code == "Threatened")~ "WT",
                                (WildSource == "No" & Threat_code == "Non-threatened")~ "CNT",
                                (WildSource == "No" & Threat_code == "Threatened")~ "CT",
                                (WildSource == "Yes" & Threat_code == "Not assessed")~ "WNA",
                                (WildSource == "No" & Threat_code == "Not assessed")~ "CNA"),
         ## Subsp to sp (apart from grey and timneh parrots - timneh are recognized species and have distinct
         ## reproductive trait values)
         Name_for_CITESdb = case_when(Name_for_CITESdb == "Amazona festiva festiva" ~ "Amazona festiva",
                                      Name_for_CITESdb == "Branta canadensis leucopareia" ~ "Branta canadensis",
                                      Name_for_CITESdb == "Buceros hydrocorax hydrocorax" ~ "Buceros hydrocorax",
                                      Name_for_CITESdb == "Falco pelegrinoides babylonicus" ~ "Falco pelegrinoides",
                                      Name_for_CITESdb == "Falco peregrinus anatum" ~ "Falco peregrinus",
                                      Name_for_CITESdb == "Falco peregrinus peregrinus" ~ "Falco peregrinus",
                                      Name_for_CITESdb == "Phoenicopterus ruber ruber" ~ "Phoenicopterus ruber",
                                      Name_for_CITESdb == "Pterocnemia pennata pennata" ~ "Pterocnemia pennata",
                                      Name_for_CITESdb == "Alisterus chloropterus mozskowskii" ~ "Alisterus chloropterus",
                                      Name_for_CITESdb == "Bubo bubo bengalensis" ~ "Bubo bubo",
                                      Name_for_CITESdb == "Eos squamata riciniata" ~ "Eos squamata",
                                      Name_for_CITESdb == "Grus canadensis pratensis" ~ "Grus canadensis",
                                      Name_for_CITESdb == "Poephila cincta cincta" ~ "Poephila cincta",
                                      Name_for_CITESdb == "Rhea americana albescens" ~ "Rhea americana",
                                      TRUE ~ (as.character(Name_for_CITESdb)))) %>%
  ## New totals at the species level
  group_by(across(c(-n))) %>%
  tally(n) %>%
  ungroup() %>%
  arrange(WildThreat, Name_for_CITESdb, Year) %>%
  ## Create the lagged volume term and remove the trades before 2000
  mutate(lagN = lag(n))  %>%
  filter(Year>-1) %>%
  group_by(Name_for_CITESdb) %>% filter(sum(n) > 0) %>% ## remove species with all 0's
  ungroup() %>%
  ## Standardise everything to the trade dataset level (not the full avian trait value data).
  mutate(SLag = (lagN - mean(lagN))/sd(lagN), 
         SMax_longevity = (log2(Max_longevity) - mean(log2(Max_longevity)))/sd(log2(Max_longevity)), 
         SLogBodymass = (log2(Bodymass) - mean(log2(Bodymass)))/sd(log2(Bodymass)),
         SFirst_repro = (log2(Age_at_first_breeding) - mean(log2(Age_at_first_breeding)))/sd(log2(Age_at_first_breeding)),
         LogMax_longevity = log2(Max_longevity), 
         LogBodymass = log2(Bodymass),
         LogFirst_repro = log2(Age_at_first_breeding),
         SYear = (Year - mean(Year))/sd(Year)) %>%
  group_by(Name_for_CITESdb, WildSource) %>%
  mutate(group_id = as.factor(cur_group_id()),
         Year = as.factor(Year))

length(unique(SLH_Final$Name_for_CITESdb)) ## 760
length(unique(SLH_Final$group_id)) ## 1520 twice the number of species as each species has two potential series (wild/captive)

## Checking for NA's
SLH_Final %>% filter(is.na(SLogBodymass))
SLH_Final %>% filter(is.na(SLag))
SLH_Final %>% filter(is.na(SYear))

#### HPC Model ####
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
SLH_Final <- SLH_Final %>%
  mutate(Name_for_phylo = str_replace(Name_for_CITESdb, " ", "_"))

## Correct the phylonames
SLH_Final_phy <- SLH_Final %>% 
  mutate(Name_for_phylo = case_when(Name_for_phylo == "Alisterus_chloropterus mozskowskii" ~ "Alisterus_chloropterus",
                                    Name_for_phylo == "Amazona_festiva festiva" ~ "Amazona_festiva",
                                    Name_for_phylo == "Amazona_xanthops" ~ "Alipiopsitta_xanthops",
                                    Name_for_phylo == "Baillonius_bailloni" ~ "Pteroglossus_bailloni",
                                    Name_for_phylo == "Barnardius_barnardi" ~ "Barnardius_zonarius",
                                    Name_for_phylo == "Branta_canadensis leucopareia" ~ "Branta_canadensis",
                                    Name_for_phylo == "Bubo_bubo bengalensis" ~ "Bubo_bubo",
                                    Name_for_phylo == "Buceros_hydrocorax hydrocorax" ~ "Buceros_hydrocorax",
                                    Name_for_phylo == "Buteo_poecilochrous" ~ "Buteo_polyosoma",
                                    Name_for_phylo == "Eos_rubra" ~ "Eos_bornea",
                                    Name_for_phylo == "Eos_squamata riciniata" ~ "Eos_squamata",
                                    Name_for_phylo == "Falco_pelegrinoides babylonicus" ~ "Falco_pelegrinoides",
                                    Name_for_phylo == "Falco_peregrinus anatum" ~ "Falco_peregrinus",
                                    Name_for_phylo == "Falco_peregrinus peregrinus" ~ "Falco_peregrinus",
                                    Name_for_phylo == "Grus_canadensis pratensis" ~ "Grus_canadensis",
                                    Name_for_phylo == "Guarouba_guarouba" ~ "Guaruba_guarouba",
                                    Name_for_phylo == "Hieraaetus_fasciatus" ~ "Aquila_fasciatus",
                                    Name_for_phylo == "Nigrita_canicapilla" ~ "Nigrita_canicapillus",
                                    Name_for_phylo == "Nyctea_scandiaca" ~ "Bubo_scandiaca",
                                    Name_for_phylo == "Otus_choliba" ~ "Megascops_choliba",
                                    Name_for_phylo == "Otus_kennicottii" ~ "Megascops_kennicottii",
                                    Name_for_phylo == "Otus_roboratus" ~ "Megascops_roboratus",
                                    Name_for_phylo == "Otus_watsonii" ~ "Megascops_watsonii",
                                    Name_for_phylo == "Phoenicopterus_ruber ruber" ~ "Phoenicopterus_ruber",
                                    Name_for_phylo == "Pionopsitta_barrabandi" ~ "Pyrilia_barrabandi",
                                    Name_for_phylo == "Platycercus_barnardi" ~ "Barnardius_zonarius",
                                    Name_for_phylo == "Platycercus_zonarius" ~ "Barnardius_zonarius",
                                    Name_for_phylo == "Poephila_cincta cincta" ~ "Poephila_cincta",
                                    Name_for_phylo == "Psittacula_calthorpae" ~ "Psittacula_calthropae",
                                    Name_for_phylo == "Psittacula_echo" ~ "Psittacula_eques",
                                    Name_for_phylo == "Psittacus_erithacus timneh" ~ "Psittacus_erithacus",
                                    Name_for_phylo == "Pterocnemia_pennata" ~ "Rhea_pennata",
                                    Name_for_phylo == "Pterocnemia_pennata pennata" ~ "Rhea_pennata",
                                    Name_for_phylo == "Rhea_americana albescens" ~ "Rhea_americana",
                                    Name_for_phylo == "Spizaetus_africanus" ~ "Aquila_africanus",
                                    Name_for_phylo == "Spizaetus_nipalensis" ~ "Nisaetus_nipalensis",
                                    Name_for_phylo == "Spizastur_melanoleucus" ~ "Spizaetus_melanoleucus",
                                    Name_for_phylo == "Streptopelia_senegalensis" ~ "Stigmatopelia_senegalensis",
                                    Name_for_phylo == "Torgos_tracheliotus" ~ "Torgos_tracheliotos",
                                    Name_for_phylo == "Afrotis_afraoides" ~ "Eupodotis_afraoides",
                                    Name_for_phylo == "Anthropoides_paradiseus" ~ "Grus_paradisea",
                                    Name_for_phylo == "Anthropoides_virgo" ~ "Grus_virgo",
                                    Name_for_phylo == "Aratinga_maculata" ~ "Aratinga_solstitialis",
                                    Name_for_phylo == "Asarcornis_scutulata" ~ "Cairina_scutulata",
                                    Name_for_phylo == "Berenicornis_comatus" ~ "Aceros_comatus",
                                    Name_for_phylo == "Bugeranus_carunculatus" ~ "Grus_carunculatus",
                                    Name_for_phylo == "Catreus_wallichii" ~ "Catreus_wallichi",
                                    Name_for_phylo == "Chlamydotis_macqueenii" ~ "Chlamydotis_undulata",
                                    Name_for_phylo == "Ciccaba_huhula" ~ "Strix_huhula",
                                    Name_for_phylo == "Eolophus_roseicapilla" ~ "Cacatua_roseicapilla",
                                    Name_for_phylo == "Hieraaetus_wahlbergi" ~ "Aquila_wahlbergi",
                                    Name_for_phylo == "Lissotis_melanogaster" ~ "Eupodotis_melanogaster",
                                    Name_for_phylo == "Lonchura_oryzivora" ~ "Padda_oryzivora",
                                    Name_for_phylo == "Lophotis_ruficrista" ~ "Eupodotis_ruficrista",
                                    Name_for_phylo == "Micronisus_gabar" ~ "Melierax_gabar",
                                    Name_for_phylo == "Otus_asio" ~ "Megascops_asio",
                                    Name_for_phylo == "Poicephalus_fuscicollis" ~ "Poicephalus_robustus",
                                    Name_for_phylo == "Ptilopsis_leucotis" ~ "Otus_leucotis",
                                    Name_for_phylo == "Pyrrhura_caeruleiceps" ~ "Pyrrhura_picta",
                                    Name_for_phylo == "Rhyticeros_plicatus" ~ "Aceros_plicatus",
                                    Name_for_phylo == "Rhyticeros_undulatus" ~ "Aceros_undulatus",
                                    Name_for_phylo == "Tauraco_schuettii" ~ "Tauraco_schuetti",
                                    Name_for_phylo == "Trichoglossus_capistratus" ~ "Trichoglossus_haematodus",
                                    Name_for_phylo == "Trichoglossus_forsteni" ~ "Trichoglossus_haematodus",
                                    Name_for_phylo == "Spizaetus_cirrhatus" ~ "Nisaetus_cirrhatus",
                                    TRUE ~ (as.character(Name_for_phylo))))

length(unique(SLH_Final_phy$Name_for_phylo)) ## 751
Phylo_names <- SLH_Final_phy %>% group_by(Name_for_phylo) %>% slice_head() %>%
  column_to_rownames(var = 'Name_for_phylo')

## pruned to now have the 758 species
phy_bird <- match.phylo.data(MCC_tree, Phylo_names)$phy

rownames(Phylo_names) %>% as.data.frame() %>% filter(!. %in% phy_bird$tip.label)
length(unique(SLH_Final_phy$Name_for_CITESdb)) ## 760
SLH_Final_phy %>% filter(n>0) %>% group_by(WildSource) %>% summarise(n_distinct(Name_for_CITESdb))

VCV_mat <- vcv(phy_bird)

write.csv(SLH_Final_phy, "../1_Life_History/Outputs/SOM/CITES_Fitting_Data.csv")

## For traded species as per the Scheffers data.
All1 <- brm(bf(n ~ WildSource*Threat_code*SYear +
                 WildSource*SLogBodymass +
                 WildSource*SFirst_repro +
                 WildSource*SMax_longevity +
                 (1+ WildSource | Year) +
                 (1 + SYear + WildSource + SYear:WildSource | Name_for_CITESdb) +
                 (1|gr(Name_for_phylo, cov = A)),
               hu ~ WildSource*Threat_code*SYear +
                 WildSource*SLogBodymass +
                 WildSource*SFirst_repro +
                 WildSource*SMax_longevity +
                 (1+ WildSource | Year) +
                 (1 + SYear + WildSource + SYear:WildSource | Name_for_CITESdb) +
                 (1|gr(Name_for_phylo, cov = A))),
            family = hurdle_negbinomial(), 
            sample_prior = TRUE,
            prior = c(
              prior(normal(0,1), "b"),
              prior(normal(0,1), "Intercept"),
              prior(normal(0,1), "b", dpar = "hu"),
              prior(normal(0,1), "Intercept", dpar = "hu"),
              prior(normal(0,1), "sd"),
              prior(normal(0,1), "sd", dpar = "hu"),
              prior(lkj(2), "cor")),
            data = SLH_Final_phy,
            data2 = list(A = VCV_mat),
            backend = "cmdstanr", threads = threading(2),
            chains = 4, iter = 3000, thin = 1, cores = 4, warmup = 1000)


library(bayestestR)
hu_sum <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% select(contains("hu_")) %>% select(!hu_Intercept)
hu_sum <- hu_sum*-1
mu_sum <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% select(!contains("hu_")) 

hu_plot <- plot(p_direction(hu_sum)) + 
  theme_bw(base_size = 12) + 
  labs(title = element_blank()) +
  theme(legend.position = "none", panel.grid = element_blank())
mu_plot <- plot(p_direction(mu_sum)) + 
  theme_bw(base_size = 12) + 
  labs(title = element_blank()) +
  theme(legend.position = "none", panel.grid = element_blank())

library(cowplot)
Sum_plot <- plot_grid(hu_plot, mu_plot,labels = c("A.", "B."), nrow = 2, align = "v")
ggsave("../1_Life_History/Outputs/SOM/ModSum.png", Sum_plot, device = "jpeg",
       units = "cm", width = 30, height = 30)
