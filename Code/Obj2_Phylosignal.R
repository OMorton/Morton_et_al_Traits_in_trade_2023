#########################
## Phylogenetic signal ##
#########################

## global options settings
options(scipen = 999) ## remove scientific notation, namely "e"  to show very small and large numbers in full.
options(na.action = "na.pass")

library(ggtree)
library(ggridges)
library(ggnewscale)

## data for 10254 bird species
Trade_dat <- read.csv("Data/CITES/Traded_Traits_Data") %>% 
  rename("phyloname" = "phyloname.x")
## CITES data fro traded species
SLH_Final_phy <- read.csv("../1_Life_History/Outputs/SOM/CITES_Fitting_Data.csv")


#### Listed ####
## Extract all listed species
Listed <- Trade_dat %>% filter(Listing == "Listed") ##1473
Listed_Phy <- Listed %>% distinct(phyloname, Order.x) ##1463

Sum <- SLH_Final_phy %>% group_by(Name_for_phylo, WildSource) %>% 
  summarise(total = sum(n)) %>%
  mutate(Presence = ifelse(total >0, 1, 0)) %>%
  select(-total) %>%
  pivot_wider(names_from = "WildSource", values_from = "Presence") %>%
  rename("Captive" = "No", "Wild" = "Yes")

colSums(Sum[,2:3]) ## 628 CB and 450 W phy species

Listed_tr <-left_join(Listed_Phy, Sum, by = c("phyloname" = "Name_for_phylo"))
colSums(Listed_tr[,3:4], na.rm = TRUE)

## test the missing species - all cases of was listed and currently is not (so in the current listings data)
test <-Listed_tr %>% filter(!is.na(Wild))
Sp_no_longer_listed <- Sum %>% filter(!Name_for_phylo %in% test$phyloname)

## Create final data
Listed_final <- Listed_tr %>% mutate(Captive = ifelse(is.na(Captive), 0, Captive),
                                     Wild = ifelse(is.na(Wild), 0, Wild))

## Make a match the phylo tree
MCC_tree <- read.nexus("Data/Vertlife_supertrees/AllBirdsMCC.nex")
phy_bird2 <- match.phylo.data(MCC_tree, column_to_rownames(Listed_final, "phyloname"))$phy

#### Phylogenetic signal ####
## Make the comp df
signal_df <- 
  caper::comparative.data(phy = phy_bird2, 
                          data = Listed_final, 
                          names.col = phyloname, 
                          na.omit = FALSE)

## Calc signal
capt_signal <- caper::phylo.d(signal_df, binvar = Captive)
plot(capt_signal)

wild_signal <- caper::phylo.d(signal_df, binvar = Wild)
plot(wild_signal)

## Pull out data for plotting - Capt
brownian <- capt_signal$Permutations$brownian
random   <- capt_signal$Permutations$random

centre <- capt_signal$Parameters$MeanBrownian
scale <- capt_signal$Parameters$MeanRandom - centre

brownian <- (brownian - centre) / scale
random   <- (random   - centre) / scale	
obs_capt      <- (capt_signal$Parameters$Observed - centre) / scale

Dat_Capt <- rbind(data.frame(Est = brownian, Model = "Brownian"), data.frame(Est = random, Model = "Model")) %>%
  mutate(Var = "Capt")

## Pull out data for plotting - wild
brownian <- wild_signal$Permutations$brownian
random   <- wild_signal$Permutations$random

centre <- wild_signal$Parameters$MeanBrownian
scale <- wild_signal$Parameters$MeanRandom - centre

brownian <- (brownian - centre) / scale
random   <- (random   - centre) / scale	
obs_wild      <- (wild_signal$Parameters$Observed - centre) / scale

Dat_Wild <- rbind(data.frame(Est = brownian, Model = "Brownian"), data.frame(Est = random, Model = "Model"))%>%
  mutate(Var = "Wild")


#### Plot phylogenetic signal ####

Capt_sig <- ggplot(Dat_Capt, aes(x = Est, y = Var, group = Model, fill = Model)) + 
  geom_density_ridges(stat = "binline", scale = .9, draw_baseline = FALSE, alpha = .5) +
  geom_vline(xintercept = c(obs_capt), colour = c("grey50"), 
             linetype = c("solid"), size = 2) +
  scale_fill_manual(values = c("dodgerblue", "black")) +
  coord_cartesian(ylim = c(1, 16), expand = FALSE) +
  #annotate("text", label = "Estimated D for captive trade", x = obs_capt - 0.1, y = 8, size = 4, colour = "grey50", angle = 90) +
  xlab("D-statistic") +
  theme_bw(base_size = 16) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", panel.grid = element_blank(), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color = 'black'))
  
Wild_sig <- ggplot(Dat_Wild, aes(x = Est, y = Var, group = Model, fill = Model)) + 
  geom_density_ridges(stat = "binline", scale = .9, draw_baseline = FALSE, alpha = .5) +
  geom_vline(xintercept = c(obs_wild), colour = c("chartreuse4"), 
             linetype = c("solid"), size = 2) +
  scale_fill_manual(values = c("dodgerblue", "black")) +
  coord_cartesian(ylim = c(1, 16), expand = FALSE) +
  #annotate("text", label = "Estimated D for captive trade", x = obs_wild - 0.1, y = 8, size = 4, colour = "chartreuse4", angle = 90) +
  xlab("D-statistic") +
  theme_bw(base_size = 16) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", panel.grid = element_blank(), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color = 'black'))


#### plot phylogeny ####
circ <- ggtree(phy_bird2, layout = "circular")

Listed_final2 <- column_to_rownames(Listed_final, "phyloname") %>% 
  mutate(Captive = ifelse(Captive == 1, "Yes", "No"),
         Wild = ifelse(Wild == 1, "Yes", "No"))

## For clarity pick out the most specious orders for labeling (species > 25)
Listed_final2 %>% group_by(Order.x) %>% tally() %>% filter(n>25)
Listed_final2 <- Listed_final2 %>% mutate(Order = case_when(Order.x == "ACCIPITRIFORMES" ~ "ACCIPITRIFORMES",
                                                            Order.x == "CAPRIMULGIFORMES" ~ "CAPRIMULGIFORMES",
                                                            Order.x == "FALCONIFORMES" ~ "FALCONIFORMES",
                                                            Order.x == "GALLIFORMES" ~ "GALLIFORMES",
                                                            Order.x == "PASSERIFORMES" ~ "PASSERIFORMES",
                                                            Order.x == "PSITTACIFORMES" ~ "PSITTACIFORMES",
                                                            Order.x == "STRIGIFORMES" ~ "STRIGIFORMES"))


## Inner ring of order
p1 <- gheatmap(circ, select(Listed_final2, Order), offset=.8, width=.1,
               colnames = FALSE, color = NA ) +
  scale_fill_manual(values=c("tan3", "black", "palegreen3", "dodgerblue", "mediumpurple2", "tomato", "khaki2"),
                    na.value = NA) +
  theme(legend.position = "none")

p2 <- p1 + new_scale_fill()

## Wild ring
p3 <-  gheatmap(p2, select(Listed_final2, Wild), offset=25, width=.2,
         colnames = FALSE, color = NA ) +
  scale_fill_manual(breaks=c("No", "Yes"), 
                    values=c("white", "chartreuse4"), name="Wild") +
  theme(legend.position = "none")

p4 <- p3 + new_scale_fill() 

## Captive ring
P5 <- gheatmap(p4, select(Listed_final2, Captive), offset=50, width=.2,
         colnames = FALSE, color = NA ) +
  scale_fill_manual(breaks=c("No", "Yes"), 
                    values=c("white", "grey50"), name="Captive") +
  theme(legend.position = "none")

#### Plotting species traits ####
Traits <- SLH_Final_phy %>% filter(n>0) %>% 
  group_by(WildSource, Name_for_CITESdb, Bodymass, Age_at_first_breeding, Max_longevity) %>%
  tally()

Summary_Traits <- Traits %>% group_by(WildSource) %>%
  summarise(meanBM = mean(Bodymass),minBM = min(Bodymass), maxBM = max(Bodymass),
            meanFR = mean(Age_at_first_breeding),minFR = min(Age_at_first_breeding), maxFR = max(Age_at_first_breeding),
            meanML = mean(Max_longevity),minML = min(Max_longevity), maxML = max(Max_longevity))

ML_Trait_plot <- ggplot(Traits, aes(x = Max_longevity, y = WildSource, fill = WildSource)) + 
  geom_density_ridges(stat = "binline", scale = .9, draw_baseline = FALSE) +
  geom_segment(data = Summary_Traits, aes(x = meanML, xend = meanML, y = c(1, 2), yend = c(1.9, 2.9)),
               colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = maxML, xend = maxML, y = c(1, 2), yend = c(1.9, 2.9)), 
               linetype = "dashed", colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = minML, xend = minML, y = c(1, 2), yend = c(1.9, 2.9)),
               linetype = "dashed", colour = "red", size = 1) +
  scale_fill_manual(values = c("grey50", "chartreuse4")) +
  coord_cartesian(expand = FALSE) +
  ylab("") +
  xlab("Maximum longevity (Years)") +
  theme_bw(base_size = 16) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
        legend.position = "none", plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color = 'black'))

FR_Trait_plot <- ggplot(Traits, aes(x = Age_at_first_breeding, y = WildSource, fill = WildSource)) + 
  geom_density_ridges(stat = "binline", scale = .9, draw_baseline = FALSE) +
  geom_segment(data = Summary_Traits, aes(x = meanFR, xend = meanFR, y = c(1, 2), yend = c(1.9, 2.9)),
               colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = maxFR, xend = maxFR, y = c(1, 2), yend = c(1.9, 2.9)), 
               linetype = "dashed", colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = minFR, xend = minFR, y = c(1, 2), yend = c(1.9, 2.9)),
               linetype = "dashed", colour = "red", size = 1) +
  scale_fill_manual(values = c("grey50", "chartreuse4")) +
  coord_cartesian(expand = FALSE) +
  ylab("") +
  xlab("First reproduction (Years)") +
  theme_bw(base_size = 16) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
        legend.position = "none", plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color = 'black'))

BM_Trait_plot <- ggplot(Traits, aes(x = Bodymass, y = WildSource, fill = WildSource)) + 
  geom_density_ridges(stat = "binline", scale = .9, draw_baseline = FALSE) +
  geom_segment(data = Summary_Traits, aes(x = meanBM, xend = meanBM, y = c(1, 2), yend = c(1.9, 2.9)),
               colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = maxBM, xend = maxBM, y = c(1, 2), yend = c(1.9, 2.9)), 
               linetype = "dashed", colour = "red", size = 1) +
  geom_segment(data = Summary_Traits, aes(x = minBM, xend = minBM, y = c(1, 2), yend = c(1.9, 2.9)),
               linetype = "dashed", colour = "red", size = 1) +
  scale_fill_manual(values = c("grey50", "chartreuse4")) +
  coord_cartesian(expand = FALSE) +
  ylab("") +
  scale_x_log10() +
  xlab("Body mass (g)") +
  theme_bw(base_size = 16) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
        legend.position = "none", plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color = 'black'))
  
  
#### Arrangement #####

Legend <- ggplot(Listed_final2, aes(Captive, Wild, fill = Order)) + geom_col()+
  scale_fill_manual(values=c("tan3", "black", "palegreen3", "dodgerblue", "mediumpurple2", "tomato", "khaki2"),
                    na.value = NA, na.translate=FALSE, 
                    labels = c("Accipitriformes", "Caprimulgiformes", "Falconiformes",
                              "Galliformes", "Passeriformes", "Psittaciformes",
                              "Strigiformes")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(), legend.position = "bottom")

library(ggpubr)

Leg <- get_legend(Legend)

Phy_plot1 <- ggarrange(P5, legend.grob = Leg, legend = "bottom", labels = "A.")

Sig_plot1 <- ggarrange(Capt_sig, Wild_sig, nrow = 1, labels = c("B.", "C."))

Trait_plot1 <- ggarrange(BM_Trait_plot, FR_Trait_plot, ML_Trait_plot, nrow = 1, labels = c("D.", "E.", "F."))


Sig_phy_plot <- ggarrange(Phy_plot1, Sig_plot1, Trait_plot1,  ncol = 1, heights = c(1, .4, .4))

ggsave("../1_Life_History/Outputs/Phy_Figure_alt.png", Sig_phy_plot, 
       width = 30, height = 30, units = "cm", device = "png", bg = "white")

