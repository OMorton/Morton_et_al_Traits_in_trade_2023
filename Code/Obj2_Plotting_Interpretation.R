######################################
#### Plotting and Model summaries ####
######################################

library(tidybayes)
library(tidyverse)

## plus model object
SLH_Final_phy <- read.csv("../1_Life_History/Outputs/SOM/CITES_Fitting_Data.csv")


##### New data ####
## ML
ML_dat1 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SMax_longevity = seq(from = min(SMax_longevity), to = max(SMax_longevity), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SFirst_repro = 0,
         SYear = (20 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "20")

ML_dat2 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SMax_longevity = seq(from = min(SMax_longevity), to = max(SMax_longevity), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SFirst_repro = 0,
         SYear = (0 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "0")

ML_dat3 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SMax_longevity = seq(from = min(SMax_longevity), to = max(SMax_longevity), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SFirst_repro = 0,
         SYear = (10 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "10")

ML_fit <- add_epred_draws(All1, newdata =  rbind(ML_dat1, ML_dat2, ML_dat3), 
                          re_formula = ~(1+ WildThreat | Year), dpar = "hu")

ML_fit_sum <- ML_fit %>% filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SMax_longevity, Year) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(Max_longevity = SMax_longevity*sd(SLH_Final_phy$LogMax_longevity) + mean(SLH_Final_phy$LogMax_longevity),
         Max_longevity = 2^Max_longevity,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

## FR
FR_dat1 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SFirst_repro = seq(from = min(SFirst_repro), to = max(SFirst_repro), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SMax_longevity = 0,
         SYear = (20 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "20")

FR_dat2 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SFirst_repro = seq(from = min(SFirst_repro), to = max(SFirst_repro), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SMax_longevity = 0,
         SYear = (0 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "0")

FR_dat3 <- SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SFirst_repro = seq(from = min(SFirst_repro), to = max(SFirst_repro), length.out = 20)) %>%
  mutate(SLogBodymass = 0,
         SMax_longevity = 0,
         SYear = (10 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "10")

FR_fit <- add_epred_draws(All1, newdata =  rbind(FR_dat1, FR_dat2, FR_dat3), 
                          re_formula = ~(1+ WildThreat | Year), dpar = "hu")

FR_fit_sum <- FR_fit %>% filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SFirst_repro, Year) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(SFirst_repro = SFirst_repro*sd(SLH_Final_phy$LogFirst_repro) + mean(SLH_Final_phy$LogFirst_repro),
         SFirst_repro = 2^SFirst_repro,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))



## BM
BM_dat1 <-SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SLogBodymass = seq(from = min(SLogBodymass), to = max(SLogBodymass), length.out = 20)) %>%
  mutate(SFirst_repro = 0,
         SMax_longevity = 0,
         SYear = (20 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "20")

BM_dat2 <-SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SLogBodymass = seq(from = min(SLogBodymass), to = max(SLogBodymass), length.out = 20)) %>%
  mutate(SFirst_repro = 0,
         SMax_longevity = 0,
         SYear = (0 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "0")

BM_dat3 <-SLH_Final_phy %>%  group_by(WildSource, Threat_code) %>%
  summarise(SLogBodymass = seq(from = min(SLogBodymass), to = max(SLogBodymass), length.out = 20)) %>%
  mutate(SFirst_repro = 0,
         SMax_longevity = 0,
         SYear = (10 - mean(as.numeric(as.character(SLH_Final_phy$Year))))/sd(as.numeric(as.character(SLH_Final_phy$Year))),
         Year = "10")

BM_fit <- add_epred_draws(All1, newdata =  rbind(BM_dat1, BM_dat2, BM_dat3), 
                          re_formula = ~(1+ WildThreat | Year), dpar = "hu")

BM_fit_sum <- BM_fit %>% filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SLogBodymass, Year) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(SLogBodymass = SLogBodymass*sd(log2(SLH_Final_phy$Bodymass)) + mean(log2(SLH_Final_phy$Bodymass)),
         SLogBodymass = 2^SLogBodymass,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))



## Year ML plot
ML_20_CB <- ML_fit_sum %>% filter(Year == 20, Source == "Captive")
ML_10_CB <- ML_fit_sum %>% filter(Year == 10, Source == "Captive")
ML_0_CB <- ML_fit_sum %>% filter(Year == 0, Source == "Captive")
ML_20_W <- ML_fit_sum %>% filter(Year == 20, Source == "Wild")
ML_10_W <- ML_fit_sum %>% filter(Year == 10, Source == "Wild")
ML_0_W <- ML_fit_sum %>% filter(Year == 0, Source == "Wild")

(ML_CB_plot <- ggplot(ML_20_CB, aes(Max_longevity, .epred, colour = Source, fill = Source)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
  geom_line(aes(), size = 1.25) +
  geom_line(data = ML_0_CB, size = 1.25, linetype = "dashed") +
  geom_line(data = ML_10_CB, size = 1.25, linetype = "longdash") +
  ylab("Volume estimates") + 
  xlab("Maximum longevity (Years)") +
  scale_colour_manual(values = c("grey50")) +
  scale_fill_manual(values = c("grey50")) +
  scale_y_log10(breaks = c(0.001, 0.1, 10), labels = c(0.001, 0.1, 10)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.position = "none", strip.text = element_blank()))

(ML_W_plot <- ggplot(ML_20_W, aes(Max_longevity, .epred, colour = Source, fill = Source)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
  geom_line(aes(), size = 1.25) +
  geom_line(data = ML_0_W, size = 1.25, linetype = "dashed") +
  geom_line(data = ML_10_W, size = 1.25, linetype = "longdash") +
  ylab("Volume estimates") + 
  xlab("Maximum longevity (Years)") +
  scale_colour_manual(values = c("chartreuse4")) +
  scale_fill_manual(values = c("chartreuse4")) +
  annotate("text", x = 70, y = 0.4, label = "2000", size = 5) +
  annotate("text", x = 70, y = 0.006, label = "2010", size = 5) +
  annotate("text", x = 70, y = 0.0003, label = "2020", size = 5) +
  scale_y_log10(breaks = c(0.00001, 0.01, 10), labels = c(0.00001, 0.01, 10)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.position = "none", strip.text = element_blank(), axis.title.y = element_blank()))





## Year BM plot
BM_20_CB <- BM_fit_sum %>% filter(Year == 20, Source == "Captive")
BM_10_CB <- BM_fit_sum %>% filter(Year == 10, Source == "Captive")
BM_0_CB <- BM_fit_sum %>% filter(Year == 0, Source == "Captive")
BM_20_W <- BM_fit_sum %>% filter(Year == 20, Source == "Wild")
BM_10_W <- BM_fit_sum %>% filter(Year == 10, Source == "Wild")
BM_0_W <- BM_fit_sum %>% filter(Year == 0, Source == "Wild")

(BM_CB_plot <- ggplot(BM_20_CB, aes(SLogBodymass , .epred, colour = Source, fill = Source)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
    geom_line(aes(), size = 1.25) +
    geom_line(data = BM_0_CB, size = 1.25, linetype = "dashed") +
    geom_line(data = BM_10_CB, size = 1.25, linetype = "longdash") +
    ylab("Volume estimates") + 
    xlab("Body mass (g)") +
    scale_colour_manual(values = c("grey50")) +
    scale_fill_manual(values = c("grey50")) +
    scale_y_log10(breaks = c(0.001, 0.1, 10), labels = c(0.001, 0.1, 10), limits = c(0.001, 100)) +
    scale_x_log10() +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank()))

(BM_W_plot <- ggplot(BM_20_W, aes(SLogBodymass , .epred, colour = Source, fill = Source)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
  geom_line(aes(), size = 1.25) +
  geom_line(data = BM_0_W, size = 1.25, linetype = "dashed") +
  geom_line(data = BM_10_W, size = 1.25, linetype = "longdash") +
  ylab("Volume estimates") + 
  xlab("Body mass (g)") +
  annotate("text", x = 10, y = 3, label = "2000", size = 5) +
  annotate("text", x = 10, y = 0.01, label = "2010", size = 5) +
  annotate("text", x = 10, y = 0.00004, label = "2020", size = 5) +
  scale_colour_manual(values = c("chartreuse4")) +
  scale_fill_manual(values = c("chartreuse4")) +
  scale_y_log10(breaks = c(0.00001, 0.01, 10), labels = c(0.00001, 0.01, 10)) +
  scale_x_log10() +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.position = "none", strip.text = element_blank(), axis.title.y = element_blank()))

ggplot(data = filter(SLH_Final_phy, WildSource == "Yes", n == 0), aes(Bodymass)) + 
  geom_histogram() + scale_x_log10()

## Year FR plot
FR_20_CB <- FR_fit_sum %>% filter(Year == 20, Source == "Captive")
FR_10_CB <- FR_fit_sum %>% filter(Year == 10, Source == "Captive")
FR_0_CB <- FR_fit_sum %>% filter(Year == 0, Source == "Captive")

FR_20_W <- FR_fit_sum %>% filter(Year == 20, Source == "Wild")
FR_10_W <- FR_fit_sum %>% filter(Year == 10, Source == "Wild")
FR_0_W <- FR_fit_sum %>% filter(Year == 0, Source == "Wild")

(FR_CB_plot <- ggplot(FR_20_CB, aes(SFirst_repro , .epred, colour = Source, fill = Source)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
    geom_line(aes(), size = 1.25) +
    geom_line(data = FR_0_CB, size = 1.25, linetype = "dashed") +
    geom_line(data = FR_10_CB, size = 1.25, linetype = "longdash") +
    ylab("Volume estimates") + 
    xlab("First reproduction (Years)") +
    scale_colour_manual(values = c("grey50")) +
    scale_fill_manual(values = c("grey50")) +
    scale_y_log10(breaks = c(0.001, 0.1, 10), labels = c(0.001, 0.1, 10)) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank()))

(FR_W_plot <- ggplot(FR_20_W, aes(SFirst_repro , .epred, colour = Source, fill = Source)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15, colour = NA) +
    geom_line(aes(), size = 1.25) +
    geom_line(data = FR_0_W, size = 1.25, linetype = "dashed") +
    geom_line(data = FR_10_W, size = 1.25, linetype = "longdash") +
    ylab("Volume estimates") + 
    xlab("First reproduction (Years)") +
    scale_colour_manual(values = c("chartreuse4")) +
    scale_fill_manual(values = c("chartreuse4")) +
    scale_y_log10(breaks = c(0.00001, 0.01, 10), labels = c(0.00001, 0.01, 10)) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "tomato", size = 1) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.y = element_blank()))


#### Figure 3 simplistic arrangement ####
library(ggpubr)

Fig3_v1 <- ggarrange(BM_CB_plot, BM_W_plot, FR_CB_plot, FR_W_plot, ML_CB_plot, ML_W_plot,
          align = "hv", ncol = 2, nrow = 3, labels = c("A.", "B.", "C.", "D.", "E.", "F."))

Fig3_v2 <- annotate_figure( Fig3_v1, left = text_grob("Estimated volumes (WOEs)", color = "black", rot = 90))

empty <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "white", colour = NA))

Figure_3 <- ggarrange(ggarrange(empty,empty, nrow = 1, ncol = 2, 
                                labels = c("Captive-sourced", "Wild-sourced"), hjust = c(-1, -1.4),
                                font.label = list(face = "bold.italic")), 
                      Fig3_v2, 
                      nrow = 2, heights = c(.03, 1))


##### hu and mu fits ####

## mu
ML_hu <- add_epred_draws(All1, newdata =  rbind(ML_dat1, ML_dat2, ML_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "hu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SMax_longevity, Year) %>% 
  mutate(hu = 1 - hu) %>%
  median_hdci(hu, .width = .9)

FR_hu <- add_epred_draws(All1, newdata =  rbind(FR_dat1, FR_dat2, FR_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "hu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SFirst_repro, Year) %>% 
  mutate(hu = 1 - hu) %>%
  median_hdci(hu, .width = .9)

BM_hu <- add_epred_draws(All1, newdata =  rbind(BM_dat1, BM_dat2, BM_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "hu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SLogBodymass, Year) %>% 
  mutate(hu = 1 - hu) %>%
  median_hdci(hu, .width = .9)

## mu
ML_mu <- add_epred_draws(All1, newdata =  rbind(ML_dat1, ML_dat2, ML_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "mu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SMax_longevity, Year) %>% 
  median_hdci(mu, .width = .9)

FR_mu <- add_epred_draws(All1, newdata =  rbind(FR_dat1, FR_dat2, FR_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "mu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SFirst_repro, Year) %>% 
  median_hdci(mu, .width = .9)

BM_mu <- add_epred_draws(All1, newdata =  rbind(BM_dat1, BM_dat2, BM_dat3), 
                         re_formula = ~(1+ WildThreat | Year), dpar = "mu") %>% 
  filter(Threat_code == "Non-threatened") %>%
  group_by(WildSource, Threat_code, SLogBodymass, Year) %>% 
  median_hdci(mu, .width = .9)

## Tidy
## Mu
ML_mu_sum <- ML_mu_sum%>% 
  mutate(Max_longevity = SMax_longevity*sd(SLH_Final_phy$LogMax_longevity) + mean(SLH_Final_phy$LogMax_longevity),
         Max_longevity = 2^Max_longevity,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

BM_mu_sum <- BM_mu_sum %>% 
  mutate(SLogBodymass = SLogBodymass*sd(log2(SLH_Final_phy$Bodymass)) + mean(log2(SLH_Final_phy$Bodymass)),
         SLogBodymass = 2^SLogBodymass,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

FR_mu_sum <- FR_mu_sum %>% 
  mutate(SFirst_repro = SFirst_repro*sd(SLH_Final_phy$LogFirst_repro) + mean(SLH_Final_phy$LogFirst_repro),
         SFirst_repro = 2^SFirst_repro,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

## Hu
ML_hu_sum <- ML_hu_sum%>% 
  mutate(Max_longevity = SMax_longevity*sd(SLH_Final_phy$LogMax_longevity) + mean(SLH_Final_phy$LogMax_longevity),
         Max_longevity = 2^Max_longevity,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

BM_hu_sum <- BM_hu_sum %>% 
  mutate(SLogBodymass = SLogBodymass*sd(log2(SLH_Final_phy$Bodymass)) + mean(log2(SLH_Final_phy$Bodymass)),
         SLogBodymass = 2^SLogBodymass,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))

FR_hu_sum <- FR_hu_sum %>% 
  mutate(SFirst_repro = SFirst_repro*sd(SLH_Final_phy$LogFirst_repro) + mean(SLH_Final_phy$LogFirst_repro),
         SFirst_repro = 2^SFirst_repro,
         Source = ifelse(WildSource == "Yes", "Wild", "Captive"))


#### Plotting Hu and Mu ####
SLH_Final_phy <- SLH_Final_phy %>% mutate(Pres = ifelse(n > 0, 1, 0))

## ML
(ML_mu <- ggplot(data = filter(SLH_Final_phy, n >0), aes(x = Max_longevity, n, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(ML_mu_sum, Year == "0"), aes(Max_longevity, mu, colour = WildSource), size = 1) +
    geom_ribbon(data = filter(ML_mu_sum, Year == "0"), aes(Max_longevity, mu, ymin = .lower, ymax = .upper, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_y_log10(breaks = c(10, 1000, 100000), labels = c(10, 1000, 100000)) +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    labs(y = expression(italic("Mu"))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

(ML_hu <- ggplot(data = filter(SLH_Final_phy), aes(x = Max_longevity, Pres, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(ML_hu_sum, Year == "0"), aes(Max_longevity, hu, colour = WildSource), size = 1) +
    geom_ribbon(data = filter(ML_hu_sum, Year == "0"), aes(Max_longevity, hu, ymin = .lower, ymax = .upper, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    labs(y = expression(italic("Hu"))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

## FR
(FR_mu <- ggplot(data = filter(SLH_Final_phy, n >0), aes(x = Age_at_first_breeding, n, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(FR_mu_sum, Year == "0"), aes(SFirst_repro, mu, colour = WildSource), size = 1) +
    geom_ribbon(data = filter(FR_mu_sum, Year == "0"), aes(SFirst_repro, mu,ymin = .lower, ymax = .upper, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_y_log10(breaks = c(10, 1000, 100000), labels = c(10, 1000, 100000)) +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    labs(y = expression(italic("Mu"))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

(FR_hu <- ggplot(data = filter(SLH_Final_phy), aes(x = Age_at_first_breeding, Pres, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(FR_hu_sum, Year == "0"), aes(SFirst_repro, hu, colour = WildSource), size = 1) +
    geom_ribbon(data = filter(FR_hu_sum, Year == "0"), aes(SFirst_repro, hu, ymin = .lower, ymax = .upper, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    labs(y = expression(italic("Hu"))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

## BM
(BM_mu <- ggplot(data = filter(SLH_Final_phy, n >0), aes(x = Bodymass, n, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(BM_mu_sum, Year == "0"), aes(SLogBodymass, mu, colour = WildSource), size =1) +
    geom_ribbon(data = filter(BM_mu_sum, Year == "0"), aes(SLogBodymass, ymin = .lower, ymax = .upper, mu, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_y_log10(breaks = c(10, 1000, 100000), labels = c(10, 1000, 100000)) +
    scale_x_log10() +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    labs(y = expression(italic("Mu"))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

(BM_hu <- ggplot(data = filter(SLH_Final_phy), aes(x = Bodymass, Pres, colour = WildSource)) + 
    geom_point(shape = 16, alpha = .02) +
    geom_line(data = filter(BM_hu_sum, Year == "0"), aes(SLogBodymass, hu, colour = WildSource), size = 1) +
    geom_ribbon(data = filter(BM_hu_sum, Year == "0"), aes(SLogBodymass, ymin = .lower, ymax = .upper, hu, fill = WildSource),
                colour = NA, alpha = .1) +
    scale_color_manual(values = c("grey50", "chartreuse4")) +
    scale_fill_manual(values = c("grey50", "chartreuse4")) +
    scale_x_log10() +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(y = expression(italic("Hu"))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", strip.text = element_blank(), axis.title.x = element_blank()))

#### Figure 3 arrangement with data + hu/mu ####
ML_HuMu <- annotate_figure(ggarrange(ML_hu, ML_mu, nrow = 2, ncol = 1, align = "hv"),
                           bottom = text_grob("", color = "black"))

ML_fits <- annotate_figure(ggarrange(ML_CB_plot, ML_W_plot, ncol = 2, labels = c("H.", "I.")),
                           left = text_grob("Estimated volumes (WOEs)", color = "black", rot = 90))

## FR
FR_HuMu <- annotate_figure(ggarrange(FR_hu, FR_mu, nrow = 2, ncol = 1, align = "hv"),
                           bottom = text_grob("", color = "black"))

FR_fits <- annotate_figure(ggarrange(FR_CB_plot, FR_W_plot, ncol = 2, labels = c("E.", "F.")),
                           left = text_grob("Estimated volumes (WOEs)", color = "black", rot = 90))


## BM
BM_HuMu <- annotate_figure(ggarrange(BM_hu, BM_mu, nrow = 2, ncol = 1, align = "hv"),
                           bottom = text_grob("", color = "black"))

BM_fits <- annotate_figure(ggarrange(BM_CB_plot, BM_W_plot, ncol = 2, labels = c("B.", "C.")),
                           left = text_grob("Estimated volumes (WOEs)", color = "black", rot = 90))

## final
all_fits <- ggarrange(BM_CB_plot, BM_W_plot, FR_CB_plot, FR_W_plot, ML_CB_plot, ML_W_plot, ncol = 2, nrow = 3,
                      labels = c("B.", "C.", "E.", "F.", "H.", "I."), align = "hv")


Figure3_humu <- ggarrange(BM_HuMu, FR_HuMu, ML_HuMu, ncol = 1, nrow = 3, align = "hv", 
                          labels = c("A.", "D.", "G."))

Figure3_points <- ggarrange(Figure3_humu, all_fits, ncol = 2, widths = c(1, 2))


ggsave("../Outputs/Figure3/Figure3_Points.png", Figure3_points, 
       width = 25, height = 25, units = "cm", device = "png", bg = "white")



#### Fixed effects - Main####

Fixef_ML <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SMax_longevity")) %>%
  summarise(Wild_mu = `SMax_longevity` + `WildSourceYes:SMax_longevity`, 
            Captive_mu = SMax_longevity,
            Wild_hu = `hu_SMax_longevity`*-1 + `hu_WildSourceYes:SMax_longevity`*-1, 
            Captive_hu = hu_SMax_longevity*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "ML")

Fixef_ML_Sum <- Fixef_ML %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)


Fixef_BM <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SLogBodymass")) %>%
  summarise(Wild_mu = `SLogBodymass` + `WildSourceYes:SLogBodymass`, 
            Captive_mu = SLogBodymass,
            Wild_hu = `hu_SLogBodymass`*-1 + `hu_WildSourceYes:SLogBodymass`*-1, 
            Captive_hu = hu_SLogBodymass*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "BM")

Fixef_BM_Sum <- Fixef_BM %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)


Fixef_AFR <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SFirst_repro")) %>%
  summarise(Wild_mu = `SFirst_repro` + `WildSourceYes:SFirst_repro`, 
            Captive_mu = SFirst_repro,
            Wild_hu = `hu_SFirst_repro`*-1 + `hu_WildSourceYes:SFirst_repro`*-1, 
            Captive_hu = hu_SFirst_repro*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "AFR")

Fixef_AFR_Sum <- Fixef_AFR %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)



#### Alt plot ####
All_coefs <- rbind(Fixef_BM_Sum, Fixef_AFR_Sum, Fixef_ML_Sum)
All_coefs_raw <- rbind(Fixef_BM, Fixef_AFR, Fixef_ML)

library(ggridges)
library(ggstance)

(Coef_plot <- ggplot(All_coefs, aes(Est,Trait, colour = Source, fill = Source)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = "grey95", colour = NA) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5, fill = "grey95", colour = NA) +
  geom_point(position = position_dodge(.5), size = 3) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), position = position_dodge(0.5), height = 0, size = .75) +
  scale_colour_manual(values = c("grey50", "chartreuse4")) +
  facet_wrap(~dpar, labeller = as_labeller(c(`hu` = "Reoccurence", `mu` = "Volume when traded"))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Parameter estimate") +
  scale_y_discrete(breaks = c("BM", "AFR", "ML"), 
                   labels = expression(beta["Body mass"], 
                                       beta["First reproduction"],
                                       beta["Longevity"])) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), axis.text.y = element_text(colour = "black", size = 22)))

#### Fixed effects - Interaction ####

## ML
Fixef_ML_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SMax_longevity")) %>%
  summarise(Wild_mu = `SYear:SMax_longevity` + `WildSourceYes:SYear:SMax_longevity`, 
            Captive_mu = `SYear:SMax_longevity`,
            Wild_hu = `hu_SYear:SMax_longevity`*-1 + `hu_WildSourceYes:SYear:SMax_longevity`*-1, 
            Captive_hu = `hu_SYear:SMax_longevity`*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "ML") 

Fixef_ML_Int_Sum <- Fixef_ML_Int %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)

## BM
Fixef_BM_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SLogBodymass")) %>%
  summarise(Wild_mu = `SYear:SLogBodymass` + `WildSourceYes:SYear:SLogBodymass`, 
            Captive_mu = `SYear:SLogBodymass`,
            Wild_hu = `hu_SYear:SLogBodymass`*-1 + `hu_WildSourceYes:SYear:SLogBodymass`*-1, 
            Captive_hu = `hu_SYear:SLogBodymass`*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "BM") 

Fixef_BM_Int_Sum <- Fixef_BM_Int %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)

## AFR
Fixef_AFR_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SFirst_repro")) %>%
  summarise(Wild_mu = `SYear:SFirst_repro` + `WildSourceYes:SYear:SFirst_repro`, 
            Captive_mu = `SYear:SFirst_repro`,
            Wild_hu = `hu_SYear:SFirst_repro`*-1 + `hu_WildSourceYes:SYear:SFirst_repro`*-1, 
            Captive_hu = `hu_SYear:SFirst_repro`*-1) %>%
  pivot_longer(everything(), names_to = "Par", values_to = "Est") %>%
  mutate(dpar = case_when(grepl("hu", Par) ~ "hu",
                          grepl("mu", Par) ~ "mu"),
         Source = case_when(grepl("Captive", Par) ~ "Captive",
                            grepl("Wild", Par) ~ "Wild"),
         Trait = "AFR") 

Fixef_AFR_Int_Sum <- Fixef_AFR_Int %>% group_by(Trait, Source, dpar) %>% median_hdci(Est, .width = .9)

#### Alt plot ####
All_coefs_Int <- rbind(Fixef_BM_Int_Sum, Fixef_AFR_Int_Sum, Fixef_ML_Int_Sum)
All_coefs_Int_raw <- rbind(Fixef_BM_Int, Fixef_AFR_Int, Fixef_ML_Int)


(Coef_plot_Int <- ggplot(All_coefs_Int, aes(Est,Trait, colour = Source, fill = Source)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = "grey95", colour = NA) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5, fill = "grey95", colour = NA) +
    geom_point(position = position_dodge(.5), size = 3) +
    geom_errorbarh(aes(xmin = .lower, xmax = .upper), position = position_dodge(0.5), height = 0, size = .75) +
    scale_colour_manual(values = c("grey50", "chartreuse4")) +
    facet_wrap(~dpar, labeller = as_labeller(c(`hu` = "Reoccurence", `mu` = "Volume when traded"))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ylab("") + xlab("Parameter estimate") +
    scale_y_discrete(breaks = c("BM", "AFR", "ML"), 
                     labels = expression(beta["(Body mass x Year)"], 
                                         beta["(First reproduction x Year)"],
                                         beta["(Longevity x Year)"])) +
    theme_bw(base_size = 20) +
    theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), axis.text.y = element_text(colour = "black", size = 22)))


#### ICC - Variance decomp ####
PPD_Fix <- posterior_predict(All1,  re_formula = NA, summary = FALSE)
var_Fixed <- apply(PPD_Fix, MARGIN = 1, FUN = stats::var)
median(var_Fixed)

PPD_Full <- posterior_predict(All1,  re_formula = NULL, summary = FALSE)
var_Full <- apply(PPD_Full, MARGIN = 1, FUN = stats::var)
median(var_Full)

Fix_var_sum <- data.frame(Var = var_Fixed/var_Full) %>% median_hdci(Var, .width = .9) %>% mutate(Comp = "Fixed")

v1 <- data.frame(x = var_Fixed) %>% median_hdci(x, .width = .9) %>% mutate(Comp = "Fixed effects only")
v2 <- data.frame(x = var_Full) %>% median_hdci(x, .width = .9) %>% mutate(Comp = "Full model")

Var_plot <- ggplot(rbind(v1, v2), aes(x, Comp, colour = Comp)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0, size = 1) +
  scale_x_log10(breaks = c(10, 10000, 10000000)) +
  xlab(expression(sigma^2 ~ "recovered")) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none", panel.grid = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black"))


#### Figure 4 Arrangement ####

library(ggpubr)

(Fig4 <- cowplot::plot_grid(Coef_plot, Coef_plot_Int, Var_plot,  
                           ncol = 1, nrow = 3, align = "v", axis = "l", rel_heights = c(1, 1,.5),
                           labels = c("A. Trait trade dynamics", 
                                      "B. Trait-time trade dynamics", 
                                      "C."),  
                           label_size = 18, hjust = c(0,0,0)))


ggsave(".../1_Life_History/Outputs/Figure4/Figure4_NOAS.png", Fig4, 
       width = 35, height = 30, units = "cm", device = "png")


#### Write out ###
library(bayestestR)

PD_ML <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SMax_longevity")) %>%
  summarise(Wild_mu = `SMax_longevity` + `WildSourceYes:SMax_longevity`, 
            Captive_mu = SMax_longevity,
            Wild_hu = `hu_SMax_longevity`*-1 + `hu_WildSourceYes:SMax_longevity`*-1, 
            Captive_hu = hu_SMax_longevity*-1) %>%
  p_direction()

PD_BM <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SLogBodymass")) %>%
  summarise(Wild_mu = `SLogBodymass` + `WildSourceYes:SLogBodymass`, 
            Captive_mu = SLogBodymass,
            Wild_hu = `hu_SLogBodymass`*-1 + `hu_WildSourceYes:SLogBodymass`*-1, 
            Captive_hu = hu_SLogBodymass*-1) %>%
  p_direction()

PD_AFR <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SFirst_repro")) %>%
  summarise(Wild_mu = `SFirst_repro` + `WildSourceYes:SFirst_repro`, 
            Captive_mu = SFirst_repro,
            Wild_hu = `hu_SFirst_repro`*-1 + `hu_WildSourceYes:SFirst_repro`*-1, 
            Captive_hu = hu_SFirst_repro*-1) %>% 
  p_direction()

Fixef_BM_Sum2 <- Fixef_BM_Sum %>% unite("Parameter", 2:3, sep = "_")
Fixef_AFR_Sum2 <- Fixef_AFR_Sum %>% unite("Parameter", 2:3, sep = "_")
Fixef_ML_Sum2 <- Fixef_ML_Sum %>% unite("Parameter", 2:3, sep = "_")

BM_out <- left_join(PD_BM, Fixef_BM_Sum2)
AFR_out <- left_join(PD_AFR, Fixef_AFR_Sum2)
ML_out <- left_join(PD_ML, Fixef_ML_Sum2)

write.csv(BM_out, "../1_Life_History/Outputs/Figure4/BM_Sum.csv")
write.csv(AFR_out, "../1_Life_History/Outputs/Figure4/AFR_Sum.csv")
write.csv(ML_out, "../1_Life_History/Outputs/Figure4/ML_Sum.csv")

## Int
PD_ML_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SMax_longevity")) %>%
  summarise(Wild_mu = `SYear:SMax_longevity` + `WildSourceYes:SYear:SMax_longevity`, 
            Captive_mu = `SYear:SMax_longevity`,
            Wild_hu = `hu_SYear:SMax_longevity`*-1 + `hu_WildSourceYes:SYear:SMax_longevity`*-1, 
            Captive_hu = `hu_SYear:SMax_longevity`*-1) %>%
  p_direction() %>% as.data.frame()

PD_BM_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SLogBodymass")) %>%
  summarise(Wild_mu = `SYear:SLogBodymass` + `WildSourceYes:SYear:SLogBodymass`, 
            Captive_mu = `SYear:SLogBodymass`,
            Wild_hu = `hu_SYear:SLogBodymass`*-1 + `hu_WildSourceYes:SYear:SLogBodymass`*-1, 
            Captive_hu = `hu_SYear:SLogBodymass`*-1) %>%
  p_direction() %>% as.data.frame()

PD_AFR_Int <- fixef(All1, summary = FALSE) %>% as.data.frame() %>% 
  select(contains("SFirst_repro")) %>%
  summarise(Wild_mu = `SYear:SFirst_repro` + `WildSourceYes:SYear:SFirst_repro`, 
            Captive_mu = `SYear:SFirst_repro`,
            Wild_hu = `hu_SYear:SFirst_repro`*-1 + `hu_WildSourceYes:SYear:SFirst_repro`*-1, 
            Captive_hu = `hu_SYear:SFirst_repro`*-1) %>%
  p_direction() %>% as.data.frame()

Fixef_BM_Int_Sum2 <- Fixef_BM_Int_Sum %>% unite("Parameter", 2:3, sep = "_")
Fixef_AFR_Int_Sum2 <- Fixef_AFR_Int_Sum %>% unite("Parameter", 2:3, sep = "_")
Fixef_ML_Int_Sum2 <- Fixef_ML_Int_Sum %>% unite("Parameter", 2:3, sep = "_")

BM_Int_out <- left_join(PD_BM_Int, Fixef_BM_Int_Sum2)
AFR_Int_out <- left_join(PD_AFR_Int, Fixef_AFR_Int_Sum2)
ML_Int_out <- left_join(PD_ML_Int, Fixef_ML_Int_Sum2)

write.csv(BM_Int_out, "../1_Life_History/Outputs/Figure4/BM_Int_Sum.csv")
write.csv(AFR_Int_out, "../1_Life_History/Outputs/Figure4/AFR_Int_Sum.csv")
write.csv(ML_Int_out, "../1_Life_History/Outputs/Figure4/ML_Int_Sum.csv")
