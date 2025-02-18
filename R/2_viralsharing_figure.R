# Sharing Model Code for Andrea.R ####

library(tidyverse)
PredictedNetwork <- readRDS("data/PredictedNetwork.rds")

PredictedNetwork %<>% as.matrix

# library(magrittr)
# vir_tab_tot_long %<>% mutate_at("spp", ~str_replace_all(.x, " ", "_"))
# 
# vir_tab_tot_long %<>% filter(spp != "Trichosurus_vulpecula")
# 
# vir_tab_tot_long$spp
# 
# PredictedNetwork[vir_tab_tot_long$spp, ]
# 
# PredictedNetwork %>% colnames %>% intersect(vir_tab_tot_long$spp)
# PredictedNetwork %>% colnames %>% setdiff(vir_tab_tot_long$spp, .)

alien_endemic_overlap <- read_csv("data/alien_endemic_overlap.csv")

alien_endemic_overlap %<>% filter(alien != "Macaca mulatta")

alien_endemic_overlap %>%
  mutate(alien = str_replace(alien, " ", "_"),
         endemic = str_replace(endemic, " ", "_")) -> alien_endemic

no_match_mammals <- PredictedNetwork %>% colnames %>% setdiff(alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[1]], .)
alien_endemic <- alien_endemic %>% filter(endemic %!in% no_match_mammals)

# PredictedNetwork[vir_tab_tot_long$spp[5], c("Marmota_marmota", "Zyzomys_argurus", "Procyon_lotor", "Canis_lupus", "Tamiasciurus_hudsonicus", "Sciurus_vulgaris")]

NVison <- PredictedNetwork[alien_endemic$alien[1], alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[1]]] 
NVison_tibble <- tibble(
  species = names(NVison),
  probability = NVison) %>% 
  mutate(Alien = "Neovison vison")

Ocuniculus <- PredictedNetwork[unique(alien_endemic$alien)[2], alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[2]]] 
Ocuniculus_tibble <- tibble(
  species = names(Ocuniculus),
  probability = Ocuniculus) %>% 
  mutate(Alien = "Oryctolagus cuniculus")
  
Leuropaeus <- PredictedNetwork[unique(alien_endemic$alien)[3], alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[3]]] 
Leuropaeus_tibble <- tibble(
  species = names(Leuropaeus),
  probability = Leuropaeus) %>% 
  mutate(Alien = "Lepus europaeus")

Lcapensis <- PredictedNetwork[unique(alien_endemic$alien)[4], alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[4]]] 
Lcapensis_tibble <- tibble(
  species = names(Lcapensis),
  probability = Lcapensis) %>% 
  mutate(Alien = "Lepus capensis")

Rtarandus <- PredictedNetwork[unique(alien_endemic$alien)[5], alien_endemic$endemic[alien_endemic$alien == alien_endemic$alien[5]]] 
Rtarandus_tibble <- tibble(
  species = names(Rtarandus),
  probability = Rtarandus) %>% 
  mutate(Alien = "Rangifer tarandus")

require(gghalves)
require("ggrepel")

rbind(Rtarandus_tibble, NVison_tibble, Ocuniculus_tibble, Leuropaeus_tibble, Lcapensis_tibble, Rtarandus_tibble) %>%
  distinct(species, Alien, probability) %>%
  mutate(labelling = ifelse(probability > 0.6, species, NA)) %>% 
  ggplot(aes(x = Alien, y = probability))+
  geom_point(position = position_jitter(width = 0.3), size = 2.5, shape = 20, alpha = 0.4, aes(color = probability))+
  geom_text_repel(aes(label = labelling), size = 2.8)+
  # geom_point(alpha = 1, position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.6), size = 2.5, shape = 20) +
  # geom_half_violin(lwd = 1, alpha = 0.5, position = position_dodge(width = 0.6), color = "#704926", fill = "#704926") +
  geom_boxplot(lwd = 1, color = "firebrick", fill = NA, outliers = F, width = 0.2)+
  stat_summary(fun.y = "mean", mapping = aes(y = probability, group = Alien), color = "black", geom = "point", size = 2.5) + #, position = position_dodge(width = 0.5)
  paletteer::scale_color_paletteer_c("ggthemes::Red")+
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Sharing probability", x = "", color = "", shape = "")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank())

ggsave("m/Alien_sharingprob_oldmodel.jpeg", width = 6, height = 6, dpi = 600)

