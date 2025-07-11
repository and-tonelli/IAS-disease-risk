# Assembing training and prediction datasets

require(tidyverse)

# Min-Max scaling function
scale_minmax <- function(column){
  
  max_val <- max(column)
  min_val <- min(column)
  
  return((column - min_val)/(max_val - min_val))
  
}

# "not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

#### Training dataset ####
# Loading viral sharing according to VIRION
Sharing_all <- read_csv("Data/Sharing_all.csv")

# Loading native-native overlap
endemic_network <- read.csv("Data/endemic_endemic_network.csv") %>% select(c(1, 2, 6))

# Loading phylo distances
phylo_dist <- read.csv("data/PHYLACINE_pairwaise_div_time_mammals.csv") %>% select(-3)

# Loading pairwise biological distances
SpeciesPairsDistances <- read_csv("Data/SpeciesPairsDistances.csv")

# Loading citation counts
cit_virus <- read.csv("Data/cit_virus.csv")
cit_virus <- rbind(cit_virus, tibble(X = "Macaca fascicularis", V1 = as.numeric(1084))) # taken from Albery et al., 2020

# Merging all distances with the sharing network
inverted_pairs <- SpeciesPairsDistances
inverted_pairs <- inverted_pairs %>%
  rename(Species1_temp = Species1,
         Species2_temp = Species2) %>%
  mutate(Species1 = Species2_temp,
         Species2 = Species1_temp) %>%
  select(-Species1_temp, -Species2_temp)

SpeciesPairsDistances_mirrored <- rbind(SpeciesPairsDistances, inverted_pairs)

species_dist_phylo <- left_join(SpeciesPairsDistances_mirrored, phylo_dist, join_by(Species1 == sp, Species2 == item2))

species_dist_phylo %<>% na.omit()

Sharing_all_inv <- Sharing_all
Sharing_all_inv <- Sharing_all_inv %>%
  rename(Species1_temp = Species1,
         Species2_temp = Species2) %>%
  mutate(Species1 = Species2_temp,
         Species2 = Species1_temp) %>%
  select(-Species1_temp, -Species2_temp)

Sharing_all_mirrored <- rbind(Sharing_all, Sharing_all_inv)

species_dist_phylo_sharing <- left_join(species_dist_phylo, Sharing_all_mirrored, join_by(Species1 == Species1, Species2 == Species2))

species_dist_phylo_sharing <- species_dist_phylo_sharing %>% 
  mutate(n_virus = ifelse(is.na(n_virus), 0, n_virus))

endemic_network_in_virion <- endemic_network %>% 
  mutate(Species1 = str_replace(Species1, "_", " "),
         Species2 = str_replace(Species2, "_", " ")) %>% 
  filter(Species1 %in% c(Sharing_all$Species1, Sharing_all$Species2) & Species2 %in% c(Sharing_all$Species1, Sharing_all$Species2)) # Only including species with at least one sharing

endemic_network_in_virion_mirror <- endemic_network_in_virion %>% 
  rename(Species1 = Species2,
         Species2 = Species1)

# Complete network (duplicated swapping Species1 and Species2)
endemic_network_mirrored <- rbind(endemic_network_in_virion_mirror, endemic_network_in_virion)

# Getting all combinations (combinations without known sharing become pseudoabsences)
grid_combinations <- expand.grid(Species1 = unique(endemic_network_in_virion$Species1),
                                 Species2 = unique(endemic_network_in_virion$Species2))

NativeDataset <- left_join(grid_combinations, species_dist_phylo_sharing, join_by(Species1 == Species1, Species2 == Species2), )
NativeDataset_fin <- left_join(NativeDataset, endemic_network_mirrored, join_by(Species1 == Species1, Species2 == Species2))
NativeDataset_fin %<>% filter(Species1 != Species2)
NativeDataset_fin$Overlapping_Cells[is.na(NativeDataset_fin$Overlapping_Cells)] <- 0

Final_df_sc <- left_join(NativeDataset_fin, cit_virus %>% rename("n_cit1" = V1), join_by(Species1 == X))
Final_df_sc <- left_join(Final_df_sc, cit_virus %>% rename("n_cit2" = V1), join_by(Species2 == X))

NativeDataset_fin <- Final_df_sc %>%
  mutate(Pair = pmin(Species1, Species2), Pair2 = pmax(Species1, Species2)) %>%
  distinct(Pair, Pair2, .keep_all = TRUE) %>%
  select(-Pair, -Pair2)

NativeDataset_fin$n_virus[is.na(NativeDataset_fin$n_virus)] <- 0

NativeDataset_fin %<>% 
  na.omit()

# Retrieving marsupials and monotremata 
combine_imputed <- read_csv("data/combine_imputed.csv")
monotremata <- c("Ornithorhynchus anatinus", "Zaglossus attenboroughi", "Zaglossus bartoni", "Tachyglossus aculeatus")
marsupials_o <- c("Microbiotheria", "Dasyuromorphia", "Peramelemorphia", "Notoryctemorphia", "Diprotodontia", "Didelphimorphia", "Paucituberculata", "Australidelphia")
marsupials <- combine_imputed %>% filter(order %in% marsupials_o) %>% pull(iucn2020_binomial)

# Removing marsupials and monotremata 
NativeDataset_fin %<>% filter(Species1 %!in% monotremata & Species2 %!in% monotremata)
NativeDataset_fin %<>% filter(Species1 %!in% marsupials & Species2 %!in% marsupials)

# Only overlapping species:
NativeDataset_fin %<>% filter(Overlapping_Cells > 0)

#### Prediction dataset (Alien-native network) ####
alien_endemic_network <- read_csv("Data/alien_endemic_network.csv")

AlienDataset <- left_join(alien_endemic_network, species_dist_phylo_sharing, join_by(Native_Species == Species1, Alien_Species == Species2))
colSums(!is.na(AlienDataset)) / nrow(AlienDataset)

AlienDataset$n_virus[is.na(AlienDataset$n_virus)] <- 0

AlienDataset %<>% 
  na.omit()

AlienDataset %<>% 
  mutate(vir_bin = factor(ifelse(n_virus > 0, 1, 0))) %>% 
  mutate(Overlapping_Cells_log = log10(1 + AlienDataset$Overlapping_Cells)) 

AlienDataset <- left_join(AlienDataset, cit_virus %>% rename("n_cit1" = V1), join_by(Native_Species == X))
AlienDataset <- left_join(AlienDataset, cit_virus %>% rename("n_cit2" = V1), join_by(Alien_Species == X))

AlienDataset <- AlienDataset %>% filter(Overlapping_Cells_log > 0)

AlienDataset$log_min_cit <- log10(1 + pmin(AlienDataset$n_cit1, AlienDataset$n_cit2))
AlienDataset$log_max_cit <- log10(1 + pmax(AlienDataset$n_cit1, AlienDataset$n_cit2))
AlienDataset$log_sum_cit <- log10(1 + rowSums(AlienDataset[, c("n_cit1", "n_cit2")]))

# Removing marsupials and monotremata
AlienDataset %<>% filter(Alien_Species %!in% monotremata & Native_Species %!in% monotremata)
AlienDataset %<>% filter(Alien_Species %!in% marsupials & Native_Species %!in% marsupials)


AlienDataset %<>% na.omit() 
AlienDataset %<>% select(-c(3, 4, 5, 6))

columns_to_scale <- names(AlienDataset)[3:11]

scale_minmax_based_on_dataset <- function(column_name){
  # Extract the vector from the reference dataset
  max_val <- max(NativeDataset_fin[[column_name]], na.rm = TRUE)
  min_val <- min(NativeDataset_fin[[column_name]], na.rm = TRUE)
  
  # Return a function that will scale another vector based on those min/max
  function(x) {
    (x - min_val) / (max_val - min_val)
  }
}

AlienDataset[columns_to_scale] <- mapply(
  function(col, name) scale_minmax_based_on_dataset(name)(col),
  AlienDataset[columns_to_scale],
  columns_to_scale,
  SIMPLIFY = FALSE
)

## Saving datasets

AlienDataset %<>% rename("Species1" = Native_Species, "Species2" = Alien_Species) %>% select(-c(4, 5, 9, 10)) %>% 
  rename("phylo_sim" = divergence_time) %>% 
  mutate(phylo_sim = 1-phylo_sim) %>% 
  rename("foraging_sim" = GowDist_eco) %>% 
  mutate(foraging_sim = 1-foraging_sim) %>%  
  rename("trait_sim_eu" = EuDist_alltraits) %>% 
  mutate(trait_sim_eu = 1-trait_sim_eu) %>% 
  rename("trait_sim_gow" = GowDist_alltraits) %>% 
  mutate(trait_sim_gow = 1-trait_sim_gow) %>% 
  rename("trait_sim_mah" = MahDist_alltraits) %>% 
  mutate(trait_sim_mah = 1-trait_sim_mah)

write.csv(AlienDataset, "Data/DatasetPred.csv", row.names = F)

NativeDataset_fin[c(3:11)] <- lapply(NativeDataset_fin[c(3:11)], scale_minmax)

NativeDataset_fin %<>% select(-c(4, 5, 9, 10)) %>% 
  rename("phylo_sim" = divergence_time) %>% 
  mutate(phylo_sim = 1-phylo_sim) %>% 
  rename("foraging_sim" = GowDist_eco) %>% 
  mutate(foraging_sim = 1-foraging_sim) %>%  
  rename("trait_sim_eu" = EuDist_alltraits) %>% 
  mutate(trait_sim_eu = 1-trait_sim_eu) %>% 
  rename("trait_sim_gow" = GowDist_alltraits) %>% 
  mutate(trait_sim_gow = 1-trait_sim_gow) %>% 
  rename("trait_sim_mah" = MahDist_alltraits) %>% 
  mutate(trait_sim_mah = 1-trait_sim_mah)

write.csv(NativeDataset_fin, "Data/DatasetMain.csv", row.names = F)

