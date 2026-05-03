# Retrieving mammalian traits and computing pairwise distances
require(tidyverse)
require(usdm)
require(magrittr)
require(readxl)
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%
theme_set(theme_bw())

# heavy files at https://drive.google.com/drive/folders/1zqe7K-6qv-6VE_uvHOypio8qH-zts5D_

mammals_with_AOH <- readRDS("Data/AOH_terrestrialspecies.rds")

trait_data_imputed <- read.csv("Data/trait_data_imputed.csv") %>% 
  filter(iucn2020_binomial %in% mammals_with_AOH)

trait_data_reported <- read.csv("Data/trait_data_reported.csv") %>% 
  filter(iucn2020_binomial %in% mammals_with_AOH)

# Checking completeness of available traits
completeness_imp <- colSums(!is.na(trait_data_imputed)) / nrow(trait_data_imputed)
completeness_rep <- colSums(!is.na(trait_data_reported)) / nrow(trait_data_reported)

completeness_summary_imp <- data.frame(Var = names(completeness_imp),
                                       Completeness_imputed = completeness_imp)

completeness_summary_rep <- data.frame(Var = names(completeness_rep),
                                       Completeness_reported = completeness_rep)

Summary_comp <- cbind(completeness_summary_imp, completeness_summary_rep[2])
Summary_comp %>% filter(Completeness_imputed >= 0.90) %>% arrange(desc(Completeness_reported))

# Selecting traits with >90% completeness and <70% imputation
key_traits <- c("adult_mass_g", "max_longevity_d", "age_first_reproduction_d", "gestation_length_d", "litter_size_n", "litters_per_year_n", "weaning_age_d")

#### Gower's distance for continuous and categorical variables ####
require(cluster)
combine_filtered_gow <- trait_data_imputed %>% 
  select(c("order", "family", "iucn2020_binomial", "phylacine_binomial", "trophic_level", "foraging_stratum", 
           "dphy_invertebrate", "dphy_vertebrate", "dphy_plant", all_of(key_traits)
           # "det_inv", "det_vend" ,"det_vect", "det_vfish", "det_vunk", "det_scav", "det_fruit", "det_nect", "det_seed", "det_plantother" # too many missing data (losing ~1500 spp)
  )) %>% 
  na.omit() %>% 
  mutate(trophic_level = factor(trophic_level, levels = 1:3, labels = c("herbivore", "omnivore", "carnivore"), ordered = T),
         foraging_stratum = factor(foraging_stratum) %>% ordered(c("M", "G", "S", "Ar", "A")))


# Just diet (diet%, trophic level, foraging stratum) --> Foraging distance
gow_dist_matrix <- daisy(combine_filtered_gow[, c(5:9)], metric = "gower") %>% as.matrix()
rownames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial
colnames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial

# check 
min(gow_dist_matrix, na.rm = TRUE)
max(gow_dist_matrix, na.rm = TRUE)

# Life history traits (+ bm) --> Traits distance
gow_dist_matrix2 <- daisy(combine_filtered_gow[, c(10:16)], metric = "gower") %>% as.matrix()
rownames(gow_dist_matrix2) <- combine_filtered_gow$iucn2020_binomial
colnames(gow_dist_matrix2) <- combine_filtered_gow$iucn2020_binomial

# check 
gow_dist_matrix2["Elephas maximus", "Rattus rattus"]
gow_dist_matrix2["Canis lupus", "Vulpes vulpes"]

min(gow_dist_matrix2, na.rm = TRUE)
max(gow_dist_matrix2, na.rm = TRUE)

#### Phylogenetic distance obtain from a mammalian phylogenetic tree extracted from PHYLACINE (Faurby et al., 2018) ####
phylo_dist <- read.csv("data/PHYLACINE_pairwaise_div_time_mammals.csv") %>% select(-3)
