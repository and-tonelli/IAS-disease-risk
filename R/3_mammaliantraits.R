# Retrieving mammalian traits and computing pairwaise distances

require(tidyverse)
require(usdm)
require(magrittr)
require(readxl)
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%
theme_set(theme_bw())

mammals_with_AOH <- readRDS("Data/AOH_terrestrialspecies.rds")

trait_data_imputed <- read.csv("Data/trait_data_imputed.csv") %>% 
  filter(iucn2020_binomial %in% mammals_with_AOH)

trait_data_reported <- read.csv("Data/trait_data_reported.csv") %>% 
  filter(iucn2020_binomial %in% mammals_with_AOH)

completeness_imp <- colSums(!is.na(trait_data_imputed)) / nrow(trait_data_imputed)
completeness_rep <- colSums(!is.na(trait_data_reported)) / nrow(trait_data_reported)

completeness_summary_imp <- data.frame(Var = names(completeness_imp),
                                       Completeness_imputed = completeness_imp)

completeness_summary_rep <- data.frame(Var = names(completeness_rep),
                                       Completeness_reported = completeness_rep)

Summary_comp <- cbind(completeness_summary_imp, completeness_summary_rep[2])
Summary_comp %>% filter(Completeness_imputed >= 0.90) %>% arrange(desc(Completeness_reported))


key_traits <- c("adult_mass_g", "max_longevity_d", "age_first_reproduction_d", "gestation_length_d", "litter_size_n", "litters_per_year_n", "weaning_age_d")

#### Euclidean distance ####
# Remove highly collinear variables
vif <- usdm::vifstep(trait_data_imputed[key_traits], th = 3)

vif@excluded

# Excluding: "age_first_reproduction_d"
filtered_combine_eu <- trait_data_imputed[c("order", "family", "iucn2020_binomial", "phylacine_binomial", key_traits)] %>%
  select(-vif@excluded) %>% 
  na.omit() # losing 181 species

# Rescaling variables
scale_minmax <- function(column){
  
  max_val <- max(column)
  min_val <- min(column)
  
  return((column - min_val)/(max_val - min_val))
  
}

filtered_combine_eu[5:10] <- lapply(filtered_combine_eu[5:10], scale_minmax)

# Computing euclidean distance
eu_dist_matrix <- stats::dist(filtered_combine_eu[, c(5:10)], method = "euclidean", diag = T) %>% as.matrix()

rownames(eu_dist_matrix) <- filtered_combine_eu$iucn2020_binomial
colnames(eu_dist_matrix) <- filtered_combine_eu$iucn2020_binomial

eu_dist_matrix["Elephas maximus", "Loxodonta africana"]
eu_dist_matrix["Elephas maximus", "Rattus norvegicus"]

#### Mahalanobis distance ####
require(StatMatch)
filtered_combine_ma <- trait_data_imputed[c("order", "family", "iucn2020_binomial", "phylacine_binomial", key_traits)] %>% 
  na.omit() # losing 183 species

mah_dist_matrix <- mahalanobis.dist(filtered_combine_ma[5:11])

rownames(mah_dist_matrix) <- filtered_combine_ma$iucn2020_binomial
colnames(mah_dist_matrix) <- filtered_combine_ma$iucn2020_binomial
mah_dist_matrix["Elephas maximus", "Loxodonta africana"]
mah_dist_matrix["Elephas maximus", "Rattus norvegicus"]
# mah_dist_matrix["Pteropus giganteus", "Pteropus vampyrus"]

eu_dist_matrix[upper.tri(eu_dist_matrix)] <- NA
mah_dist_matrix[upper.tri(mah_dist_matrix)] <- NA

min(eu_dist_matrix, na.rm = TRUE)
max(eu_dist_matrix, na.rm = TRUE)
min(mah_dist_matrix, na.rm = TRUE)
max(mah_dist_matrix, na.rm = TRUE)

eu_pos <- which(eu_dist_matrix != 0 & !is.na(eu_dist_matrix), arr.ind = TRUE)
eu_values <- cbind(eu_pos, value = eu_dist_matrix[eu_pos])

mah_pos <- which(mah_dist_matrix != 0 & !is.na(mah_dist_matrix), arr.ind = TRUE)
mah_values <- cbind(mah_pos, value = mah_dist_matrix[mah_pos])

colnames(eu_values) <- c("row", "col", "value_eu")
colnames(mah_values) <- c("row", "col", "value_mah")

# Euclidean distance distribution
mah_values %>%
  ggplot()+
  geom_histogram(aes(x = value_mah), fill = "steelblue4")


# Mahalanobis distance distribution
eu_values %>% 
  ggplot()+
  geom_histogram(aes(x = value_eu), fill = "steelblue4")


# mahalanobis ~ euclidean (first 500)
cbind(eu_values[1:500, ], mah_values[1:500, 3]) %>% 
  data.frame() %>% 
  ggplot(aes(x = value_eu, y = V4), size = 3)+
  geom_point()+
  labs(x = "Euclidean distance", y = "Mahalanobis distance")+
  theme_bw()

#### Gower's distance for categorical variables ####
require(cluster)
combine_filtered_gow <- trait_data_imputed %>% 
  select(c("order", "family", "iucn2020_binomial", "phylacine_binomial", "trophic_level", "foraging_stratum", 
           "dphy_invertebrate", "dphy_vertebrate", "dphy_plant", all_of(key_traits)
           # "det_inv", "det_vend" ,"det_vect", "det_vfish", "det_vunk", "det_scav", "det_fruit", "det_nect", "det_seed", "det_plantother" # too many missing data (losing ~1500 spp)
  )) %>% 
  na.omit() %>% 
  mutate(trophic_level = factor(trophic_level, levels = 1:3, labels = c("herbivore", "omnivore", "carnivore"), ordered = T),
         foraging_stratum = factor(foraging_stratum) %>% ordered(c("M", "G", "S", "Ar", "A")))


# Just diet (diet%, trophic level, foraging stratum)
gow_dist_matrix <- daisy(combine_filtered_gow[, c(5:9)], metric = "gower") %>% as.matrix()
rownames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial
colnames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial

gow_dist_matrix["Elephas maximus", "Loxodonta africana"]
gow_dist_matrix["Elephas maximus", "Rattus norvegicus"]
gow_dist_matrix["Canis lupus", "Vulpes lagopus"]

min(gow_dist_matrix, na.rm = TRUE)
max(gow_dist_matrix, na.rm = TRUE)

# Life history traits (+ bm)
gow_dist_matrix2 <- daisy(combine_filtered_gow[, c(10:16)], metric = "gower") %>% as.matrix()
rownames(gow_dist_matrix2) <- combine_filtered_gow$iucn2020_binomial
colnames(gow_dist_matrix2) <- combine_filtered_gow$iucn2020_binomial

gow_dist_matrix2["Elephas maximus", "Loxodonta africana"]
gow_dist_matrix2["Elephas maximus", "Rattus norvegicus"]
gow_dist_matrix2["Canis lupus", "Vulpes lagopus"]

min(gow_dist_matrix2, na.rm = TRUE)
max(gow_dist_matrix2, na.rm = TRUE)

