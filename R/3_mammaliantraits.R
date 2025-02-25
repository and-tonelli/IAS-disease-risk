# Retrieving mammalian traits and computing pairwaise distances

require(tidyverse)
require(usdm)
require(magrittr)
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%

combine_imputed <- read_csv("data/combine_imputed.csv")
vector_of_terrestrialAOH <- readRDS("data/vector_of_terrestrialAOH.rds")

# filtering out aquatic species
combine_imputed %>% filter((marine == 1 | freshwater == 1) & iucn2020_binomial %!in% vector_of_terrestrialAOH) %>% pull(iucn2020_binomial) -> aquatic

key_traits <- c("adult_mass_g", "max_longevity_d", "age_first_reproduction_d", "gestation_length_d", "litter_size_n", "litters_per_year_n", "interbirth_interval_d", "weaning_age_d", "generation_length_d")

combine_imputed %>% 
  filter(iucn2020_binomial %!in% aquatic) %>% 
  select(c(1, 2, 5, 6, key_traits)) %>% 
  na.omit() %>% 
  data.frame() -> filtered_combine

#### Starting with euclidean distance ####
# Remove highly collinear variables
vif <- usdm::vifstep(filtered_combine[-c(1:4)], th = 3)

vif@excluded

# Excluding: "generation_length_d"      "interbirth_interval_d"    "age_first_reproduction_d" "max_longevity_d" 
filtered_combine_eu <- filtered_combine %>% select(-vif@excluded)

# Rescaling variables
scale_minmax <- function(column){
  
  max_val <- max(column)
  min_val <- min(column)
    
  return((column - min_val)/(max_val - min_val))
  
}

filtered_combine_eu2 <- filtered_combine_eu
filtered_combine_eu2[5:9] <- lapply(filtered_combine_eu2[5:9], scale_minmax)

# Actually computing euclidean distance
eu_dist_matrix <- stats::dist(filtered_combine_eu2[, c(5:9)], method = "euclidean", diag = T) %>% as.matrix()

rownames(eu_dist_matrix) <- filtered_combine_eu2$iucn2020_binomial
colnames(eu_dist_matrix) <- filtered_combine_eu2$iucn2020_binomial

eu_dist_matrix["Elephas maximus", "Loxodonta africana"]

#### Trying Mahalanobis distance ####
require(StatMatch)
mah_dist_matrix <- mahalanobis.dist(filtered_combine[5:13])

rownames(mah_dist_matrix) <- filtered_combine$iucn2020_binomial
colnames(mah_dist_matrix) <- filtered_combine$iucn2020_binomial
mah_dist_matrix["Pteropus giganteus", "Pteropus vampyrus"]
mah_dist_matrix["Elephas maximus", "Loxodonta africana"]
# eu_dist_matrix["Pteropus giganteus", "Pteropus vampyrus"]

eu_dist_matrix[upper.tri(eu_dist_matrix)] <- NA
mah_dist_matrix[upper.tri(mah_dist_matrix)] <- NA
eu_dist_matrix[1:5, 1:5]

min(eu_dist_matrix, na.rm = TRUE)
max(eu_dist_matrix, na.rm = TRUE)

eu_pos <- which(eu_dist_matrix != 0 & !is.na(eu_dist_matrix), arr.ind = TRUE)
eu_values <- cbind(eu_pos, value = eu_dist_matrix[eu_pos])

mah_pos <- which(mah_dist_matrix != 0 & !is.na(mah_dist_matrix), arr.ind = TRUE)
mah_values <- cbind(mah_pos, value = mah_dist_matrix[mah_pos])

colnames(eu_values) <- c("row", "col", "value_eu")
colnames(mah_values) <- c("row", "col", "value_mah")

cbind(eu_values[1:500, ], mah_values[1:500, 3]) %>% 
  data.frame() %>% 
  ggplot(aes(x = value_eu, y = V4), size = 3)+
  geom_point()+
  labs(x = "Euclidean distance", y = "Mahalanobis distance")+
  theme_bw()

mah_dist_matrix[404, 1]

min(mah_dist_matrix, na.rm = TRUE)
max(mah_dist_matrix, na.rm = TRUE)

rownames(mah_dist_matrix)[404]
colnames(mah_dist_matrix)[1]

#### Gower's distance for categorical variables ####
require(cluster)
combine_filtered_gow <- combine_imputed %>% 
  filter(iucn2020_binomial %!in% aquatic) %>% 
  # filter(iucn2020_binomial %in% colnames(eu_dist_matrix)) %>% 
  select(c(1, 2, 5, 6, "trophic_level", "foraging_stratum")) %>% 
  na.omit() %>% 
  mutate(trophic_level = factor(trophic_level, levels = 1:3, labels = c("herbivore", "omnivore", "carnivore"), ordered = T),
         foraging_stratum = factor(foraging_stratum) %>% ordered(c("M", "G", "S", "Ar", "A")))


gow_dist_matrix <- daisy(combine_filtered_gow[, c(5, 6)], metric = "gower") %>% as.matrix()
rownames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial
colnames(gow_dist_matrix) <- combine_filtered_gow$iucn2020_binomial

gow_dist_matrix["Abrocoma bennettii", "Abditomys latidens"]
min(gow_dist_matrix, na.rm = TRUE)
max(gow_dist_matrix, na.rm = TRUE)

#### Mahalanobis distance sim ####
require(StatMatch)
tib_sim <- tibble(x = c(1, 1.5, 1.5, 2, 2.3, 3, 3, 3, 3.4, 4, 6.2, 6.5, 4.2, 4.9, 5, 4.5),
                  y = c(1.8, 3, 3.3, 3.8, 5, 6, 5.8, 5.7, 7, 7.6, 11.5, 13, 13, 14, 10, 9.3))

tib_sim %>% 
  ggplot()+
  geom_point(aes(x = x, y = y), size = 3, color = "steelblue4")+
  theme_bw()

dist_matrix_sim <- mahalanobis.dist(tib_sim[c(1, 2)])
dist_matrix_sim_eu <- stats::dist(tib_sim[c(1, 2)], "euclidean", diag = T)

as.matrix(dist_matrix_sim_eu)[13, 14]
as.matrix(dist_matrix_sim)[13, 14]

as.matrix(dist_matrix_sim_eu)[5, 6]
as.matrix(dist_matrix_sim)[5, 6]

