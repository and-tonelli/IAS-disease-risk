# Associate mammal species with the countries where they have been found
require(rnaturalearth)
require(tidyverse)
require(sf)
require(terra)
library(exactextractr)
library(purrr)
library(Matrix)
require(stringr)

path_to_aoh <- "MAMMALS_5km_final"

# Take all rasters
r_list <- list.files(path_to_aoh, pattern = ".tif", full.names = TRUE, recursive = TRUE)
n <- length(r_list)

# Read rasters
rasters <- lapply(r_list, rast)
print("Finished loading rasters")

# Rasters to matrices (1s and 0s -- NAs as 0s)
to_sparse <- function(r) {
  v <- values(r, mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
}

options(future.globals.maxSize = 5 * 1024^3)

sparse_list <- lapply(rasters, to_sparse)
names(sparse_list) <- sapply(rasters, function(x) names(x))
names(sparse_list) <- str_replace(names(sparse_list), "ex_", "")

cell_counts <- sapply(sparse_list, function(m) sum(m != 0))

print("Finished converting rasters to matrices")

# Function to calculate overlap between two matrices
calculate_overlap <- function(species_pair) {
  sp_1 <- species_pair[1]
  sp_2 <- species_pair[2]
  
  sparse1 <- sparse_list[[sp_1]]
  sparse2 <- sparse_list[[sp_2]]
  
  overlap <- sum((sparse1 != 0) & (sparse2 != 0))
  union_cells <- sum((sparse1 != 0) | (sparse2 != 0))
  
  if (overlap > 0) {
    return(data.frame(
      Species1 = sp_1,
      Species2 = sp_2,
      Cells_Species1 = cell_counts[[sp_1]],
      Cells_Species2 = cell_counts[[sp_2]],
      Union_Cells = union_cells,
      Overlapping_Cells = overlap
    ))
  } else {
    return(NULL)
  }
}

# All unique species combinations
species_combinations <- t(combn(names(sparse_list), 2))

# Filter index combinations that do not overlap (non-overlapping extents)
# Loading extent lists (saved as vectors)
ExtentList_native <- readRDS("ExtentList_native.rds")
ExtentList_native <- lapply(ExtentList_native, function(coords) {
  ext(coords[1], coords[2], coords[3], coords[4]) # needed to transform them back to SpatExtent
})

BlankMatrix <- matrix(NA, nrow = length(ExtentList_native), ncol = length(ExtentList_native))
rownames(BlankMatrix) <- names(ExtentList_native)
colnames(BlankMatrix) <- names(ExtentList_native)

for(i in 1:length(ExtentList_native)){
  
  # ExtentList_native[[i]] %>% print
  
  for(j in i:length(ExtentList_native)){
    
    BlankMatrix[i, j] <- ifelse(is.null(terra::intersect(ExtentList_native[[i]], ExtentList_native[[j]])), 0, 1)
    
  }
  
}

species_combinations_copy <- cbind(species_combinations, rep(NA, nrow(species_combinations)))

for (n in seq(1:nrow(species_combinations))){
  sp_1 <- as.character(species_combinations[n, 1])
  sp_2 <- as.character(species_combinations[n, 2])
  
  if (BlankMatrix[sp_1, sp_2] == 1) {
    
    species_combinations_copy[n, 3] <- 1
    
  } else {
    
    species_combinations_copy[n, 3] <- 0
    
  }
}

filtered_combinations <- species_combinations_copy[species_combinations_copy[, 3] == 1, c(1, 2)]
print(paste0("Finished filtering index combinations: ", nrow(filtered_combinations)))

# Combine results into one dataframe
overlap_df <- NULL
for (k in 1:nrow(filtered_combinations)){
  
  one_overlap <- calculate_overlap(filtered_combinations[k, ])
  overlap_df <- rbind(overlap_df, one_overlap)
  print(paste0(k, "/", nrow(filtered_combinations)))  
  write.csv(overlap_df, "Data/endemic_endemic_network.csv", row.names = F)
  
}
