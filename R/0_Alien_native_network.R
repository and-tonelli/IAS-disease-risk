library(tidyverse)
library(terra)
library(Matrix)
require(tools)
require(purrr)
require(stringr)

# Heavy files at https://drive.google.com/drive/folders/1zqe7K-6qv-6VE_uvHOypio8qH-zts5D_
path_to_aoh <- "MAMMALS_5km_final"
path_to_alien_aoh <- "AOH_invasive/tiffs1975_5Km"

# Loading extent lists
load("ExtentList_native.RData")
ExtentList_native
load("ExtentList_alien.RData")
ExtentList_alien

set1_paths <- list.files(path_to_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE) # Native species
set2_paths <- list.files(path_to_alien_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE) # Alien species

# Read rasters
rasters1 <- lapply(set1_paths, rast)
rasters2 <- lapply(set2_paths, rast)

# Function to convert to alue matrix
to_sparse <- function(r) {
  v <- values(r, mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
}

sparse_list1 <- lapply(rasters1, to_sparse)
names(sparse_list1) <- sapply(rasters1, function(x) names(x))
names(sparse_list1) <- str_replace(names(sparse_list1), "ex_", "")

sparse_list2 <- lapply(rasters2, to_sparse)
names(sparse_list2) <- sapply(rasters2, function(x) names(x))

cell_counts1 <- sapply(sparse_list1, function(m) sum(m != 0))
cell_counts2 <- sapply(sparse_list2, function(m) sum(m != 0))

# All combinations
species_combinations <- expand.grid(names(sparse_list1), names(sparse_list2))
colnames(species_combinations) <- c("Native_sp", "Alien_sp")

# Overlap function
calculate_overlap <- function(species_pair) {
  sp_n <- species_pair[1] #native
  sp_a <- species_pair[2] #alien
  
  sparse1 <- sparse_list1[[sp_n]]
  sparse2 <- sparse_list2[[sp_a]]
  
  overlap <- sum((sparse1 != 0) & (sparse2 != 0))
  union_cells <- sum((sparse1 != 0) | (sparse2 != 0))
  
  result <- data.frame(
    Native_Species = sp_n,
    Alien_Species = sp_a,
    Cells_Native = cell_counts1[[sp_n]],
    Cells_Alien = cell_counts2[[sp_a]],
    Union_Cells = union_cells,
    Overlapping_Cells = overlap
  )
  
  return(result)
}

# Empty matric for all native-alien combinations
BlankMatrix <- matrix(NA, nrow = length(ExtentList_native), ncol = length(ExtentList_alien))
rownames(BlankMatrix) <- names(ExtentList_native)
colnames(BlankMatrix) <- names(ExtentList_alien)

# Change matrix values. If extents overlap -> 1, otherwise -> 0
for(i in 1:length(ExtentList_native)){
  
  for(j in 1:length(ExtentList_alien)){
    
    # ExtentList_alien[[j]] %>% print
    
    BlankMatrix[i, j] <- ifelse(is.null(intersect(ExtentList_native[[i]], ExtentList_alien[[j]])), 0, 1)
    
  }
  
}

# Add one row to "species_combinations" to flag if the extents overlap or not
species_combinations_copy <- cbind(species_combinations, ext_overlap = rep(NA, nrow(species_combinations)))

# Update species_combinations_copy based on the extent overlap matrix
for (n in seq(1:nrow(species_combinations))){
  n_sp <- as.character(species_combinations[n, 1])
  a_sp <- as.character(species_combinations[n, 2])
  
  if (BlankMatrix[n_sp, a_sp] == 1) {
    
    species_combinations_copy[n, 3] <- 1
    
  } else {
    
    species_combinations_copy[n, 3] <- 0
    
  }
}

# Get final combination table to loop over
filtered_combinations <- species_combinations_copy[species_combinations_copy[, 3] == 1, c(1, 2)]
print(nrow(filtered_combinations))

print(paste0("Number of combinations after filtering: ", nrow(filtered_combinations)))


# Applying function to the final combination table
overlap_df <- NULL
for (k in 1:nrow(filtered_combinations)){
  
  one_overlap <- calculate_overlap(filtered_combinations[k, ])
  overlap_df <- rbind(overlap_df, one_overlap)
  
  print(paste0(k, "/", nrow(filtered_combinations)))  
  write.csv(overlap_df, "Data/alien_endemic_network1975.csv", row.names = F) # saving csv at each iteration
  
}

