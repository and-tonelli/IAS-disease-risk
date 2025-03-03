library(terra)
library(Matrix)
library(dplyr)

# heavy files at https://drive.google.com/drive/folders/1zqe7K-6qv-6VE_uvHOypio8qH-zts5D_
path_to_aoh <- "MAMMALS_5km_final"
path_to_alien_aoh <- "AOH_invasive/tiffs1975_5Km"

set1_paths <- list.files(path_to_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)  # Native species
set2_paths <- list.files(path_to_alien_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)  # Alien species

species_names1 <- gsub(".tif$", "", basename(set1_paths))
species_names2 <- gsub(".tif$", "", basename(set2_paths))

rasters1 <- lapply(set1_paths, rast)
rasters2 <- lapply(set2_paths, rast)

# Convert to sparse matrix
to_sparse <- function(r) {
  v <- values(r, mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
}

cell_counts1 <- sapply(sparse_list1, function(m) sum(m != 0))
cell_counts2 <- sapply(sparse_list2, function(m) sum(m != 0))

# All combinatiosn
index_combinations <- expand.grid(seq_along(species_names1), seq_along(species_names2))

# Calculate overlap function
calculate_overlap <- function(index_pair) {
  i <- as.numeric(index_pair[1])  
  j <- as.numeric(index_pair[2])
  
  sparse1 <- to_sparse(rasters1[[i]])
  sparse2 <- to_sparse(rasters2[[j]])
  
  overlap <- sum((sparse1 != 0) & (sparse2 != 0))
  union_cells <- sum((sparse1 != 0) | (sparse2 != 0))
  
  result <- data.frame(
    Native_Sp = species_names1[i],
    Alien_Sp = species_names2[j],
    Cells_Native = sum(sparse1 != 0),
    Cells_Alien = sum(sparse2 != 0),
    Union_Cells = union_cells,
    Overlap_Cells = overlap
  )
  
  return(result)
}

options(future.globals.maxSize = 5 * 1024^3)
overlap_list <- list()

# Parallel (not working - so setting workers = 1 for now)
plan(multisession, workers = 1) 

overlap_list <- future_lapply(seq_len(nrow(index_combinations)), function(k) {
  index_pair <- index_combinations[k, ] # go throughall combinations
  result <- calculate_overlap(index_pair) 
  
  return(result)
}, future.seed = TRUE)

overlap_df <- do.call(rbind, overlap_list)

write.csv(overlap_df, "Data/alien_native_network1975.csv", row.names = F)

