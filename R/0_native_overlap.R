#### Native-native overlap ####
library(terra)
library(Matrix)
library(dplyr)
require(tools)
require(future.apply)
require(purrr)
require(tools)
require(stringr)

path_to_aoh <- "..."


set1_paths <- list.files(path_to_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE) # Native species
set2_paths <- list.files(path_to_aoh, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE) # Native species 

# Read rasters
rasters1 <- lapply(set1_paths, rast)
rasters2 <- lapply(set2_paths, rast)

# Function to convert to sparse matrix
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
names(sparse_list2) <- str_replace(names(sparse_list2), "ex_", "")


cell_counts1 <- sapply(sparse_list1, function(m) sum(m != 0))
cell_counts2 <- sapply(sparse_list2, function(m) sum(m != 0))

# All combinations
species_combinations <- expand.grid(names(sparse_list1), names(sparse_list2))
colnames(species_combinations) <- c("Species1", "Species2")

# Overlap function
calculate_overlap <- function(species_pair) {
  sp_n <- as.character(species_pair$Species1) #native
  sp_a <- as.character(species_pair$Species2) #alien
  
  # sparse1 <- to_sparse(rasters1[[i]])
  # sparse2 <- to_sparse(rasters2[[j]])
  
  sparse1 <- sparse_list1[[sp_n]]
  sparse2 <- sparse_list2[[sp_a]]
  
  overlap_mask <- (sparse1 != 0) & (sparse2 != 0)
  
  overlap <- sum(overlap_mask)
  
  union_cells <- sum((sparse1 != 0) | (sparse2 != 0))
  
  if (overlap > 0) {
    # Get cell numbers of overlapping cells
    overlap_indices <- which(overlap_mask)
    
    # Use raster1 as reference
    ref_raster <- rasters1[[which(names(sparse_list1) == sp_n)]]
    
    # Convert cell indices to spatial coordinates
    xy <- terra::xyFromCell(ref_raster, overlap_indices)
    
    centroid_x <- mean(xy[, 1])
    centroid_y <- mean(xy[, 2])
  } else {
    centroid_x <- NA
    centroid_y <- NA
  }
  
  result <- data.frame(
    Native_Species1 = sp_n,
    Native_Species2 = sp_a,
    Cells_Species1 = cell_counts1[[sp_n]],
    Cells_Species2 = cell_counts2[[sp_a]],
    Union_Cells = union_cells,
    Overlapping_Cells = overlap,
    Centroid_X = centroid_x,
    Centroid_Y = centroid_y
    
  )
  
  return(result)
}

# Loading extent lists (saved as vectors)
ExtentList_native <- readRDS("PRJ_IAS/data/ExtentList_native.rds")
ExtentList_native <- lapply(ExtentList_native, function(coords) {
  ext(coords[1], coords[2], coords[3], coords[4]) # needed to transform them back to SpatExtent
})

ExtentList_alien <- readRDS("PRJ_IAS/data/ExtentList_native.rds") 
ExtentList_alien <- ExtentList_alien[names(ExtentList_alien) %in% names(sparse_list2)]
ExtentList_alien <- lapply(ExtentList_alien, function(coords) {
  ext(coords[1], coords[2], coords[3], coords[4])
})

BlankMatrix <- matrix(NA, nrow = length(ExtentList_native), ncol = length(ExtentList_alien))
rownames(BlankMatrix) <- names(ExtentList_native)
colnames(BlankMatrix) <- names(ExtentList_alien)

for(i in 1:length(ExtentList_native)){
  
  for(j in 1:length(ExtentList_alien)){
    
    BlankMatrix[i, j] <- ifelse(is.null(terra::intersect(ExtentList_native[[i]], ExtentList_alien[[j]])), 0, 1)
    
  }
  
}

species_combinations_copy <- cbind(species_combinations, ext_overlap = rep(NA, nrow(species_combinations)))

for (n in seq(1:nrow(species_combinations))){
  n_sp <- as.character(species_combinations[n, 1])
  a_sp <- as.character(species_combinations[n, 2])
  
  if (BlankMatrix[n_sp, a_sp] == 1) {
    
    species_combinations_copy[n, 3] <- 1
    
  } else {
    
    species_combinations_copy[n, 3] <- 0
    
  }
}

filtered_combinations <- species_combinations_copy[species_combinations_copy[, 3] == 1, c(1, 2)]
print(nrow(filtered_combinations))

print(paste0("Finished filtering index combinations: ", nrow(filtered_combinations)))


overlap_df <- NULL
for (k in 1:nrow(filtered_combinations)){
  
  one_overlap <- calculate_overlap(filtered_combinations[k, ])
  overlap_df <- rbind(overlap_df, one_overlap)
  print(paste0(k, "/", nrow(filtered_combinations)))  
  write.csv(overlap_df, "data/native_native_overlap3.csv", row.names = F)
  
}

