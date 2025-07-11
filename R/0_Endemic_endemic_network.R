# Associate mammal species with the countries where they have been found
require(rnaturalearth)
require(tidyverse)
require(sf)
require(terra)
library(exactextractr)
library(purrr)

raster_files <- list.files("MAMMALS_5km_sum", pattern = ".tif", full.names = TRUE, recursive = TRUE)

res_deg <- 5 / 110.574  
blank_r <- rast(
  xmin = -180, xmax = 180,
  ymin = -90, ymax = 90,
  resolution = res_deg,
  crs = "EPSG:4326",
  vals = NA
)

Mammals_GADMcountries <- read_csv("Mammals_GADMcountries_final.csv")

endemic_species <- unique(Mammals_GADMcountries$Binomial)

# Load rasters and country polygons
sf_use_s2(FALSE)

tab_overlaps <- NULL

for (e_sp1 in endemic_species){
  
  setwd("/data/")
  path_endemic1 <- raster_files[grep(str_replace(e_sp1, " ", "_"), raster_files)]
  endemic_rast1 <- rast(paste0(path_endemic1))
  
  if (minmax(endemic_rast1)[2] < 2500){
    
    endemic_rast1[endemic_rast1[] < minmax(endemic_rast1)[2]] <- NA
    endemic_rast1[endemic_rast1[] == minmax(endemic_rast1)[2]] <- 1
    
  }
  
  else {
    
    endemic_rast1[endemic_rast1[] < 2500] <- NA
    endemic_rast1[endemic_rast1[] >= 2500] <- 1
  
  }
  
  endemic_countries <- unique(Mammals_GADMcountries %>% 
                              filter(Binomial == e_sp1) %>% 
                              pull(GID_0))
  
  country_mammals <- Mammals_GADMcountries %>% filter(GID_0 %in% endemic_countries) %>% pull(Binomial)
  # country_mammals <- str_replace(country_mammals, "_", " ")
    
    for (e_sp2 in country_mammals){
      
      setwd("/data/")
      path_endemic2 <- raster_files[grep(e_sp2, raster_files)]
      endemic_rast2 <- rast(paste0(path_endemic2))
      endemic_countries2 <- Mammals_GADMcountries %>% filter(GID_0 %in% endemic_countries, Binomial == e_sp2) %>% pull(GID_0)
      
      if (minmax(endemic_rast2)[2] < 2500){
        
        endemic_rast2[endemic_rast2[] < minmax(endemic_rast2)[2]] <- NA
        endemic_rast2[endemic_rast2[] == minmax(endemic_rast2)[2]] <- 1
        
      }
      
      else {
        
        endemic_rast2[endemic_rast2[] < 2500] <- NA
        endemic_rast2[endemic_rast2[] >= 2500] <- 1
        
      }
      
      overlap <- terra::mask(endemic_rast1, endemic_rast2) 
      union <- endemic_rast1+endemic_rast2
      
      ncell <- length(overlap[!is.na(overlap)])
      ncell_tot <- length(union[!is.na(union)])
      
      if (ncell > 0){
        
        tab_overlaps <- rbind(tab_overlaps, tibble(species1 = e_sp1,
                                                   species2 = e_sp2,
                                                   tot_cells = ncell_tot,
                                                   overlap_cells = ncell,
                                                   overlap_prop = overlap_cells/tot_cells,
                                                   countries_overlap = intersect(endemic_countries, endemic_countries2) |> paste(collapse = ", "),
                                                   ))
        
        write.csv(tab_overlaps, "endemic_endemic_overlap.csv", row.names = F)
        
      }
      
    }
    
  }

