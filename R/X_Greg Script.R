
# X_Greg Script.R ####

library(tidyverse); library(fs); library(igraph);
library(ggraph); library(tidygraph);
library(magrittr)

DataList <- "Data" %>% dir_ls(regex = "csv") %>% map(read.csv)

DataList %>% map(head)

DataList %>% 
  last %>% select(Host = spp, Country = countries) %>% 
  table() %>% 
  graph_from_incidence_matrix() %>% 
  plot

BipartiteGraph <- 
  DataList$`Data/Mammals_GADMcountries.csv` %>% 
  # last %>% 
  select(Host = species, Country = countries) %>% 
  mutate_at("Country", ~str_split(.x, ", ")) %>% 
  unnest(Country) %>% 
  table() %>% 
  graph_from_incidence_matrix()

Unipartite <- 
  BipartiteGraph %>% 
  bipartite.projection() %>% 
  extract2("proj1") %>% 
  as_tbl_graph

giant.component <- function(graph, ...) {
  
  cl <- clusters(graph, ...)
  
  subgraph(graph, which(cl$membership == which.max(cl$csize)-1)-1)
  
}

giant.component <- function(g, ...) {
  
  graphs <- decompose.graph(g)
  largest <- which.max(sapply(graphs, vcount))
  
  return(graphs[[largest]])
  
}

library(cowplot)

theme_set(theme_cowplot())

Unipartite %>% as_tbl_graph() %>% activate(nodes) %>% 
  slice(1:1000) %>%
  ggraph() + 
  geom_node_point() +
  coord_fixed()

SF <- read_sf("Data/DAMA Shapefiles")

SF

Invasives <- 
  "Data/DAMA Shapefiles" %>% 
  list.files %>% str_split("[.]") %>% 
  map_chr(1) %>% unique

FocalSpecies <- Invasives[1]

dir_create("Data/DAMA Individual Species")

for(FocalSpecies in Invasives){
  
  print(FocalSpecies)
  
  dir_create(paste0("Data/DAMA Individual Species/", FocalSpecies))
  
  file.copy(paste0("Data/DAMA Shapefiles/", FocalSpecies, c(".dbf", ".prj", ".shp", ".shx")),
            paste0("Data/DAMA Individual Species/", FocalSpecies, "/", FocalSpecies, c(".dbf", ".prj", ".shp", ".shx")))
  
}

ShapefileList <- list()

for(FocalSpecies in Invasives){
  
  print(FocalSpecies)
  
  Shapefile <- read_sf(paste0("Data/DAMA Individual Species/", FocalSpecies))
  
  ShapefileList[[FocalSpecies]] <- Shapefile
  
}

ShapeFile <- 
  ShapefileList %>% 
  bind_rows()

ShapeFile[2,] %>% plot


library(rnaturalearth)
library(rnaturalearthdata)

World <- ne_countries(scale = "medium", returnclass = "sf")

World %<>% filter(!continent == "Antarctica")

ShapeFile[1:25,] %>% 
  ggplot() + 
  # geom_sf(data = World) %>%
  geom_sf(colour = "grey")

ShapeFile[1:250 + 250,] %>% 
  ggplot() + 
  # geom_sf(data = World) %>%
  geom_sf(#colour = "grey", 
    aes(colour = Binomial)) +
  facet_wrap(~Binomial)

ShapeFile[1:250,] %>% 
  ggplot() + 
  # geom_sf(data = World) %>%
  geom_sf(#colour = "grey", 
    aes(colour = Binomial)) +
  facet_wrap(~Binomial, scales = "free")

ggsave("AllRanges.jpeg")

