
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

# Demonstrating filtering ####

VIRION_DATA %>% 
  graph_from_data_frame(directed = F) %>% 
  as_tbl_graph %>% 
  activate(edges) %>% 
  filter(Detection_Antibodies == 1)

VIRION_DATA[c(1, 10)] %>% 
  graph_from_data_frame(directed = F) %>% 
  as_tbl_graph %>% 
  activate(edges) %>% 
  # mutate(EdgeWeight = 1) %>% dplyr::select(EdgeWeight) %>% 
  group_by(to) %>%
  simplify(edge.attr.comb = "sum") %>% 
  filter(countries == "CHN")

plot_data <- VIRION_DATA[1:1000, c(1, 10)] %>%
  rename(from = Harmonised_host, to = VirusFamily) %>%
  as.data.frame()

# Create the graph
graph <- graph_from_data_frame(plot_data, directed = FALSE) %>% 
  as_tbl_graph()

BipartiteGraph <- 
  VIRION_DATA %>%
  filter(HostClass == "mammalia")
select(Host = Harmonised_host, VirusFamily) %>% 
  table() %>% 
  graph_from_incidence_matrix()

# Unipartite graph species-species (edge is sharing)
BipartiteGraph %>% 
  bipartite.projection() %>% 
  extract2("proj1") %>% 
  as_tbl_graph %>% 
  ggraph(layout = "kk") + 
  geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) + 
  geom_node_point(size = 1, color = "firebrick") + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3) + 
  theme_void()

# Graphs of target viral families
vir_tab <- read_csv("Data/vir_tab_long_final.csv")

alien_viruses <- vir_tab %>% pull(virus) %>% unique
target_fams <- VIRION_DATA %>% filter(Virus %in% alien_viruses) %>% pull(VirusFamily) %>% unique
target_genera <- VIRION_DATA %>% filter(Virus %in% alien_viruses & !is.na(VirusGenus)) %>% pull(VirusGenus) %>% unique

VIRION_DATA <- read.csv("Data/Virion_mammals.csv")

VIRION_DATA %>%
  filter(HostClass == "mammalia", VirusGenus %in% target_genera) %>% 
  select(Host, VirusGenus) %>% 
  table() %>% 
  graph_from_incidence_matrix() %>% 
  bipartite.projection() %>% 
  extract2("proj1") %>% 
  as_tbl_graph %>% # weight is how many viral genera they share
  activate(edges) %>% 
  ggraph(layout = "kk") + 
  geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) + 
  geom_node_point(size = 1, color = "firebrick") + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3) + 
  theme_void()

host_graph <- VIRION_DATA %>%
  filter(HostClass == "mammalia", VirusGenus %in% target_genera) %>% 
  select(Host, VirusGenus) %>% 
  table() %>% 
  graph_from_incidence_matrix() %>% 
  bipartite.projection() %>% 
  extract2("proj1")

host_gen_tab <- VIRION_DATA %>%
  filter(HostClass == "mammalia", VirusGenus %in% target_genera) %>% 
  select(Host, VirusGenus) %>% 
  table()

edges <- get.edgelist(host_graph)
virus_genera <- apply(edges, 1, function(e) {
  shared_gen <- intersect(
    names(which(host_gen_tab[e[1], ] > 0)),
    names(which(host_gen_tab[e[2], ] > 0))
  )
  paste(shared_gen, collapse = ", ")
})

# Store Vir Gen in the edge
E(host_graph)$ViralGen <- virus_genera


final_host_graph <- as_tibble(as_data_frame(host_graph)) %>%
  separate_rows(ViralGen, sep = ", ") %>% 
  graph_from_data_frame %>% 
  as_tbl_graph

edge_counts <- final_host_graph %>%
  activate(edges) %>%
  as_tibble() %>%
  group_by(ViralGen) %>%
  summarise(n_edges = n(), .groups = "drop")

# Join edge_counts to graph
final_host_graph <- final_host_graph %>%
  activate(edges) %>%
  left_join(edge_counts, by = "ViralGen") %>%
  activate(edges) %>%
  mutate(ViralGen = paste0(ViralGen, " (", n_edges, " edges)"))

set.seed(157)
ggraph(final_host_graph, layout = "kk") + #fr
  geom_edge_link(edge_alpha = 0.2, color = "steelblue4") +
  geom_node_point(size = 0.5, color = "grey", alpha = 0.5) +
  # geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  facet_wrap(~ViralGen, scales = "free", nrow = 4, ncol = 5) +
  theme(strip.text = element_text(
    size = 10))+
  theme_void()

# ggsave("m/ViralGeneraEdges.jpeg", width = 9, height = 6.5, dpi = 600)


# Looking at the filtered rasters ####

library(raster)

RasterList <- 
  "Data/MAMMALS_5km_final/Mammals_carnivora" %>% 
  dir_ls() %>% 
  map(raster)

names(RasterList) <- 
  "Data/MAMMALS_5km_final/Mammals_carnivora" %>% 
  list.files

RasterList[[1]] %>% plot
RasterList[[2]] %>% plot

RasterList$Panthera_leo.tif %>% plot

RasterList %>% 
  map(values)

RasterSums <- 
  RasterList %>% 
  map(values) %>% 
  map_dbl(sum)

# Trying to eliminate non-extent overlaps ####

library(sf)

SFList <- 
  "Data/DAMA Individual Species" %>% 
  dir_ls() %>% 
  map(st_read)

names(SFList) <- 
  "Data/DAMA Individual Species" %>% 
  list.files %>% 
  str_remove_all(".tif$")

library(ggregplot)

ExtentGet
GetExtent

ExtentList <- 
  SFList %>% 
  map(raster::extent) %>% 
  map(~data.frame(XMin = .x@xmin,
                  XMax = .x@xmax,
                  YMin = .x@ymin,
                  YMax = .x@ymax)) %>% 
  bind_rows(.id = "Species")






