
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
path_to_aoh <- "AOH_alien1975_5Km"

# Take all rasters
r_list <- list.files(path_to_aoh, pattern = ".tif", full.names = TRUE, recursive = TRUE)
require(tools)
alien_spp <- file_path_sans_ext(basename(r_list)) %>% str_replace("_", " ")

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
  as_tbl_graph %>%
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

host_graph_vir_lev <- VIRION_DATA %>%
  filter(HostClass == "mammalia",
         VirusFamily %in% target_fams) %>% 
  distinct(Host, Virus) %>% 
  table() %>% 
  graph_from_incidence_matrix() %>% 
  bipartite.projection() %>% 
  extract2("proj1")

host_gen_tab <- VIRION_DATA %>%
  filter(HostClass == "mammalia", VirusGenus %in% target_genera) %>% 
  distinct(Host, VirusGenus) %>% 
  table()

host_vir_tab <- VIRION_DATA %>%
  filter(HostClass == "mammalia", VirusFamily %in% target_fams) %>% 
  distinct(Host, Virus) %>% 
  table()

edges <- get.edgelist(host_graph)
virus_genera <- apply(edges, 1, function(e) {
  shared_gen <- intersect(
    names(which(host_gen_tab[e[1], ] > 0)),
    names(which(host_gen_tab[e[2], ] > 0))
  )
  paste(shared_gen, collapse = ", ")
})

edges <- get.edgelist(host_graph_vir_lev)
virus <- apply(edges, 1, function(e) {
  shared_vir <- intersect(
    names(which(host_vir_tab[e[1], ] > 0)),
    names(which(host_vir_tab[e[2], ] > 0))
  )
  paste(shared_vir, collapse = ", ")
})

VirGenTab <- VIRION_DATA %>% 
  distinct(Virus, VirusGenus)

VirFamTab <- VIRION_DATA %>% 
  distinct(Virus, VirusFamily)

# Store Vir Gen in the edge
E(host_graph_vir_lev)$Virus <- virus

final_host_graph <- as_tibble(as_data_frame(host_graph_vir_lev)) %>%
  separate_rows(Virus, sep = ", ") %>% 
  merge(., VirFamTab) %>%
  relocate(from, to) %>% 
  mutate(alien = ifelse(from %in% alien_spp | to %in% alien_spp, 1, 0)) %>% 
  graph_from_data_frame %>% 
  as_tbl_graph


edge_counts <- final_host_graph %>%
  activate(edges) %>%
  as_tibble() %>%
  group_by(VirusFamily) %>%
  summarise(n_edges = n(), .groups = "drop")

# Join edge_counts to graph
final_host_graph <- final_host_graph %>%
  activate(edges) %>%
  left_join(edge_counts, by = "VirusFamily") %>%
  activate(edges) %>%
  mutate(ViralFam = paste0(VirusFamily, " (", n_edges, " edges)"))

set.seed(123)
ggraph(final_host_graph, layout = "kk") + #fr
  geom_node_point(size = 0.5, color = "grey", alpha = 0.5) +
  geom_edge_link(edge_alpha = 0.2, aes(color = factor(alien))) +
  scale_edge_color_manual(values = c("1" = "firebrick", "0" = "steelblue4")) +
  # geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  facet_wrap(~ViralFam, scales = "free", nrow = 3, ncol = 6) +
  theme(strip.text = element_text(
    size = 10))+
  theme_void()+
  theme(legend.position = "none")

ggsave("m/ViralFamiliesEdges.jpeg", width = 10, height = 5, dpi = 600)


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

BlankMatrix <- matrix(NA, nrow = nrow(ExtentList), ncol = nrow(ExtentList))

i <- 1

for(i in 1:nrow(ExtentList)){
  
  ExtentList[i,] %>% print
  
  for(j in i:nrow(ExtentList)){
    
    ifelse(
      
      ExtentList[i,"YMin"] > ExtentList[j, "YMax"], 0,
      
      ifelse(
        
        ExtentList[i,"XMin"] > ExtentList[j, "XMax"], 0,
        
      )
      
    )
    
  }
  
}


# Alien endemic overlaps (13/94 aliens)
alien_endemic_network1975_2 <- read_csv("Data/alien_endemic_network1975_2.csv")

alien_endemic_network1975_2 %<>% 
  filter(Overlapping_Cells != 0)

alien_endemic_network1975_2 %>% 
  relocate(Alien_Species, Native_Species) %>% 
  graph_from_data_frame(directed = F) %>% 
  as_tbl_graph %>% 
  activate(nodes) %>%
  mutate(alien = ifelse(name %in% alien_endemic_network1975_2$Alien_Species, 1, 0)) %>% 
  ggraph(layout = "kk")+
  geom_edge_link(edge_alpha = 0.2, color = "grey80") +
  geom_node_point(size = 1, alpha = 0.5, aes(color = factor(alien))) +
  scale_color_manual(values = c("1" = "firebrick", "0" = "steelblue4")) +
  geom_node_text(aes(filter = alien == 1, label = name), repel = TRUE, size = 2)+
  theme_void()+
  theme(legend.position = "none")

ggsave("m/Alien_endemic_Net_2.jpeg", width = 12, height = 8, dpi = 600)

# Looking at sharing script etc #####

alien_endemic_network1975_2 <- read_csv("Data/alien_endemic_network.csv")

alien_endemic_network1975_2 %<>% 
  filter(Overlapping_Cells != 0)

alien_endemic_network1975_2 %>% 
  relocate(Alien_Species, Native_Species) %>% 
  graph_from_data_frame(directed = F) %>% 
  as_tbl_graph %>% 
  activate(nodes) %>%
  mutate(alien = ifelse(name %in% alien_endemic_network1975_2$Alien_Species, 1, 0)) %>% 
  ggraph(layout = "kk")+
  geom_edge_link(edge_alpha = 0.2, color = "grey80") +
  geom_node_point(size = 1, alpha = 0.5, aes(color = factor(alien))) +
  scale_color_manual(values = c("1" = "firebrick", "0" = "steelblue4")) +
  geom_node_text(aes(filter = alien == 1, label = name), repel = TRUE, size = 2)+
  theme_void()+
  theme(legend.position = "none")

# Making some models ####

SharingEdgeLists <- 
  "Data" %>% dir_ls(regex = "Sharing") %>% map(read.csv)

SharingEdgeLists %>% map(nrow)

library(ggregplot)

SharingEdgeLists %>% map(1) %>% map(nunique)

SharingEdgeLists[[1]] %>% write.csv("AAAAAA2.csv")

SharingEdgeLists[[1]] %>% 
  dplyr::select(1:2) %>% 
  table() %>% as.matrix %>% object.size

SpatialNetwork <- read.csv("Data/endemic_endemic_network.csv")

SpatialNetwork$Species2 %>% table

SpatialNetwork

DyadSharingEdgeList <- 
  SharingEdgeLists[[1]] %>% 
  rename_all(~str_replace(.x, "1$", "A") %>% str_replace("2$", "B")) %>% 
  bind_rows(SharingEdgeLists[[1]] %>% 
              rename_all(~str_replace(.x, "1$", "B") %>% str_replace("2$", "A"))) %>% 
  mutate_at(vars(matches("Species")), ~str_replace(.x, " ", "_"))

DyadSpatialNetwork <- 
  SpatialNetwork %>% 
  rename_all(~str_replace(.x, "1$", "A") %>% str_replace("2$", "B")) %>% 
  bind_rows(SpatialNetwork %>% 
              rename_all(~str_replace(.x, "1$", "B") %>% str_replace("2$", "A")))

FullEdgeList <- 
  DyadSpatialNetwork %>% 
  inner_join(DyadSharingEdgeList, by = c("SpeciesA", "SpeciesB"))

VarNames <- FullEdgeList %>% 
  dplyr::select(3:ncol(.)) %>% 
  names

AdjList <- 
  VarNames %>% 
  map(function(a){
    
    FullEdgeList %>% graph_from_data_frame(directed = F) %>% 
      get.adjacency(attr = a, sparse = F)
    
  })# %>% map(~.x/2)

AdjList[[1]]

FullDataFrame <- 
  AdjList %>% 
  map(reshape2::melt) %>% 
  reduce(~full_join(.x, .y, by = c("Var1", "Var2")))

names(FullDataFrame)[3:7] <- VarNames

FullDataFrame %>% dim

FullDataFrame %>% filter(n_virus == 0) %>% 
  select_if(is.numeric) %>% 
  colSums

TestDF <- 
  FullDataFrame %>% 
  filter(n_virus>0)

TestDF %>% 
  ggplot(aes(Overlapping_Cells, n_virus)) +
  geom_point()

library(cowplot)

theme_set(theme_cowplot())

TestDF %>% 
  ggplot(aes(Overlapping_Cells, n_virus)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10()

TestDF %>% 
  ggplot(aes(Overlapping_Cells, n_virus)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_log10() +
  scale_x_log10()

library(mgcv)

BAM1 <- bam(data = TestDF,
            n_virus ~ Overlapping_Cells,
            family = nb()
            
)

BAM1 %>% summary

BAM1 %>% plot

BAM2 <- bam(data = TestDF,
            # n_virus ~ s(Overlapping_Cells),
            n_virus ~ t2(Overlapping_Cells, Phylogeny),
            family = nb()
            
)

BAM2 %>% summary

BAM2 %>% plot

BAM2 <- bam(data = TestDF,
            # n_virus ~ s(Overlapping_Cells),
            n_virus ~ t2(Overlapping_Cells, Phylogeny),
            family = nb()
            
)

BAM2 %>% summary

BAM2 %>% plot

# Multi-membership random effect ####

TestDF %<>% arrange(Var1, Var2)

FullSpecies <- TestDF %>% dplyr::select(1:2) %>% unlist %>%  unique %>% sort

TestDF %<>% mutate_at(c("Var1", "Var2"), ~factor(.x, levels = FullSpecies))

MZ1 <- model.matrix( ~ Var1 - 1, data = TestDF) %>% as.matrix
MZ2 <- model.matrix( ~ Var2 - 1, data = TestDF) %>% as.matrix

SppMatrix = MZ1 + MZ2

TestDF$Spp <- SppMatrix

ParaPen <- 
  list(Spp = list(rank = nunique(FullSpecies), 
                  diag(nunique(FullSpecies))))

BAMMM <- bam(n_virus ~ s(Overlapping_Cells) + Spp,
             data = TestDF, 
             family = nb(),
             paraPen = ParaPen
)



