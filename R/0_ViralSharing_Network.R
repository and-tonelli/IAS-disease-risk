# Making the viral sharing network

library(tidygraph)
require(tidyverse)
require(magrittr)

Virion_complete <- read_csv("Data/VirionMammals_ICTV.csv")

# Sharing of any virus
Virion_complete %>%
  select(Host, Virus) %>% 
  distinct() %>% 
  table() %>% 
  graph_from_biadjacency_matrix() -> bipartite_graph

projections <- bipartite_projection(bipartite_graph)
host_proj <- projections$proj1
weights <- E(host_proj)$weight
edge_list <- as_edgelist(host_proj)

df_virion <- data.frame(
  Species1 = edge_list[, 1],
  Species2 = edge_list[, 2],
  Weight = weights,
  stringsAsFactors = FALSE)

colnames(df_virion) <- c("Species1", "Species2", "n_virus")
write.csv(df_virion, "Data/Sharing_all.csv", row.names = F) # Species that share at least one virus with others