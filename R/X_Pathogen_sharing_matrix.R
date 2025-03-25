library(igraph);
library(ggraph);
library(tidygraph)
require(tidyverse); require(magrittr)

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
  stringsAsFactors = FALSE
)

colnames(df_virion) <- c("Species1", "Species2", "n_virus")
write.csv(df_virion, "Data/Sharing_all.csv", row.names = F) # Species that share at least one virus with others


# Sharing of any virus of target ("alien") viral families
vir_tab <- read_csv("Data/vir_tab_long_final.csv")
target_family <- Virion_complete %>% filter(Virus %in% vir_tab$virus) %>% pull(VirusFamily) %>% unique

for (vir in target_family){
  
  Virion_complete %>%
    filter(VirusFamily == vir) %>% 
    select(Host, Virus) %>% 
    distinct() %>% 
    table() %>% 
    graph_from_incidence_matrix() -> bipartite_graph
  
  projections <- bipartite.projection(bipartite_graph)
  host_proj <- projections$proj1
  weights <- E(host_proj)$weight
  edge_list <- get.edgelist(host_proj)
  
  df_virion <- data.frame(
    Species1 = edge_list[, 1],
    Species2 = edge_list[, 2],
    Weight = weights,
    stringsAsFactors = FALSE
  )
  colnames(df_virion) <- c("Species1", "Species2", "n_virus")
  
  write.csv(df_virion, paste0("Data/Sharing_", vir, ".csv"), row.names = F)
  
  
}


# PLOT STUFF
# connected_by_vir <- lapply(V(bipartite_graph)[type == F], function(v) {
#   neighbors(bipartite_graph, v)$name
# })
# 
# edge_list <- as_data_frame(host_proj, what = "edges") %>%
#   rowwise() %>%
#   mutate(
#     shared_viruses = list(intersect(connected_by_vir[[from]], connected_by_vir[[to]]))
#   ) %>%
#   ungroup()
# 
# g <- as_tbl_graph(host_proj) %>%
#   activate(edges) %>%
#   mutate(virus_name = edge_list$shared_viruses) %>% 
#   activate(nodes) %>% 
#   mutate(bat = ifelse(name %in% bats, 1, 0))
# 
# 
# g %>% 
#   ggraph(layout = "kk") + 
#   geom_edge_arc(edge_alpha = 0.4, color = "steelblue4") +
#   geom_node_point(size = 1, aes(color = factor(bat))) + 
#   # geom_node_text(aes(label = name), repel = TRUE, size = 3) + 
#   # scale_edge_color_manual(values = c("Betacoronavirus" = "#364C54FF",
#   #                                    "Alphacoronavirus" = "#94475EFF",
#   #                                    "Gammacoronavirus" = "#E5A11FFF",
#   #                                    "Unknown" = "grey80",
#   #                                    "Both" = "black")) +
#   scale_color_manual(values = c("1" = "firebrick", "0" = "black"))+
#   theme_void()
# 
# activate(nodes) %>%
#   mutate(alien = ifelse(name %in% alien_endemic_network$Alien_Species, 1, 0),
#          to_label = ifelse(name %in% mostconnected_spp, 1, 0)) %>% 
#   ggraph()+ #layout = "kk"
#   geom_edge_arc(edge_alpha = 0.4, color = "grey80") +
#   geom_node_point(size = 1, alpha = 0.5, aes(color = factor(alien))) +
#   scale_color_manual(values = c("1" = "firebrick", "0" = "#04225CFF")) +
#   geom_node_text(aes(filter = to_label == 1, label = name), repel = TRUE, size = 2)+
#   theme_void()+
#   theme(legend.position = "none")
# 
# host_graph <- Virion_complete %>%
#   filter(VirusGenus %in% target_genera) %>% 
#   select(Host, VirusGenus) %>% 
#   table() %>% 
#   graph_from_incidence_matrix() %>% 
#   bipartite.projection() %>% 
#   extract2("proj1")
# 
# host_graph_vir_lev <- VIRION_DATA %>%
#   filter(HostClass == "mammalia",
#          VirusFamily %in% target_fams) %>% 
#   distinct(Host, Virus) %>% 
#   table() %>% 
#   graph_from_incidence_matrix() %>% 
#   bipartite.projection() %>% 
#   extract2("proj1")
# 
# host_gen_tab <- VIRION_DATA %>%
#   filter(HostClass == "mammalia", VirusGenus %in% target_genera) %>% 
#   distinct(Host, VirusGenus) %>% 
#   table()
# 
# host_vir_tab <- VIRION_DATA %>%
#   filter(HostClass == "mammalia", VirusFamily %in% target_fams) %>% 
#   distinct(Host, Virus) %>% 
#   table()
# 
# edges <- get.edgelist(host_graph)
# virus_genera <- apply(edges, 1, function(e) {
#   shared_gen <- intersect(
#     names(which(host_gen_tab[e[1], ] > 0)),
#     names(which(host_gen_tab[e[2], ] > 0))
#   )
#   paste(shared_gen, collapse = ", ")
# })
# 
# edges <- get.edgelist(host_graph_vir_lev)
# virus <- apply(edges, 1, function(e) {
#   shared_vir <- intersect(
#     names(which(host_vir_tab[e[1], ] > 0)),
#     names(which(host_vir_tab[e[2], ] > 0))
#   )
#   paste(shared_vir, collapse = ", ")
# })
# 
# VirGenTab <- VIRION_DATA %>% 
#   distinct(Virus, VirusGenus)
# 
# VirFamTab <- VIRION_DATA %>% 
#   distinct(Virus, VirusFamily)
# 
# # Store Vir Gen in the edge
# E(host_graph_vir_lev)$Virus <- virus
# 
# final_host_graph <- as_tibble(as_data_frame(host_graph_vir_lev)) %>%
#   separate_rows(Virus, sep = ", ") %>% 
#   merge(., VirFamTab) %>%
#   relocate(from, to) %>% 
#   mutate(alien = ifelse(from %in% alien_spp | to %in% alien_spp, 1, 0)) %>% 
#   graph_from_data_frame %>% 
#   as_tbl_graph
# 
# 
# edge_counts <- final_host_graph %>%
#   activate(edges) %>%
#   as_tibble() %>%
#   group_by(VirusFamily) %>%
#   summarise(n_edges = n(), .groups = "drop")
# 
# # Join edge_counts to graph
# final_host_graph <- final_host_graph %>%
#   activate(edges) %>%
#   left_join(edge_counts, by = "VirusFamily") %>%
#   activate(edges) %>%
#   mutate(ViralFam = paste0(VirusFamily, " (", n_edges, " edges)"))
# 
# set.seed(123)
# ggraph(final_host_graph, layout = "kk") + #fr
#   geom_node_point(size = 0.5, color = "grey", alpha = 0.5) +
#   geom_edge_link(edge_alpha = 0.2, aes(color = factor(alien))) +
#   scale_edge_color_manual(values = c("1" = "firebrick", "0" = "steelblue4")) +
#   # geom_node_text(aes(label = name), repel = TRUE, size = 3) +
#   facet_wrap(~ViralFam, scales = "free", nrow = 3, ncol = 6) +
#   theme(strip.text = element_text(
#     size = 10))+
#   theme_void()+
#   theme(legend.position = "none")
# 
# ggsave("m/ViralFamiliesEdges.jpeg", width = 10, height = 5, dpi = 600)