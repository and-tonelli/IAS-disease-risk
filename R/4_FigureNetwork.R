library(tidyverse)
library(ggraph)
library(tidygraph)
library(scales)
require(ggforce)
require(igraph)

FinalPreds <- readRDS("Data/FinalPreds.rds")

# student t with 9 df to estimate uncertainty around the median n of pred edges
library(distributions3)
T_9 <- StudentsT(df = 9) # 10 model replicates, 9 df

combined_df %>% #from 2_Modelling.R
  filter(vir_bin == 0) %>% 
  group_by(Species2, id) %>% 
  summarise(sum_events = sum(pred_bin)) %>% 
  group_by(Species2) %>% 
  summarise(median_events = round(median(sum_events), 0),
            sd_events = sd(sum_events),
            V = round(1.96*sd_events/sqrt(10), 0),
            V2 = round(quantile(T_9, 1 - 0.05 / 2) * sd(sum_events) / sqrt(10), 0),
            low_events = median_events - V2,
            up_events = median_events + V2) -> SummarisedN

#### Network Figure ####
# Top 20% alien sharers
quantile(probs = 0.80, SummarisedN$median_events)
top_alien_sp <- SummarisedN %>% filter(median_events >= 13) %>% pull(Species2)

# Alien edges
DatasetA <- FinalPredictions %>% 
  merge(summarized_data %>% rename(biogeographical_realm1 = biogeographical_realm), 
        by.x = "Species1", by.y = "iucn2020_binomial") %>%
  merge(summarized_data %>% rename(biogeographical_realm2 = biogeographical_realm), 
        by.x = "Species2", by.y = "iucn2020_binomial") %>%
  filter(pred_bin == 1) %>% 
  select(c(1, 2, 12, 13)) %>% 
  filter(Species2 %in% top_alien_sp) %>% 
  mutate(Species1 = ifelse(Species1 %in% top_alien_sp, str_replace(Species1, " ", "_"), Species1)) %>% 
  mutate(Alien = TRUE)

# Native edges
DatasetM <- DatasetMainModel %>%
  merge(summarized_data %>% rename(biogeographical_realm1 = biogeographical_realm), 
        by.x = "Species1", by.y = "iucn2020_binomial") %>%
  merge(summarized_data %>% rename(biogeographical_realm2 = biogeographical_realm), 
        by.x = "Species2", by.y = "iucn2020_binomial", suffixes = c("1", "2")) %>%
  filter(vir_bin == 1,
         Species1 %in% DatasetA$Species1,
         Species2 %in% DatasetA$Species1) %>% 
  mutate(Species1 = ifelse(Species1 %in% top_alien_sp, str_replace(Species1, " ", "_"), Species1),
         Species2 = ifelse(Species2 %in% top_alien_sp, str_replace(Species2, " ", "_"), Species2)) %>% # These are top ten aliens but also native that share with top ten aliens (changing label to treat the native counterpart as different node)
  select(c(1, 2, 16, 17)) %>% 
  mutate(Alien = FALSE)

# Shorten spp names for labels
itc <- sub("^(\\w)\\w*\\s+(\\w+)$", "\\1. \\2", top_alien_sp)

# Combine edges
edges_all <- rbind(DatasetM, DatasetA)

unique(edges_all$biogeographical_realm2)

# Build graph
graph_obj <- graph_from_data_frame(edges_all, directed = FALSE) %>%
  as_tbl_graph()

# Add node info: biogeographic realm + alien node flag
graph_obj <- graph_obj %>%
  activate(nodes) %>%
  mutate(
    biogeographical_realm = summarized_data$biogeographical_realm[match(name, summarized_data$iucn2020_binomial)],
    `Alien species` = name %in% DatasetA$Species2  # Only flag Species2 from DatasetA
  )

V(graph_obj)

# Reattach edge alien flags
graph_obj <- graph_obj %>%
  activate(edges) %>%
  mutate(
    `Alien species` = edge_attr(graph_obj, "Alien species") 
  )

realms <- graph_obj %>% 
  activate(nodes) %>%
  pull(biogeographical_realm) %>%
  unique()

# Make palette
realm_colors <- c("#406FB6", "#65BFBB", "#ADDA54", "grey60", "#FFD787", "#F39E6C") #0A3DCB
names(realm_colors) <- rev(str_sort(realms))

final_colors <- c("Alien species" = "#C41919", realm_colors) # Add red for alien

# Plot
install.packages("graphlayouts", dep = T)
require(graphlayouts)

# Grouping nodes by realm cluster
bb <- layout_as_backbone(graph_obj, keep = 0.4)

E(graph_obj)$col <- FALSE
E(graph_obj)$col[bb$backbone] <- TRUE

# Network Figure
ggraph(graph_obj, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
  # Alien edges red, others grey
  geom_edge_link(aes(color = Alien),  alpha = 0.1, show.legend = FALSE) +
  scale_edge_color_manual(values = c("TRUE" = "#C41919", "FALSE" = "grey75")) +
  
  # Nodes colored by realm or red if alien
  geom_node_point(aes(color = ifelse(`Alien species`, "Alien species", biogeographical_realm)),
                  size = 2, alpha = 0.85) +
  scale_color_manual(values = final_colors) +
  
  # Labels only for top alein nodes
  geom_node_text(aes(label = ifelse(`Alien species`, itc[match(name, top_alien_sp)], "")),
                 repel = TRUE, size = 4.5, fontface = "italic", color = "white", max.overlaps = 50,
                 segment.color = "white",           
                 segment.size = 0.4,
                 min.segment.length = 0.2,
                 max.time = 1) +
  
  theme_void() +
  theme(
    legend.position = c(0.8, 0.22),
    plot.background = element_rect(fill = "#2e2e2e", color = NA),
    panel.background = element_rect(fill = "#2e2e2e", color = NA),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white", size = 16)
  )+
  labs(color = "")

#### Sharing Map ####
Alien_distances <- readRDS("Data/Alien_distances2.rds")

FinalPredictions
  
SharingEncounters <- left_join(Alien_distances %>% group_by(Binomial) %>% 
                                 summarise(mean_dist = mean(Distance)), Fin_predictions, by = c("Binomial"= "Species2"))

CompleteCentroids <- read_csv("Data/CompleteCentroids.csv")
CompleteCentroids$Species1 <- str_replace(CompleteCentroids$Species1, "_", " ")

CompletePoints <- merge.data.frame(CompleteCentroids, Alien_distances[c(1, 10, 11)] %>% distinct(), by = 1) %>% 
  rename("LonA" = Centroid_X,
         "LatA" = Centroid_Y)

# Selecting only novel links
CompletePoints %>% 
  merge.data.frame(., SharingEncounters[c(1, 3:7)], by.x = c(1, 2), by.y = c(1, 2)) %>% 
  filter(pred_bin == 1, vir_bin == 0) %>% 
  group_by(Alien, LonA, LatA, LonN, LatN, Landmass) %>% 
  summarise(n = n()) -> link_data

bbox <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))

ggplot() +
  geom_sf(data = bbox, fill = "#2e2e2e", color = "#2e2e2e")+
  geom_sf(data = world, fill = "gray45", color = "gray45", linewidth = 0.2) +
  
  # Arcs
  geom_curve(
    data = link_data,
    aes(x = LonA, y = LatA, xend = LonN, yend = LatN),
    curvature = 0.4,
    color = "sandybrown",
    alpha = 0.30,
    linewidth = 0.2
  ) +
  
  # Origin nodes (x)
  geom_point(
    data = link_data,
    aes(x = LonN, y = LatN),
    color = "sandybrown",
    fill = "sandybrown",
    size = 2,
    alpha = 0.9,
    shape = 4, #21
    stroke = 0.6
  ) +
  
  # Destination nodes
  geom_point(
    data = link_data,
    aes(x = LonA, y = LatA, fill = n),
    shape = 21,
    size = 3,
    color = "gray10",
    stroke = 0.7,
    alpha = 0.7
  ) +
  
  # Red palette
  paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold", name = "Sharing\nevents\n", breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) +
  
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = NULL) +
  # facet_wrap(~ order)+
  theme(
    panel.background = element_rect(fill = "#2e2e2e", color = NA), ##1e1e1e
    plot.background = element_rect(fill = "#2e2e2e", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.text = element_text(color = "white", size = 18),
    legend.title = element_text(color = "white", size = 20),
    plot.title = element_text(color = "white"),
    plot.subtitle = element_text(color = "white")
  )+
  theme(legend.position = "right",
        legend.box = "right",
        legend.background = element_rect(fill = "#2e2e2e", color = NA),
        legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(2.5, 'cm'),
        legend.text = element_text(size = 20),
        legend.ticks = element_blank(),
        strip.text = element_text(size = 14),
        plot.margin = margin(1, 1, 1, 1, "cm")) + # top right bottom left
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
         size = guide_legend(title.position = "top", title.hjust = 0.5)) 

