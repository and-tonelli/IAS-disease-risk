require(tidyverse)
require(magrittr)
require(tidymodels)
require(plotly)
require(processx)

traits_comb <- combn(c("Overlapping_Cells_log", "foraging_sim", "phylo_sim", "trait_sim_gow"), 2)

names_tab <- c("Overlapping_Cells_log" = "Geographic overlap",
               "foraging_sim" = "Foraging similarity",
               "phylo_sim" = "Phylogenetic similarity",
               "trait_sim_gow" = "Trait similarity")

# Palettes
heat_palette <- paletteer::paletteer_c("grDevices::Heat", n = 100)
heat_colorscale <- lapply(seq_along(heat_palette), function(i) {
  list((i - 1)/(length(heat_palette) - 1), as.character(heat_palette[i]))
})

inferno_palette <- paletteer::paletteer_c("grDevices::Inferno", n = 100) 
inferno_colorscale <- lapply(seq_along(inferno_palette), function(i) {
  list((i - 1)/(length(inferno_palette) - 1), inferno_palette[i])
})

final_fit_models <- grep("final_fit_", names(.GlobalEnv), value=TRUE)


# Loop over predictor combinations to obtain mean 3D partial dependence plots (Fig. 2)
for (i in seq(1:6)){
  
  preds_for_data <- NULL
  
  var1 <- traits_comb[1, i]
  var2 <- traits_comb[2, i]
  
  # Create homogenous grid for prediction (var1 and var2 can change,other predictors are held at median)
  MyData <- expand.grid(VAR1 = seq(from = if (var1 == "Overlapping_Cells_log"){0.301} else {0}, to = if (var1 == "Overlapping_Cells_log"){6.444} else {1}, by = if (var1 == "Overlapping_Cells_log"){0.120451} else {0.02}), VAR2 = seq(from = if (var2 == "Overlapping_Cells_log"){0.301} else {0}, to = if (var2 == "Overlapping_Cells_log"){6.444} else {1}, by = if (var2 == "Overlapping_Cells_log"){0.120451} else {0.02}))
  MyData$phylo_sim <- rep(median(DatasetMainModel$phylo_sim), nrow(MyData)) 
  MyData$log_sum_cit <- rep(median(DatasetMainModel$log_sum_cit), nrow(MyData))
  MyData$trait_sim_gow <- rep(median(DatasetMainModel$trait_sim_gow), nrow(MyData))
  MyData$Overlapping_Cells_log <- rep(median(DatasetMainModel$Overlapping_Cells_log), nrow(MyData))
  MyData$phylo_sim <- rep(median(DatasetMainModel$phylo_sim), nrow(MyData))
  MyData$foraging_sim <- rep(median(DatasetMainModel$foraging_sim), nrow(MyData))
  MyData %<>% select(-c(var1, var2))
  colnames(MyData)[1] <- var1
  colnames(MyData)[2] <- var2
  
  # Loop over model formulations
  for (m in final_fit_models){
    
    final_fit_m <- get(m)
    
    preds1 <- predict(final_fit_m, MyData, type = "prob")
    preds_for_data <- cbind(preds_for_data, preds1$.pred_1)
    
  }
  
  # Get mean and upper and lower bound of 95% CI
    lower_bounds <- apply(preds_for_data, 1, quantile, probs = 0.025)
    upper_bounds <- apply(preds_for_data, 1, quantile, probs = 0.975)
  
    preds_for_data %<>% rowMeans()
    
    MyData <- cbind(probability = preds_for_data, MyData)
    MyData <- cbind(lower = lower_bounds, MyData)
    MyData <- cbind(upper = upper_bounds, MyData)
    
    
    probability_matrix <- matrix(preds_for_data, nrow = 51, ncol = 51)
    probability_low <- matrix(lower_bounds, nrow = 51, ncol = 51)
    probability_up <- matrix(upper_bounds, nrow = 51, ncol = 51)
    
    x_bins <- cut(DatasetMainModel %>% select(var1) %>% .[,1], breaks = 51)
    y_bins <- cut(DatasetMainModel %>% select(var2) %>% .[,1], breaks = 51)
    
    freq_tab <- table(x_bins, y_bins)
    z_matrix <- as.matrix(freq_tab)
    
    # First Surface: probability matrix
    fig1 <- plot_ly(
      x = unique(MyData %>% select(var1) %>% .[,1]),
      y = unique(MyData %>% select(var2) %>% .[,1]),
      z = t(probability_matrix),
      surfacecolor = t(probability_matrix),
      cmin = 0,
      cmax = 1,
      colorscale = inferno_colorscale,
      type = "surface",
      showscale = F,
      opacity = 1,
      colorbar = list(title = "Sharing probability", x = 1),
      name = "SP"
    )
    
    # Second Surface: data density at bottom
    fig2 <- plot_ly(
      x = unique(MyData %>% select(var1) %>% .[,1]),
      y = unique(MyData %>% select(var2) %>% .[,1]),
      z = matrix(-0.2, nrow = 51, ncol = 51),
      surfacecolor = t(log10(1 + z_matrix)),
      cmin =0,
      cmax = 4,
      colorscale = heat_colorscale,
      type = "surface",
      opacity = 1,
      showscale = F,
      colorbar = list(title = "Log(1 + count)", x = 1),
      name = "H"
    )
    
    fig_low <- plot_ly(
      x = unique(MyData %>% select(var1) %>% .[,1]),
      y = unique(MyData %>% select(var2) %>% .[,1]),
      z = t(probability_low),
      surfacecolor = t(probability_low),
      cmin = 0,
      cmax = 1,
      colorscale = inferno_colorscale,
      type = "surface",
      showscale = F,
      opacity = 0.7,
      colorbar = list(title = "Sharing probability", x = 1),
      name = "SP1"
    )
    
    fig_up <- plot_ly(
      x = unique(MyData %>% select(var1) %>% .[,1]),
      y = unique(MyData %>% select(var2) %>% .[,1]),
      z = t(probability_up),
      surfacecolor = t(probability_up),
      cmin = 0,
      cmax = 1,
      colorscale = inferno_colorscale,
      type = "surface",
      showscale = F,
      opacity = 0.7,
      colorbar = list(title = "Sharing probability", x = 1),
      name = "SP2"
    )
    
    final_fig <- subplot(fig1, fig2, shareX = TRUE, shareY = TRUE, titleX = TRUE, titleY = TRUE) %>%
      layout(
        scene = list(
          xaxis = list(title = names_tab[var1][[1]],
                       gridcolor='rgb(255, 255, 255)',
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          yaxis = list(title = names_tab[var2][[1]],
                       gridcolor='rgb(255, 255, 255)',
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          zaxis = list(title = "Sharing probability", range = c(-0.2, 1),
                       gridcolor='rgb(255, 255, 255)',
                       showticklabels = FALSE, 
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          aspectmode = "manual",
          aspectratio = list(x = 1, y = 1, z = 1),  
          camera = list(eye = list(x = -1.6, y = -1.6, z = 0.2))
        ),
        # paper_bgcolor = 'rgb(12,163,135)',
        plot_bgcolor = 'rgb(12,163,135)'
      )
    
    assign(paste0("3D_partial_", i), final_fig)
    
    htmlwidgets::saveWidget(as_widget(final_fig), paste0("m/", names(names_tab[var1]), "_x_", names(names_tab[var2]), ".html"))
    
    # Supplementary Figure with 95% CI
    uplow_fig <- subplot(fig_low, fig_up) %>%
      layout(
        scene = list(
          xaxis = list(title = names_tab[var1][[1]],
                       gridcolor='rgb(255, 255, 255)',
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          yaxis = list(title = names_tab[var2][[1]],
                       gridcolor='rgb(255, 255, 255)',
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          zaxis = list(title = "Sharing probability", range = c(0, 1),
                       gridcolor='rgb(255, 255, 255)',
                       showticklabels = FALSE, 
                       zerolinecolor='rgb(255, 255, 255)',
                       showbackground=TRUE,
                       backgroundcolor='rgb(230, 230,230)'
          ),
          aspectmode = "manual",
          aspectratio = list(x = 1, y = 1, z = 1),  
          camera = list(eye = list(x = -1.6, y = -1.6, z = 0.2))
        ),
        # paper_bgcolor = 'rgb(12,163,135)',
        plot_bgcolor = 'rgb(12,163,135)'
      )
    
    assign(paste0("3D_partial_CI_", i), uplow_fig)
    
    # Saving as interactive HTML plot
    htmlwidgets::saveWidget(as_widget(uplow_fig), paste0("m/", names(names_tab[var1]), "_x_", names(names_tab[var2]), "_CI.html"))
    
}

require(magick)

listfilespng <- list.files(pattern = "_x_.*\\.png$", recursive = T)
listfilesCI <- listfilespng[c(2, 4, 6, 8, 10, 12)]
listfilesmean <- listfilespng[-c(2, 4, 6, 8, 10, 12)]

trimmed_imgsCI <- lapply(listfilesCI, function(x) {
  image_trim(image_read(x))
})

# # 4. Combine into a 3x2 grid (use image_append by row)
# # First, group into two rows of 3 images
# scaled_imgs_row1 <- lapply(trimmed_imgsCI[1:3], function(img) image_scale(img, "x500"))
# row1 <- image_append(do.call(c, scaled_imgs_row1))
# 
# scaled_imgs_row2 <- lapply(trimmed_imgsCI[4:6], function(img) image_scale(img, "x500"))
# row2 <- image_append(do.call(c, scaled_imgs_row2))
# 
# final_grid <- image_append(c(row1, row2), stack = TRUE)
# image_write(final_grid, "m/combined_grid.png")


# labels letters and so on


# legend
# Palettes
# heat_palette <- paletteer::paletteer_d("PNWColors::Sunset2", n = 5)
# heat_colorscale <- lapply(seq_along(heat_palette), function(i) {
#   list((i - 1)/(length(heat_palette) - 1), as.character(heat_palette[i]))
# })

inferno_palette <- rev(paletteer::paletteer_d("rcartocolor::BluYl", n = 7) )
inferno_colorscale <- lapply(seq_along(inferno_palette), function(i) {
  list((i - 1)/(length(inferno_palette) - 1), inferno_palette[i])
})

# inferno_palette <- paletteer::paletteer_c("grDevices::Inferno", n = 100) 
# inferno_colorscale <- lapply(seq_along(inferno_palette), function(i) {
#   list((i - 1)/(length(inferno_palette) - 1), inferno_palette[i])
# })
# 
# inferno_palette <- rev(paletteer::paletteer_d("RColorBrewer::BuPu", n = 9)) 
# inferno_colorscale <- lapply(seq_along(inferno_palette), function(i) {
#   list((i - 1)/(length(inferno_palette) - 1), inferno_palette[i])
# })

heat_palette <- rev(paletteer::paletteer_d("MoMAColors::Ernst", n = 7))
heat_colorscale <- lapply(seq_along(heat_palette), function(i) {
  list((i - 1)/(length(heat_palette) - 1), as.character(heat_palette[i]))
})


legend_plot <- plot_ly(
  z = matrix(seq(0, 4, length.out = 10), ncol = 1),
  type = "heatmap",
  colorscale = heat_colorscale,
  showscale = TRUE,
  colorbar = list(title = "Log(1 + count)", len = 1,
                  tickfont = list(size=20))
) %>% layout(
  margin = list(l = 0, r = 0, t = 0, b = 0),
  xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
  yaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
  legend = list(font = list(size = 1, color = "black"))
)

legend_plot2 <- plot_ly(
  z = matrix(seq(0, 1, length.out = 10), ncol = 1),
  type = "heatmap",
  colorscale = inferno_colorscale,
  showscale = TRUE,
  colorbar = list(title = "Sharing probability", len = 1,
                  tickfont = list(size=20))
) %>% layout(
  margin = list(l = 0, r = 0, t = 0, b = 0),
  xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
  yaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
  legend = list(font = list(size = 1, color = "black"))
)

htmlwidgets::saveWidget(legend_plot, "legendRug.html")
htmlwidgets::saveWidget(legend_plot2, "legendProb.html")

legendProb <- image_read("legendProb.png")
legendRug <- image_read("legendRug.png")

# Crop a specific region:
# Syntax: "widthxheight+x_offset+y_offset"
# Example: Crop 800x600 area starting from (left=100, top=50)
croppedP <- image_crop(legendProb, geometry = "500x1000+1290+0")
# croppedP <- image_trim(croppedP)

croppedR <- image_crop(legendRug, geometry = "300x1000+1335+0")
# croppedR <- image_trim(croppedR)

labels <- letters[1:6]

trimmed_scaled_labeled <- lapply(seq_along(listfilesCI), function(i) {
  img <- image_read(listfilesCI[i])
  img <- image_trim(img)
  img <- image_annotate(
    img,
    text = labels[i],
    size = 30,
    gravity = "northwest",
    location = "+40+0",
    color = "black"
  )
  img <- image_scale(img, "x500")
  return(img)
})

sapply(trimmed_scaled_labeled, class)

# Combine into rows
row1 <- image_append(do.call(c, trimmed_scaled_labeled[1:3]))
row2 <- image_append(do.call(c, trimmed_scaled_labeled[4:6]))

# Combine into final grid
final_grid <- image_append(c(row1, row2), stack = TRUE)

# mean plots
trimmed_scaled_labeledMean <- lapply(seq_along(listfilesmean), function(i) {
  img <- image_read(listfilesmean[i])
  img <- image_trim(img)
  img <- image_annotate(
    img,
    text = labels[i],
    size = 30,
    gravity = "northwest",
    location = "+40+0",
    color = "black"
  )
  img <- image_scale(img, "x500")
  return(img)
})

sapply(trimmed_scaled_labeledMean, class)

# Combine into rows
row1m <- image_append(do.call(c, trimmed_scaled_labeledMean[1:3]))
row2m <- image_append(do.call(c, trimmed_scaled_labeledMean[4:6]))

final_gridmean <- image_append(c(row1m, row2m), stack = TRUE)

# Legend part
spacer <- image_blank(width = 20, height = 729, color = "white")
Legends<- image_append(c(croppedP, spacer, croppedR), stack = FALSE)

h_grid <- image_info(final_grid)$height
h_legends <- image_info(Legends)$height

pad_total <- h_grid - h_legends
pad_top <- floor(pad_total / 2)
pad_bottom <- ceiling(pad_total / 2)

# Create top and bottom spacers
top_pad <- image_blank(width = image_info(Legends)$width, height = pad_top, color = "white")
bottom_pad <- image_blank(width = image_info(Legends)$width, height = pad_bottom, color = "white")

# Stack legend with padding
Legends_centered <- image_append(c(top_pad, Legends, bottom_pad), stack = TRUE)
Legend_centered <- image_append(c(top_pad, croppedP, bottom_pad), stack = TRUE)

# final_gridL <- image_append(c(final_grid, spacer, Legends_centered), stack = F)
final_gridL <- image_append(c(final_gridmean, spacer, Legends_centered), stack = F)
final_gridL <- image_append(c(final_grid, spacer, croppedP), stack = F)

# Save
image_write(final_gridL, "m/combined_grid_mean.jpeg", quality = 100)

