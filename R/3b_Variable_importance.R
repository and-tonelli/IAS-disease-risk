# Variable importance
# Improvement in accuracy contributed by a feature to the branches it is used in

varimp_df_all <- NULL

# Loop over model formulations 
for (m in final_fit_models){
  
  final_fit_m <- get(m)
  
  print(paste0("Currently at model ", which(final_fit_models == m), "/10"))
  
  fitted_mod <- extract_fit_parsnip(final_fit_m)
  fitted_xgb <- fitted_mod$fit
  
  varimp_df <- xgb.importance(model = fitted_xgb) %>% 
    mutate(model = which(final_fit_models == m))
  
  varimp_df_all <- rbind(varimp_df_all, varimp_df)
  
}

# Summarise variable importance
importance_summary <- varimp_df_all %>%
  as_tibble() %>%
  group_by(Feature) %>%
  summarise(
    median_gain = median(Gain),
    ymin = quantile(Gain, 0.25),
    ymax = quantile(Gain, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(median_gain)) %>%
  mutate(Feature = factor(Feature, levels = Feature))  # maintain order

# Plot importance
ggplot(importance_summary, aes(x = Feature, y = median_gain)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), linewidth = 2, color = "gray25") +
  geom_point(size = 4, color = "steelblue") +
  theme_minimal() +
  labs(
    x = "",
    y = "Variable importance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())+
  scale_x_discrete(labels = names_tab)+
  ylim(c(0, 0.7))
