## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# get vector of sample_ids to loop through
sample_id <- TMM_normalised$sample_id

# get tibble of other metadata
category_scores <- TMM_normalised %>% 
  select(sample_id, plate_number, treatment, cluster)



#load dataframe that created list with shared non chaninging ABs

cell_state_model <- read_csv2(file = model)
cell_state_model$X12 <- NULL
cell_state_model$X13 <- NULL
cell_state_model$X14 <- NULL

# get categories to run through
categories <- factor(colnames(cell_state_model[2:length(cell_state_model)]), levels = colnames(cell_state_model[2:length(cell_state_model)]))

# transform TMM_normalised into long format
TMM_normalised <- TMM_normalised %>% 
  pivot_longer(names_to = "antibody", values_to = "intensity", 5:length(TMM_normalised))


#transform cell state model into long format
cell_state_model <- cell_state_model %>% 
  pivot_longer(names_to = "type", values_to = "effect", 2:length(cell_state_model))

cell_state_model$effect <- cell_state_model$effect %>% 
  recode("pos" = 1L, "neg" = -1L)

# calculations

for (category in categories) {
  # get antibodies that affect category score
  print(paste("calculating score for", category))
  
  b <- cell_state_model %>% 
    filter(type == category) %>% 
    filter(effect == -1L | effect == 1L)
  
  antibodies <- b$antibody
  
  effect_size <- b$effect
  
  # filter on the antibodies that are playing a role in this category
  TMM_temp <- TMM_normalised %>% 
    filter(antibody %in% antibodies)
  
  # initialise vector score to store scores in
  score <- vector()
  
  # select a cell to calculate on
  for (cell in sample_id) {
    c <- TMM_temp %>% 
      filter(sample_id == cell)
    
    # reverse the score for negative effects
    c$adjusted_score <- c$intensity * effect_size
    
    # calculate the score for this cell
    sum_score <- sum(c$adjusted_score)
    
    #add score to the vector score
    score <- append(score, sum_score)
  }
  
  # rescale score to 0:10 scale
  score <- scales::rescale(score, to = c(0, 10))
  
  # add score to the tibble category scores
  category_scores[category] <- score
  
}


if (save_table == "yes") {
  write_tsv(category_scores , path = paste(data_dir, "/cell_state_scores_all_treatments.tsv", sep = ""))
}

# transform category_scores into long format for plotting
category_scores <- category_scores %>% 
  pivot_longer(names_to = "category", values_to = "score", 5:length(category_scores))

category_scores_plot <- category_scores %>% 
  group_by(treatment, category) %>% 
  summarise(score = mean(score))

# set categories with the right levels
category_scores_plot$category <- factor(category_scores_plot$category, levels = categories)

# transform clusters into factors
category_scores_plot$treatment <- factor(category_scores_plot$treatment)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_cell_state_scores_cell_state_scores.pdf", sep = ""),
      height = 4,
      width = 6)
}


plot_categories <- ggplot(data = category_scores_plot, aes(y=category, x=treatment, col=score, size=score)) +
  geom_point() +
  scale_colour_gradientn(colours = c("darkblue", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "darkred"),
                         limits = c(0,10))+
  scale_size_continuous() +
  xlab("Treatment") +
  ylab("") +
  labs(size = "Cell state score", col = "Cell state score") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"))

plot_categories

if (save_figure == "yes") {
  dev.off()
  print("figures have been saved")
}

## FOR VEH VS ip70S6K

category_scores_filtered <- category_scores %>% 
  filter(treatment == "EGF" | treatment == "ip70S6K_EGF")

# set categories with the right levels
category_scores_filtered$category <- factor(category_scores_filtered$category, levels = categories)

# transform clusters into factors
category_scores_filtered$treatment <- factor(category_scores_filtered$treatment)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/violin_plots_cell_state_scores_veh_vs_ip70S6K.pdf", sep = ""),
      height = 3,
      width = 3)
}

plot_list <- vector('list', length(categories))
counter <- 1

for(variable in categories) {
  TMM_temp <- category_scores_filtered %>% 
    filter(category == variable)
  
  #create separate data frames for each treatment
  table_treatment_1 <- TMM_temp %>% 
    filter(treatment == "EGF") %>% 
    pull(score)
  
  table_treatment_2 <- TMM_temp %>% 
    filter(treatment == "ip70S6K_EGF") %>% 
    pull(score)
  
  p_val <- ks.test(table_treatment_1, table_treatment_2)$p.value
  
  p_val <- p.adjust(p_val, method = "BH")
  
  if(p_val <= 0.05) {
    sig <- "*"
  }
  if(p_val <= 0.005) {
    sig <- "**"
  }
  if(p_val <= 0.0005) {
    sig <- "***"
  }
  if(p_val <= 0.00005) {
    sig <- "****"
  }
  
  plot <- ggviolin(data = TMM_temp, x = "treatment", y = "score", col = "treatment",
           add = "jitter", add.params = list(alpha = 0.5)) +
    geom_text(aes(x = 1.5, y = 9),
              label = sig,
              size = 8) +
    scale_colour_manual(values = palette_plots) +
    facet_wrap(facets = "category", ncol = 5) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom") +
    ylab("Cell state score") +
    xlab("Treatment") +
    coord_cartesian(ylim = c(0,10))
  
  print(plot)
  
  plot_list[[counter]] <- plot
  
  counter <- counter + 1
}


if(save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/violin_plots_cell_state_scores_veh_vs_ip70S6K_combined.pdf", sep = ""),
      height = 6,
      width = 15)
}

plot_grid(plotlist = plot_list, ncol = 5)

if(save_figure == "yes") {
  dev.off()
}







## FOR VEH VS iRSK

category_scores_filtered <- category_scores %>% 
  filter(treatment == "EGF" | treatment == "iRSK_EGF")

# set categories with the right levels
category_scores_filtered$category <- factor(category_scores_filtered$category, levels = categories)

# transform clusters into factors
category_scores_filtered$treatment <- factor(category_scores_filtered$treatment)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/violin_plots_cell_state_scores_veh_vs_iRSK.pdf", sep = ""),
      height = 3,
      width = 3)
}

plot_list <- vector('list', length(categories))
counter <- 1

for(variable in categories) {
  TMM_temp <- category_scores_filtered %>% 
    filter(category == variable)
  
  #create separate data frames for each treatment
  table_treatment_1 <- TMM_temp %>% 
    filter(treatment == "EGF") %>% 
    pull(score)
  
  table_treatment_2 <- TMM_temp %>% 
    filter(treatment == "iRSK_EGF") %>% 
    pull(score)
  
  p_val <- ks.test(table_treatment_1, table_treatment_2)$p.value
  
  p_val <- p.adjust(p_val, method = "BH")
  
  if(p_val <= 0.05) {
    sig <- "*"
  }
  if(p_val <= 0.005) {
    sig <- "**"
  }
  if(p_val <= 0.0005) {
    sig <- "***"
  }
  if(p_val <= 0.00005) {
    sig <- "****"
  }
  
  plot <- ggviolin(data = TMM_temp, x = "treatment", y = "score", col = "treatment",
                   add = "jitter", add.params = list(alpha = 0.5)) +
    geom_text(aes(x = 1.5, y = 9),
              label = sig,
              size = 8) +
    scale_colour_manual(values = c(palette_plots[1], palette_plots[3])) +
    facet_wrap(facets = "category", ncol = 5) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom") +
    ylab("Cell state score") +
    xlab("Treatment") +
    coord_cartesian(ylim = c(0,10))
  
  print(plot)
  
  plot_list[[counter]] <- plot
  
  counter <- counter + 1
}


if(save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/violin_plots_cell_state_scores_veh_vs_iRSK_combined.pdf", sep = ""),
      height = 6,
      width = 15)
}

plot_grid(plotlist = plot_list, ncol = 5)

if(save_figure == "yes") {
  dev.off()
}
