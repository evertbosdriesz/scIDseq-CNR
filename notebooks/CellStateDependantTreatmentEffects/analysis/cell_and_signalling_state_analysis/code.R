## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# filter out only EGF
TMM_normalised <- TMM_normalised %>% 
  filter(treatment == "EGF")

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
  write_tsv(category_scores , path = paste(data_dir, "/cell_state_scores.tsv", sep = ""))
}

# transform category_scores into long format for plotting
category_scores <- category_scores %>% 
  pivot_longer(names_to = "category", values_to = "score", 5:length(category_scores))

category_scores_plot <- category_scores %>% 
  group_by(cluster, category) %>% 
  summarise(score = mean(score))

# set categories with the right levels
category_scores_plot$category <- factor(category_scores_plot$category, levels = categories)

# transform clusters into factors
category_scores_plot$cluster <- factor(category_scores_plot$cluster)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_cell_state_scores_cell_state_scores.pdf", sep = ""),
      height = 4,
      width = 6)
}


plot_categories <- ggplot(data = category_scores_plot, aes(y=category, x=cluster, col=score, size=score)) +
  geom_point() +
  scale_colour_gradientn(colours = c("darkblue", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "darkred"),
                         limits = c(0,10))+
  scale_size_continuous() +
  xlab("Cluster") +
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


# open cycling markers
cycling_markers <- read_csv2(file = cycling_markers)

cycling_markers <- as_factor(cycling_markers$antibody)

# filter on cycling markers
TMM_temp <- TMM_normalised %>% 
  filter(antibody %in% cycling_markers)

TMM_temp <- TMM_temp %>%
  pivot_wider(names_from = 'antibody', values_from = 'intensity')

for (ab in cycling_markers) {
  TMM_temp[ab] <- scales::rescale(as_vector(TMM_temp[ab]), to = c(0, 10))
}

TMM_temp <- TMM_temp %>%
  pivot_longer(names_to = "antibody", values_to = "intensity", 5:11)


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/violin_plots_cycling_markers_cell_state_scores.pdf", sep = ""),
      height = 5,
      width = 5.5)
}

for (ab in cycling_markers) {
  TMM_temp1 <- TMM_temp %>% 
    filter(antibody == ab)
  
  # plot violin plots of intensities per cluster
  violin_plot <- ggviolin(data = TMM_temp1, x = "cluster", y = "intensity", col = "cluster", add = "jitter") +
    facet_wrap(facets = "antibody") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white", colour = "black", line = "solid"),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = NA)) +
    theme(strip.background = element_blank(),
          axis.text = element_text(colour = "black")) +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = palette_plots) +
    ylab("Score") +
    xlab("Cluster")
  
  plot(violin_plot)

}

if (save_figure == "yes") {
  dev.off()
}

TMM_temp <- TMM_temp %>% 
  group_by(cluster, antibody) %>% 
  summarise(intensity = mean(intensity))

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_cycling_markers_cell_state_scores.pdf", sep = ""),
      height = 2.5,
      width = 6)
}

plot_cycling_markers <- ggplot(data = TMM_temp, aes(y=antibody, x=factor(cluster), col=intensity, size=intensity)) +
  geom_point(show.legend = TRUE) +
  scale_colour_gradientn(colours = c("darkblue", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "darkred"),
                         limits = c(0,10))+
  scale_size_continuous() +
  xlab("Cluster") +
  ylab("") +
  labs(title = "Score for cell state categories", 
       col = "Cycling marker score",
       size = "Cycling marker score") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(plot.title = element_blank()) +
  theme(legend.text = element_text(size = 9),
        axis.text = element_text(colour = "black"))

plot_cycling_markers

if (save_figure == "yes") {
  dev.off()
}

# combine two plots for easier layout in illustratotr

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_cell_state_scores.pdf", sep = ""),
      height = 6.1,
      width = 6)
}

plot_grid(plot_cycling_markers, plot_categories, ncol = 1, align = "v", rel_heights = c(7,10))

if (save_figure == "yes") {
  dev.off()
}

# open data frame with umap coordinates
umap_information <- read_tsv(file = umap_coord)

# select meta data
umap_information <- umap_information %>% 
  select(1:6)

# filter on EGF to find coordinates for the scored cells
umap_EGF <- umap_information %>% 
  filter(treatment == "EGF")

category_scores <- category_scores %>% 
  pivot_wider(names_from = "category", values_from = "score")

umap_category_plot <- merge(umap_EGF, category_scores, by = c("sample_id", "plate_number", "treatment", "cluster"))

umap_category_plot <- umap_category_plot %>% 
  pivot_longer(names_to = "category", values_to = "score", 7:16)

umap_category_plot$category <- as_factor(umap_category_plot$category)
levels(umap_category_plot$category) <- categories

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_scores_seperate_pages.pdf", sep = ""),
      height = 3,
      width = 3.5)
}

for (type in categories) {
  TMM_temp <- umap_category_plot %>% 
    filter(category == type)
  
  umap <- ggplot(TMM_temp) + 
    geom_point(data = umap_information, aes(x=umap_1, y=umap_2), colour="grey",pch=1, size=1) +
    geom_point(aes(x=umap_1, y=umap_2, color=score), show.legend = TRUE) + 
    geom_point(aes(x=umap_1, y=umap_2), color = "black", pch=1, show.legend = FALSE) + 
    facet_wrap(facets = "category", ncol = 5) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom") +
    scale_colour_gradientn(colours = c("darkblue", "#0000CC", "white", "#FF0000", "darkred"), 
                           space = "Lab", values = c(0,
                                                     0.495,
                                                     0.505,
                                                     1))
  
  plot(umap)
}

if (save_figure == "yes") {
  dev.off()
}


# load comma separated file into dataframe

TMM_copy<- read_csv("R/01_data_in_use/TMM_normalised_with_clusters_matrix.csv")

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# get tibble of other metadata
all_scores <- TMM_copy %>% 
  select(sample_id, plate_number, treatment, cluster)

# get sample ids to loop through
sample_id <- TMM_copy$sample_id

for (ab in colnames(TMM_copy[5:length(TMM_copy)])) {
  TMM_copy[ab] <- scales::rescale(as_vector(TMM_copy[ab]), to = c(0, 10))
}

# transform TMM_normalised into long format
TMM_copy <- TMM_copy %>% 
  pivot_longer(names_to = "antibody", values_to = "intensity", 5:length(TMM_copy))

# create tibble to put normalised category scores in
category_scores_average_normalised <- tibble()

# calculate signalling markers scores for all treatments
for (category in categories) {
  # get antibodies that affect category score
  print(paste("calculating score for", category))
  
  b <- cell_state_model %>% 
    filter(type == category) %>% 
    filter(effect == -1L | effect == 1L)
  
  antibodies <- b$antibody
  
  effect_size <- b$effect
  
  # filter on the antibodies that are playing a role in this category
  TMM_temp <- TMM_copy %>% 
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
  
  # place score in dataframe for later use
  all_scores[category] <- score
}

for (cluster_no in c(1:length(unique(all_scores$cluster)))) {
  cluster_df <- all_scores %>% 
    filter(cluster == cluster_no)
  cluster_df_long <- cluster_df %>% 
    gather(key = 'category', value = 'score', -sample_id, -plate_number, -treatment, -cluster)
  for (type in sort(colnames(all_scores[5:length(all_scores)]))) {
    a <- cluster_df_long %>% 
      filter(treatment == "EGF") %>% 
      filter(category == type) %>% 
      summarise(mean = mean(score))
    a <- as.numeric(a)
    cluster_df[type] <- cluster_df[type] / a
  }
  category_scores_average_normalised <- rbind(category_scores_average_normalised, cluster_df)
}

category_scores_average_normalised <- category_scores_average_normalised %>% 
  pivot_longer(names_to = "category", values_to = "score", 5:length(category_scores_average_normalised))

category_scores_average_normalised$category <- as_factor(category_scores_average_normalised$category)

levels(category_scores_average_normalised$category) <- categories

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/boxplots_fold_changes_upon_treatment.pdf", sep = ""),
      height = 4,
      width = 5)
}

for (type in categories) {
  TMM_temp <- category_scores_average_normalised %>% 
    filter(category == type)
  
  boxplot <- ggboxplot(data = TMM_temp, x = "cluster", 
            y = "score", col = "treatment", add = "jitter",
            add.params = list(size = 1.6, alpha = 0.5)) +
    facet_wrap(facets = "category", ncol=2, scales = "free") +
    scale_colour_manual(values = palette_plots) +
    coord_cartesian(ylim = c(0,2)) +
    ylab("Fold change") +
    xlab("Cluster") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom")
  
  plot(boxplot)
}




if (save_figure == "yes") {
  dev.off()
}


# calculate average fold change per cluster and significance of that change

matrix_fold_changes <- tibble()

for (type in sort(colnames(all_scores[5:length(all_scores)]))) {
  df_temp <- category_scores_average_normalised %>% 
    filter(category == type)
  
  p_value <- vector()
  
  for (i in c(1:length(unique(category_scores_average_normalised$cluster)))) {
    
    df_temp_cluster <- df_temp %>% 
      filter(cluster == i)
    
    EGF <- df_temp_cluster %>% 
      filter(treatment == "EGF")
    
    EGF <- EGF$score
    
    iRSK <- df_temp_cluster %>% 
      filter(treatment == "iRSK_EGF")
    
    iRSK <- iRSK$score
    
    ip70S6K <- df_temp_cluster %>% 
      filter(treatment == "ip70S6K_EGF")
    
    ip70S6K <- ip70S6K$score
    
    a <- t.test(EGF, EGF)
    b <- t.test(EGF, ip70S6K)
    c <- t.test(EGF, iRSK)
    
    p_value <- append(p_value, a[["p.value"]])
    p_value <- append(p_value, b[["p.value"]])
    p_value <- append(p_value, c[["p.value"]])
  }
  
  
  df_temp <- df_temp %>% 
    group_by(category, cluster, treatment) %>% 
    summarise(fold_change = mean(score))
  
  df_temp$p_value <- p_value
  
  matrix_fold_changes <- rbind(matrix_fold_changes, df_temp)
}

matrix_fold_changes$minus_log10_p_value = -log10(matrix_fold_changes$p_value)
matrix_fold_changes$log2_fold_change = log2(matrix_fold_changes$fold_change)

# filter out treated groups
matrix_fold_changes <- matrix_fold_changes %>% 
  filter(treatment == "ip70S6K_EGF" | treatment == "iRSK_EGF")

matrix_fold_changes$category <- as_factor(matrix_fold_changes$category)
levels(matrix_fold_changes$category) <- categories

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_fold_changes_upon_treatment.pdf", sep = ""),
      height = 4,
      width = 8)
}

ggplot(data = matrix_fold_changes, 
       aes(y = category, 
           x = factor(cluster), 
           col = log2_fold_change, 
           size = minus_log10_p_value)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("#00008B", "#0000CC", "grey", "#FF0000", "#8B0000"),
                         values = scales::rescale(c(-1,0,1)),
                         limits = c(-1, 1)) +
  scale_size_continuous(limits = c(0,15)) +
  facet_wrap(facets = "treatment") +
  labs(title= "Cell state scores",
       size = "-log10 p-value",
       colour = "log2 fold change") +
  xlab("Cluster") +
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank())+
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom")


if (save_figure == "yes") {
  dev.off()
}

if (save_table == "yes") {
  write_tsv(matrix_fold_changes, path = paste(data_dir, "/fold_changes_cell_state_scores.tsv", sep = ""))
}