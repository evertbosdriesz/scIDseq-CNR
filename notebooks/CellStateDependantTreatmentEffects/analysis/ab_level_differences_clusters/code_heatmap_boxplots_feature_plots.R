## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL


#load dataframe that created list with shared non chaninging ABs

non_significant_abs <- read_csv(file = abs_no_fold_change)

# create factor out of the antibody column
non_significant_abs <- as.factor(non_significant_abs$antibody)


# create new data frame to only use the non significant the significant abs
TMM_ns <- TMM_normalised %>%
  dplyr::select(all_of(non_significant_abs))

# make a matrix where you can compare each cell with all other cells

TMM_ns <- as.matrix(TMM_ns)



## Do UMAP and add to the visualisation data set
# set seed for umap (seed for leiden is not via R)
set.seed(1)

# # Run umap
umap_out <- umap(TMM_ns)


umap_1 <- umap_out$layout[,1]
umap_2 <- umap_out$layout[,2]

TMM_normalised$umap_1 <- umap_1
TMM_normalised$umap_2 <- umap_2

# check if treatments are distributed equally over all clusters

TMM_normalised$cluster <- as.factor(TMM_normalised$cluster)

df_check <- TMM_normalised %>% 
  group_by(cluster, treatment) %>% 
  summarise(count = n())

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/information_concencus_cluster.pdf", sep = ""),
      width=20, height=10)
}

ggplot(data = df_check, aes(x = cluster, y = count, fill = treatment)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Amount of treatment per cluster") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  ylab("Fraction of cells") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = palette_plots)


ggplot(TMM_normalised, aes(x = cluster, fill = treatment)) +
  geom_bar() +
  scale_fill_manual(values = palette_plots) +
  ggtitle("Amount of treatment per cluster") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  ylab("Fraction of cells") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())

#FOR clusters on umap

#create scatter plot of umap
umap_cluster <- ggplot(TMM_normalised) + 
    geom_point(aes(x=umap_1, y=umap_2, color=cluster), size = 5) + 
    ggtitle("UMAP on Cell State Markers") +
    labs(subtitle="Clustered on cell state marker antibodies") + 
    theme(plot.title = element_text(size=12, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    scale_colour_manual(values = palette_plots)

print(umap_cluster)
  
#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=treatment), size = 5) + 
  ggtitle("UMAP on Cell State Markers") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)

#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=plate_number), size = 5) + 
  ggtitle("UMAP on Cell State Markers") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)

# create new data frame to only use the non significant the significant abs
TMM_normalised_not_significant <- TMM_normalised %>%
  dplyr::select(sample_id, plate_number, treatment, umap_1, umap_2, cluster, all_of(non_significant_abs))

TMM_normalised_significant <- TMM_normalised %>% 
  dplyr::select(sample_id, plate_number, treatment, umap_1, umap_2, !all_of(non_significant_abs))

TMM_heatmap <- TMM_normalised_not_significant %>% 
  dplyr::select(!one_of(c("cluster","plate_number", "treatment", "sample_id", "umap_1", "umap_2")))

rownames(TMM_heatmap) <- TMM_normalised$sample_id

TMM_heatmap <- data.matrix(TMM_heatmap)


# create heatmap

# set the colour palette for the different plots
palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")


ht1 <- Heatmap(TMM_heatmap, 
               name = "Intensity",
               column_title = "Cell State Marker Antibody",
               split = TMM_normalised$cluster,
               cluster_columns = TRUE,
               cluster_rows = TRUE,
               column_names_side = "top",
               column_names_gp = gpar(fontsize = 10),
               show_row_names = FALSE,
               show_row_dend = FALSE,
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))


draw(ht1, heatmap_legend_side = "bottom")


#create heatmap of significant genes
TMM_heatmap <- TMM_normalised_significant %>% 
  dplyr::select(!one_of(c("cluster","plate_number", "treatment", "sample_id", "umap_1", "umap_2")))

rownames(TMM_heatmap) <- TMM_normalised$sample_id

TMM_heatmap <- data.matrix(TMM_heatmap)

ht2 <- Heatmap(TMM_heatmap, 
             name = "Intensity",
             column_title = "Signalling Marker Antibody",
             split = TMM_normalised$cluster,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             column_names_side = "top",
             column_names_gp = gpar(fontsize = 10),
             show_row_names = FALSE,
             show_row_dend = FALSE,
             heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))

draw(ht2, heatmap_legend_side = "bottom")

draw(ht1 + ht2, heatmap_legend_side = "bottom")


# transform data frame into long format
TMM_normalised_ns_long <- TMM_normalised_not_significant %>% 
  gather(antibody, intensity, 7:length(TMM_normalised_not_significant))

levels(TMM_normalised_ns_long$antibody) <- levels(TMM_normalised_ns_long)

for (ab in sort(non_significant_abs)) {
  TMM_temp <- TMM_normalised_ns_long %>% 
    filter(antibody == ab)
    
    # What is the effect of the treatment on the value ?
    model=lm( TMM_temp$intensity ~ TMM_temp$cluster )
    ANOVA=aov(model)
    
    # Tukey test to study each pair of treatment :
    TUKEY <- TukeyHSD(x=ANOVA, 'TMM_temp$cluster', conf.level=0.95)
    
    LABELS <- generate_label_df(TUKEY , "TMM_temp$cluster")
  
  
    plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "intensity", color = "cluster", 
                         add = "jitter", legend = "none") +
      geom_text(data = LABELS, aes(x = treatment, y = max(TMM_temp$intensity) + 0.1 * max(TMM_temp$intensity)), label = LABELS$Letters) +
      rotate_x_text(angle = 45)+
      theme_light() +
      labs(title = "Normalised antibody intensities per cluster",
           caption = "antibody used for clustering") +
      ylab("z-score normalised antibody intensity") +
      facet_wrap(facets = vars(antibody), ncol = 1) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme_light() +
      theme(plot.title = element_text(size=12, hjust = 0.5),
            plot.subtitle = element_text(size=10),
            plot.caption = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            axis.line = element_line(colour = "black"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank()) +
      theme(strip.text.x = element_text(size = 10, color = "black")) +
      scale_colour_manual(values = palette_plots)
  
  plot_b <- #create scatter plot of umap
    ggplot(TMM_temp) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity), size = 5) + 
    ggtitle("UMAP on Cell State Markers") +
    labs(caption = "antibody used for clustering") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_text(size=12, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 10, color = "black")) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) - 0.05,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))),
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) + 0.05,
                                                     1))        
  
  plot_list <- list(plot_a, plot_b, umap_cluster)
  
  print(plot_grid(plotlist = plot_list))
}

TMM_normalised_s_long <- TMM_normalised_significant %>% 
  gather(antibody, intensity, 7:length(TMM_normalised_significant))

significant_abs <- as.vector(unique(TMM_normalised_s_long$antibody))

for (ab in sort(significant_abs)) {
  TMM_temp <- TMM_normalised_s_long %>% 
    filter(antibody == ab)
  
  # What is the effect of the treatment on the value ?
  model=lm( TMM_temp$intensity ~ TMM_temp$cluster )
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'TMM_temp$cluster', conf.level=0.95)
  
  LABELS <- generate_label_df(TUKEY , "TMM_temp$cluster")
  
  plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "intensity", color = "cluster", 
                      add = "jitter", legend = "none") +
    geom_text(data = LABELS, aes(x = treatment, y = max(TMM_temp$intensity) + 0.1 * max(TMM_temp$intensity)), label = LABELS$Letters) +
        rotate_x_text(angle = 45)+
    theme_light() +
    labs(title = "Normalised antibody intensities per cluster",
         caption = "antibody has not been used for clustering") +
    ylab("z-score normalised antibody intensity") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_text(size=12, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 10, color = "black")) +
    scale_colour_manual(values = palette_plots)
  
  plot_c <- ggboxplot(TMM_temp, x = "cluster", y = "intensity", color = "treatment", 
                      add = "jitter", legend = "none") +
    rotate_x_text(angle = 45)+
    theme_light() +
    labs(title = "Normalised antibody intensities per cluster",
         caption = "antibody has not been used for clustering") +
    ylab("z-score normalised antibody intensity") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_text(size=12, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 10, color = "black")) +
    scale_colour_manual(values = palette_plots)
  
  plot_b <- #create scatter plot of umap
    ggplot(TMM_temp) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity), size = 5) + 
    ggtitle("UMAP on Cell State Markers") +
    labs(caption = "antibody has not been used for clustering") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_text(size=12, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 10, color = "black")) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) - 0.05,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))),
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) + 0.05,
                                                     1))        
  
  plot_list <- list(plot_a, plot_c, plot_b, umap_cluster)
  
  print(plot_grid(plotlist = plot_list))
}

#message at the end

if (save_figure == "yes") {
  print("figures have been saved")
  # save figures to pdf
  dev.off()
}
