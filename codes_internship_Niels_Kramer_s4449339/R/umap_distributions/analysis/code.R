#load libraries
library(tidyverse)
library(yarrr)

#load in data frame

TMM_normalised <- read_tsv("R/01_data_in_use/TMM_normalised_z_transformed_clusters_with_umap.tsv")

TMM_normalised$cluster <- as_factor(TMM_normalised$cluster)




#FOR clusters on umap

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_clusters_with_legend.pdf", sep = ""),
      width=3, height=3)
}

#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=cluster), size = 3) + 
  ggtitle("") +
  labs(subtitle="", col = "cluster") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +  
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_clusters.pdf", sep = ""),
      width=3, height=3)
}

#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=cluster), size = 1, show.legend = FALSE) + 
  ggtitle("") +
  labs(subtitle="", col = "cluster") + 
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}

# for treatments on UMAP

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_treatment_with_legend.pdf", sep = ""),
      width=3, height=3)
}

#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=treatment), size = 3) + 
  ggtitle("") +
  labs(subtitle="", col = "treatment") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +  
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_treatments.pdf", sep = ""),
      width=3, height=3)
}

#create scatter plot of umap
ggplot(TMM_normalised) + 
  geom_point(aes(x=umap_1, y=umap_2, color=treatment), size = 1, show.legend = FALSE) + 
  ggtitle("") +
  labs(subtitle="", col = "treatment") + 
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/feature_plots_with_legend.pdf", sep = ""),
      width=3, height=3)
}

TMM_normalised_long <- TMM_normalised %>% 
  pivot_longer(names_to = "antibody", values_to = "intensity", 7:length(TMM_normalised))

for(ab in sort(unique(TMM_normalised_long$antibody))){
  TMM_temp <- TMM_normalised_long %>% 
    filter(antibody == ab)
  
  plot <- ggplot(TMM_temp) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity)) + 
    facet_wrap(facets = "antibody", ncol = 4) +
    labs(fill = "Intensity") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
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
    # theme(plot.title = element_text(size=10),
    #       plot.subtitle = element_text(size=10),
    #       legend.text = element_text(size = 8),
    #       legend.title = element_text(size = 10),
    #       panel.border = element_blank(), 
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.background = element_blank(),
    #       axis.line = element_line(colour = "black")) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) - 0.05,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))),
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) + 0.05,
                                                     1))
  
  print(plot)
}


if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/feature_plots.pdf", sep = ""),
      width=3, height=3)
}

for(ab in sort(unique(TMM_normalised_long$antibody))){
  TMM_temp <- TMM_normalised_long %>% 
    filter(antibody == ab)
  
  plot <- ggplot(TMM_temp) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity), show.legend = FALSE) + 
    facet_wrap(facets = "antibody", ncol = 4) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
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
    # theme(plot.title = element_text(size=10),
    #       plot.subtitle = element_text(size=10),
    #       legend.text = element_text(size = 8),
    #       legend.title = element_text(size = 10),
    #       panel.border = element_blank(), 
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.background = element_blank(),
    #       axis.line = element_line(colour = "black")) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) - 0.05,
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))),
                                                     1 - (max(TMM_temp$intensity) * 1 / (max(TMM_temp$intensity) - min(TMM_temp$intensity))) + 0.05,
                                                     1))
  
  print(plot)
}



if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/distribution_graph.pdf", sep = ""),
      width=5, height=3)
}

# check distrubtion of cells over clusters
df_check <- TMM_normalised %>% 
  group_by(cluster, treatment) %>% 
  summarise(count = n())

# create bar graph
ggplot(TMM_normalised, aes(x = factor(cluster), fill = treatment)) +
  geom_bar() +
  scale_fill_manual(values = palette_plots) +
  ylab("cell count") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
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

# create bar graph
ggplot(TMM_normalised, aes(x = factor(cluster), fill = treatment)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = palette_plots) +
  ylab("Fraction of cells") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
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

if (save_figure == "yes") {
  dev.off()
}

