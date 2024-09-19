## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# filter on treated abs (exclude no EGF)

TMM_normalised <- TMM_normalised %>%
  filter(treatment == treatment_1 | treatment == treatment_2 | treatment == treatment_3)

# open the out_put file to add later add the new clusters to
TMM_normalised_clusters <- TMM_normalised[1:3]

# open the final output file to add the weighed clusters to
TMM_final_cluster <- TMM_normalised

#place sample_id as row names
TMM_normalised <- TMM_normalised %>%
  remove_rownames %>%
  column_to_rownames("sample_id")

#load dataframe that created list with shared non chaninging ABs

non_significant_abs <- read_csv(file = abs_no_fold_change)

# create factor out of the antibody column
non_significant_abs <- as.factor(non_significant_abs$antibody)


# create new data frame to only use the non significant the significant abs
TMM_normalised <- TMM_normalised %>%
  select(all_of(non_significant_abs))

# make a matrix where you can compare each cell with all other cells

TMM_normalised <- as.matrix(TMM_normalised)


# transpose matrix for the nearest neighbour calculation
cor_matrix <- t(TMM_normalised)

# calculate 30 nearest neighbours of cell
snn <- RANN::nn2(t(cor_matrix), k = no_of_nearest_neighbours)$nn.idx

# initialise matrix where all values are 0 (meaning no connection between cells)
adjacency_matrix <- matrix(0L, ncol(cor_matrix), ncol(cor_matrix))


# rename columns and rows to correspond to sample_id
rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- colnames(cor_matrix)

# replace value in empty adjacency matrix for the nearest neighbours to 1 (meaning a connection)
for(ii in 1:ncol(cor_matrix)) {
  adjacency_matrix[ii,colnames(cor_matrix)[snn[ii,]]] <- 1L
}

# check that rows add to k
sum(adjacency_matrix[1,]) == no_of_nearest_neighbours
table(apply(adjacency_matrix, 1, sum))

for (i in c(1:no_repeats)) {
  partition <- leiden(adjacency_matrix, resolution_parameter = res, 
                      n_iterations = no_iter)
  
  # count cells per group
  print(table(partition))
  new_column <- paste("cluster_round", i, sep = "_")
  TMM_normalised_clusters[,new_column] <- partition
} 

no_of_clusters <- c()

for (i in c(1:no_repeats)){
  no_of_clusters <- append(no_of_clusters, 
         nrow(unique(TMM_normalised_clusters[,paste("cluster_round", i, sep = "_")])))
}

sort(no_of_clusters)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/information_clusters_after_", no_repeats, "_repeats_of_leiden_filtered.pdf", sep = ""),
      width=20, height=10)
}

ggplot(data = tibble(no_of_clusters), aes(x = no_of_clusters)) +
  geom_histogram(binwidth = 1, fill = palette_plots[1:length(unique(no_of_clusters))]) +
  labs(title = paste("Frequency of amount of clusters in ", no_repeats, " repeats of leiden", sep = "")) +
  xlab("number of clusters") +
  theme_light() +
  my_theme +
  theme(axis.text.x = element_text(angle = 0))

no_of_clusters <- as.numeric(names(sort(table(no_of_clusters),decreasing=TRUE)[1]))

# delete the cluster columns that don't contain the right number of clusters
for (i in c(1:no_repeats)){
  column <- paste("cluster_round", i, sep = "_")
  if (nrow(unique(TMM_normalised_clusters[,column])) != no_of_clusters){
    TMM_normalised_clusters[,column] <- NULL
  }
}



if (save_table == "yes") {
  write.csv(TMM_normalised_clusters, paste(data_dir, "/clusters_after_", no_repeats, "_repeats_of_leiden_filtered.csv", sep = "")) 
}

rownames(TMM_normalised_clusters) <- TMM_normalised_clusters$sample_id

TMM_clusters <- TMM_normalised_clusters %>% 
  select(!one_of(c("plate_number", "treatment", "sample_id")))



TMM_clusters <- data.matrix(TMM_clusters)

# make heatmap of cluster distribution over Leiden iterations
Heatmap(TMM_clusters, 
        col = rev(rainbow(max(no_of_clusters))),
        name = "cluster number",
        column_title = "cluster round",
        show_row_names = FALSE,
        show_column_names = FALSE
)

# calculate heatmap and clusters with cooccurences
cor_clusters <- cooccur(TMM_clusters)

Heatmap(cor_clusters, 
        name = "co-occurence distance",
        column_title = "Correlation between assigned clusters per cell based on co-occurences",
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE
)


dist_clusters <- hclust(dist(cor_clusters))

clusters_weighed <- cutree(dist_clusters, k = no_of_clusters)

TMM_final_cluster$cluster <- clusters_weighed





TMM_final_cluster <- TMM_final_cluster %>% 
  select(sample_id, plate_number, treatment, cluster, everything())


if (save_table == "yes") {
  write.csv(TMM_final_cluster, paste(data_dir, "/TMM_normalised_weighed_clusters_after_", no_repeats, "_repeats_of_leiden_filtered.csv", sep = "")) 
}



## Do UMAP and add to the visualisation data set
# set seed for umap (seed for leiden is not via R)
set.seed(1)

# # Run umap
umap_out <- umap(TMM_normalised)


umap_1 <- umap_out$layout[,1]
umap_2 <- umap_out$layout[,2]

TMM_final_cluster$umap_1 <- umap_1
TMM_final_cluster$umap_2 <- umap_2

# check if treatments are distributed equally over all clusters

TMM_final_cluster$cluster <- as.factor(TMM_final_cluster$cluster)

df_check <- TMM_final_cluster %>% 
  group_by(cluster, treatment) %>% 
  summarise(count = n())

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


ggplot(TMM_final_cluster, aes(x = cluster, fill = treatment)) +
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
umap_cluster <- ggplot(TMM_final_cluster) + 
    geom_point(aes(x=umap_1, y=umap_2, color=cluster), size = 5) + 
    ggtitle("UMAP on Cell State Markers") +
    labs(subtitle="Clustered on cell state marker antibodies") + 
    scale_x_continuous(breaks = seq(0, no_of_clusters, 1)) +
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
ggplot(TMM_final_cluster) + 
  geom_point(aes(x=umap_1, y=umap_2, color=treatment), size = 5) + 
  ggtitle("UMAP on Cell State Markers") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  scale_x_continuous(breaks = seq(0, no_of_clusters, 1)) +
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
ggplot(TMM_final_cluster) + 
  geom_point(aes(x=umap_1, y=umap_2, color=plate_number), size = 5) + 
  ggtitle("UMAP on Cell State Markers") +
  labs(subtitle="Clustered on cell state marker antibodies") + 
  scale_x_continuous(breaks = seq(0, no_of_clusters, 1)) +
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
TMM_normalised_not_significant <- TMM_final_cluster %>%
  select(sample_id, plate_number, treatment, umap_1, umap_2, cluster, all_of(non_significant_abs))

TMM_normalised_significant <- TMM_final_cluster %>% 
  select(sample_id, plate_number, treatment, umap_1, umap_2, !all_of(non_significant_abs))

TMM_heatmap <- TMM_normalised_not_significant %>% 
  select(!one_of(c("cluster","plate_number", "treatment", "sample_id", "umap_1", "umap_2")))

rownames(TMM_heatmap) <- TMM_final_cluster$sample_id

TMM_heatmap <- data.matrix(TMM_heatmap)


# create heatmap

# set the colour palette for the different plots
palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

palette_right_order <- c(palette_plots[8], palette_plots[5], palette_plots[2], palette_plots[3], palette_plots[1], palette_plots[7], palette_plots[4], palette_plots[6])

ht1 <- Heatmap(TMM_heatmap, 
               name = "Intensity",
               column_title = "Cell State Marker Antibody",
               split = TMM_final_cluster$cluster,
               cluster_rows = TRUE,
               column_names_side = "top",
               column_names_gp = gpar(fontsize = 10),
               show_row_names = FALSE,
               show_row_dend = FALSE,
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))


draw(ht1, heatmap_legend_side = "bottom")


#create heatmap of significant genes
TMM_heatmap <- TMM_normalised_significant %>% 
  select(!one_of(c("cluster","plate_number", "treatment", "sample_id", "umap_1", "umap_2")))

rownames(TMM_heatmap) <- TMM_final_cluster$sample_id

TMM_heatmap <- data.matrix(TMM_heatmap)

ht2 <- Heatmap(TMM_heatmap, 
             name = "Intensity",
             column_title = "Signalling Marker Antibody",
             split = TMM_final_cluster$cluster,
             cluster_rows = TRUE,
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
  
    plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "intensity", color = "cluster", 
                         add = "jitter", legend = "none") +
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
  
  plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "intensity", color = "cluster", 
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
