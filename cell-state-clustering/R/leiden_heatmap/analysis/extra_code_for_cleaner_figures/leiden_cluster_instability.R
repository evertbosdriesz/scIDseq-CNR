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
  pdf(file=paste(figure_dir, "/different_cluster_nos_after_leiden.pdf", sep = ""),
      width=3, height=2)
}

ggplot(data = tibble(no_of_clusters), aes(x = factor(no_of_clusters))) +
  geom_histogram(binwidth = 1, fill = palette_plots[1:length(unique(no_of_clusters))]) +
  xlab("Number of clusters") +
  ylab("Count") +
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
        strip.background = element_blank())

if (save_figure == "yes") {
  dev.off()
}

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

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/heatmap_cluster_distribution.pdf", sep = ""),
      width=5, height=3.5)
}

# make heatmap of cluster distribution over Leiden iterations
ht1 <- Heatmap(TMM_clusters, 
        col = rev(rainbow(max(no_of_clusters))),
        name = "Cluster",
        column_title = "",
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(legend_direction = "horizontal"))

draw(ht1)

ht1 <- Heatmap(TMM_clusters, 
               col = rev(rainbow(max(no_of_clusters))),
               name = "Cluster",
               column_title = "",
               show_row_names = FALSE,
               show_column_names = FALSE,
               show_heatmap_legend = FALSE)

draw(ht1)

if (save_figure == "yes") {
  dev.off()
}