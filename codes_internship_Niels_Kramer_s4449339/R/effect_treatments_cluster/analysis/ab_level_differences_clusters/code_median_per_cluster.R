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
  dplyr::select(sample_id, plate_number, treatment, cluster, all_of(non_significant_abs))

TMM_s <- TMM_normalised %>% 
  dplyr::select(sample_id, plate_number, treatment, cluster, !all_of(non_significant_abs))

# transform data frame into long format
TMM_normalised_ns_long <- TMM_ns %>% 
  gather(antibody, intensity, 5:length(TMM_ns))

levels(TMM_normalised_ns_long$antibody) <- levels(TMM_normalised_ns_long)

matrix_z_score_ns <- tibble()

for (ab in sort(non_significant_abs)) {
  TMM_temp <- TMM_normalised_ns_long %>% 
    filter(antibody == ab)
  TMM_temp <- TMM_temp %>% 
    group_by(antibody, cluster) %>% 
    summarise(intensity = median(intensity))
  

  matrix_z_score_ns <- rbind(matrix_z_score_ns, TMM_temp)
}

TMM_normalised_s_long <- TMM_s %>% 
  gather(antibody, intensity, 5:length(TMM_s))

significant_abs <- as.vector(unique(TMM_normalised_s_long$antibody))

matrix_z_score_s <- tibble()

for (ab in sort(significant_abs)) {
  TMM_temp <- TMM_normalised_s_long %>% 
    filter(antibody == ab)
  TMM_temp <- TMM_temp %>% 
    group_by(antibody, cluster) %>% 
    summarise(intensity = median(intensity))
  
  
  matrix_z_score_s <- rbind(matrix_z_score_s, TMM_temp)
}

# transform matrices to wider for the correlation between antibodies and 
# transform in numeric matrices
matrix_z_score_ns <- pivot_wider(matrix_z_score_ns, 
                                 names_from = cluster,
                                 names_prefix = "C",
                                 values_from = intensity)

#ungroup
matrix_z_score_ns <- ungroup(matrix_z_score_ns)
matrix_z_score_ns$type <- "Cell state marker"

#dplyr::select all columns with median intensities per cluster
TMM_heatmap <- matrix_z_score_ns %>% 
  dplyr::select(-antibody, -type)

#rename row names to antibody name
rownames(TMM_heatmap) <- matrix_z_score_ns$antibody


# transform in matrix for heatmap calculations
TMM_heatmap <- data.matrix(TMM_heatmap)

#dplyr::select all columns with median intensities per cluster
TMM_type <- matrix_z_score_ns %>% 
  dplyr::select(type)

#rename row names to antibody name
rownames(TMM_type) <- matrix_z_score_ns$antibody


# transform in matrix for heatmap calculations
TMM_type <- as.matrix(TMM_type) 


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/heatmap_median_z_score_per_cluster_cell_state_markers.pdf", sep = ""),
      height = 5,
      width = 5)
}

#create heatmap
ht1 <- Heatmap(TMM_heatmap, 
               name = "Z-score normalised intensity",
               column_title = "Cell State Marker Antibody",
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               column_names_side = "top",
               column_names_gp = gpar(fontsize = 10),
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 10),
               show_row_names = FALSE,
               show_row_dend = TRUE,
               row_dend_side = "left",
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))

# get row names
index <- row_order(ht1)

# get indexes of names
names <- row.names(TMM_heatmap)

# create vector to later reorder factor antibodies on
order_ns_abs <- vector()

# loop through the indexes
for (i in index){
  order_ns_abs <- append(order_ns_abs, names[i])
}

hta <- Heatmap(TMM_type, 
               name = "Antibody type",
               column_title = "",
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               show_column_names = TRUE,
               column_names_side = "top",
               show_row_names = TRUE,
               show_row_dend = FALSE,
               row_dend_side = "right",
               heatmap_legend_param = list(legend_direction = "horizontal", 
                                           legend_width = unit(5, "cm")),
               col = palette_plots[7])


draw(ht1 + hta, heatmap_legend_side = "bottom")

if (save_figure == "yes") {
  # save figures to pdf
  dev.off()
}

# transform matrices to wider for the correlation between antibodies and 
# transform in numeric matrices
matrix_z_score_s <- pivot_wider(matrix_z_score_s, 
                                 names_from = cluster,
                                 names_prefix = "C",
                                 values_from = intensity)

#ungroup
matrix_z_score_s <- ungroup(matrix_z_score_s)
matrix_z_score_s$type <- "Signalling state marker"

#dplyr::select all columns with median intensities per cluster
TMM_heatmap <- matrix_z_score_s %>% 
  dplyr::select(-antibody, -type)

#rename row names to antibody name
rownames(TMM_heatmap) <- matrix_z_score_s$antibody

# transform in matrix for heatmap calculations
TMM_heatmap <- data.matrix(TMM_heatmap)

#dplyr::select all columns with median intensities per cluster
TMM_type <- matrix_z_score_s %>% 
  dplyr::select(type)

#rename row names to antibody name
rownames(TMM_type) <- matrix_z_score_s$antibody

# transform in matrix for heatmap calculations
TMM_type <- as.matrix(TMM_type) 


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/heatmap_median_z_score_per_cluster_signalling_state_markers.pdf", sep = ""),
      height = 10,
      width = 5)
}

#create heatmap
ht2 <- Heatmap(TMM_heatmap, 
               name = "Z-score normalised intensity",
               column_title = "Signalling State Marker Antibody",
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               column_names_side = "top",
               column_names_gp = gpar(fontsize = 10),
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               show_row_names = TRUE,
               show_row_dend = TRUE,
               row_dend_side = "left",
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))

# get row names
index <- row_order(ht2)

# get indexes of names
names <- row.names(TMM_heatmap)

# create vector to later reorder factor antibodies on
order_s_abs <- vector()

# loop through the indexes
for (i in index){
  order_s_abs <- append(order_s_abs, names[i])
}

htb <- Heatmap(TMM_type, 
               name = "Antibody type",
               column_title = "",
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               show_column_names = TRUE,
               column_names_side = "top",
               show_row_names = TRUE,
               show_row_dend = FALSE,
               row_dend_side = "right",
               heatmap_legend_param = list(legend_direction = "horizontal", 
                                           legend_width = unit(5, "cm")),
               col = palette_plots[8])


draw(ht2 + htb, heatmap_legend_side = "bottom")

if (save_figure == "yes") {
  # save figures to pdf
  dev.off()
}


## place the heatmaps in vertical order

# combine z-score matrices
matrix_z_score <- rbind(matrix_z_score_ns[1:(length(matrix_z_score_ns))], matrix_z_score_s[1:(length(matrix_z_score_s))])

#ungroup
matrix_z_score <- ungroup(matrix_z_score)

#dplyr::select all columns with median intensities per cluster
TMM_heatmap <- matrix_z_score %>% 
  dplyr::select(-antibody, -type)

#rename row names to antibody name
rownames(TMM_heatmap) <- matrix_z_score$antibody

# transform in matrix for heatmap calculations
TMM_heatmap <- data.matrix(TMM_heatmap)

#dplyr::select all columns with median intensities per cluster
TMM_type <- matrix_z_score %>% 
  dplyr::select(type)

#rename row names to antibody name
rownames(TMM_type) <- matrix_z_score$antibody

# transform in matrix for heatmap calculations
TMM_type <- as.matrix(TMM_type) 


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/heatmap_median_z_score_per_cluster.pdf", sep = ""),
      height = 15,
      width = 5.5)
}

#create heatmap
ht3 <- Heatmap(TMM_heatmap, 
               name = "Z-score normalised intensity",
               column_title = "Median z-score normalised antibody intensity",
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               column_names_side = "top",
               column_names_gp = gpar(fontsize = 10),
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               show_row_names = TRUE,
               show_row_dend = TRUE,
               row_dend_side = "left",
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")))

# get row names
index <- row_order(ht3)

# get indexes of names
names <- row.names(TMM_heatmap)

# create vector to later reorder factor antibodies on
order_abs <- vector()

# loop through the indexes
for (i in index){
  order_abs <- append(order_abs, names[i])
}

htc <- Heatmap(TMM_type, 
               name = "Antibody type",
               column_title = "",
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               show_column_names = TRUE,
               column_names_side = "top",
               show_row_names = TRUE,
               show_row_dend = FALSE,
               row_dend_side = "right",
               heatmap_legend_param = list(legend_direction = "horizontal", 
                                           legend_width = unit(5, "cm")),
               col = palette_plots[7:8])


draw(ht3 + htc, heatmap_legend_side = "bottom")

if (save_figure == "yes") {
  # save figures to pdf
  dev.off()
}

# transform matrices into long format in order to plot easily
matrix_z_score_ns <- matrix_z_score_ns %>% 
  pivot_longer(names_to = "cluster",values_to = "intensity", cols = 2:10)

matrix_z_score_s <- matrix_z_score_s %>% 
  pivot_longer(names_to = "cluster",values_to = "intensity", cols = 2:10)

matrix_z_score <- matrix_z_score %>% 
  pivot_longer(names_to = "cluster",values_to = "intensity", cols = 2:10)


# reorder antibody levels based on clustering in heatmap
matrix_z_score_ns$antibody <- factor(matrix_z_score_ns$antibody, levels = rev(order_ns_abs))
matrix_z_score_s$antibody <- factor(matrix_z_score_s$antibody, levels = rev(order_s_abs))
matrix_z_score$antibody <- factor(matrix_z_score$antibody, levels = rev(order_abs))


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_median_z_score_per_cluster_cell_state_markers_clusters_unordered.pdf", sep = ""),
      height = 5,
      width = 5)
}

a <- ggplot(data = matrix_z_score_ns, 
            aes(y = antibody, 
                x = cluster, 
                col = intensity)) + 
  geom_point(size = 4) + 
  scale_colour_gradientn(colours = c("#00008B", "#0000CC", "grey", "#FF0000", "#8B0000"),
                         limits = c(-2.5,2.5))+
  scale_size_continuous() +
  labs(title= "Cell State Marker Antibodies") +
  xlab("Cluster") +
  ylab("") +
  labs(col = "z-score \nnormalised \nintensity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

plot(a)

if (save_figure == "yes") {
  # save figures to pdf
  dev.off()
}


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_median_z_score_per_cluster_signalling_state_markers_clusters_unordered.pdf", sep = ""),
      height = 15,
      width = 5)
}

b <- ggplot(data = matrix_z_score_s, 
            aes(y = antibody, 
                x = cluster, 
                col = intensity)) + 
  geom_point(size = 4) + 
  scale_colour_gradientn(colours = c("#00008B", "#0000CC", "grey", "#FF0000", "#8B0000"),
                         limits = c(-2.5,2.5))+
  scale_size_continuous() +
  labs(title= "Signalling State Marker Antibodies") +
  xlab("Cluster") +
  ylab("") +
  labs(col = "z-score \nnormalised \nintensity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

plot(b)

if (save_figure == "yes") {
  # save figures to pdf
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_median_z_score_per_cluster_clusters_unordered.pdf", sep = ""),
      height = 18,
      width = 5)
}

plot_grid(a, b, ncol = 1, rel_heights = c(20, 49))

#message at the end
if (save_figure == "yes") {
  print("figures have been saved")
  # save figures to pdf
  dev.off()
}

legend <- get_legend(a)

if (save_table == "yes") {
  write_tsv(matrix_z_score, path=paste(data_dir, "/median_z_scores_per_ab_per_cluster.tsv", sep = ""))
  write_tsv(matrix_z_score_s, path=paste(data_dir, "/median_z_scores_signalling_markers_per_ab_per_cluster.tsv", sep = ""))
  write_tsv(matrix_z_score_ns, path=paste(data_dir, "/median_z_scores_cell_state_markers_per_ab_per_cluster.tsv", sep = ""))
  }
