## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment1 | treatment == "EGF")

#load dataframe that created list with shared non chaninging ABs

non_significant_abs <- read_csv(file = abs_no_fold_change)

# create factor out of the antibody column
non_significant_abs <- as.factor(non_significant_abs$antibody)


# create new data frame to only use the non significant the significant abs
TMM_ns <- TMM_normalised %>%
  select(sample_id, plate_number, treatment, cluster, all_of(non_significant_abs))

TMM_s <- TMM_normalised %>% 
  select(sample_id, plate_number, treatment, cluster, !all_of(non_significant_abs))

# transform data frame into long format
TMM_normalised_ns_long <- TMM_ns %>% 
  gather(antibody, intensity, 5:length(TMM_ns))

levels(TMM_normalised_ns_long$antibody) <- levels(TMM_normalised_ns_long)

matrix_fold_changes_ns <- tibble()

for (ab in sort(non_significant_abs)) {
  TMM_temp <- TMM_normalised_ns_long %>% 
    filter(antibody == ab)
  
  p_value <- vector()
  
  for (i in c(1:length(unique(TMM_temp$cluster)))) {
    
    TMM_temp_cluster <- TMM_temp %>% 
      filter(cluster == i)
    
    EGF <- TMM_temp_cluster %>% 
      filter(treatment == "EGF")
    
    EGF <- EGF$intensity
    
    inhibitor <- TMM_temp_cluster %>% 
      filter(treatment == treatment1)
    
    inhibitor <- inhibitor$intensity
    
    
    a <- t.test(EGF, EGF)
    b <- t.test(EGF, inhibitor)
   
    
    p_value <- append(p_value, a[["p.value"]])
    p_value <- append(p_value, b[["p.value"]])
  }
  
  
  TMM_temp <- TMM_temp %>% 
    group_by(antibody, cluster, treatment) %>% 
    summarise(fold_change = mean(intensity))
  
  TMM_temp$p_value <- p_value
  
  matrix_fold_changes_ns <- rbind(matrix_fold_changes_ns, TMM_temp)
}

TMM_normalised_s_long <- TMM_s %>% 
  gather(antibody, intensity, 5:length(TMM_s))

significant_abs <- as.vector(unique(TMM_normalised_s_long$antibody))

matrix_fold_changes_s <- tibble()

for (ab in sort(significant_abs)) {
  TMM_temp <- TMM_normalised_s_long %>% 
    filter(antibody == ab)
  
  p_value <- vector()
  
  for (i in c(1:length(unique(TMM_temp$cluster)))) {
    
    TMM_temp_cluster <- TMM_temp %>% 
      filter(cluster == i)
    
    EGF <- TMM_temp_cluster %>% 
      filter(treatment == "EGF")
    
    EGF <- EGF$intensity
    
    inhibitor <- TMM_temp_cluster %>% 
      filter(treatment == treatment1)
    
    inhibitor <- inhibitor$intensity
    
    
    a <- t.test(EGF, EGF)
    b <- t.test(EGF, inhibitor)
    
    
    p_value <- append(p_value, a[["p.value"]])
    p_value <- append(p_value, b[["p.value"]])
  }
  
  
  TMM_temp <- TMM_temp %>% 
    group_by(antibody, cluster, treatment) %>% 
    summarise(fold_change = mean(intensity))
  
  TMM_temp$p_value <- p_value
  
  matrix_fold_changes_s <- rbind(matrix_fold_changes_s, TMM_temp)
}

matrix_fold_changes_ns$minus_log10_p_value = -log10(matrix_fold_changes_ns$p_value)
matrix_fold_changes_ns$log2_fold_change = log2(matrix_fold_changes_ns$fold_change)

matrix_fold_changes_ns <- matrix_fold_changes_ns %>% 
  filter(treatment == treatment1)

# open the order for the factors antibody and cluster
order_abs_ns <- read_tsv("C:/Users/niels/Documents/R/01_data_in_use/order_ns_abs.tsv")

#transform into vector
order_abs_ns <- as_vector(order_abs_ns$order_ns_abs)

# transform antibody and cluster into factor with right levels
matrix_fold_changes_ns$antibody <- factor(matrix_fold_changes_ns$antibody, levels = rev(order_abs_ns))


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_fold_changes_per_treatment_cell_state_markers_for_", treatment1, ".pdf", sep = ""),
      height = 5,
      width = 4)
}

a <- ggplot(data = matrix_fold_changes_ns, 
       aes(y = antibody, 
           x = factor(cluster), 
           col = log2_fold_change, 
           size = minus_log10_p_value)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("#00008B", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "#8B0000"),
                         values = rescale(c(-1,0,1)),
                         limits = c(-2, 2)) +
  scale_size_continuous(limits = c(0,15)) +
  facet_wrap(facets = "treatment") +
  labs(title= "Cell State Marker Antibodies",
       size = "-log10 p-value",
       colour = "log2 fold change") +
  xlab("Cluster") +
  ylab("") +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
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
  rremove("legend")

plot(a)

if (save_figure == "yes") {
  dev.off()
}

matrix_fold_changes_s$minus_log10_p_value = -log10(matrix_fold_changes_s$p_value)
matrix_fold_changes_s$log2_fold_change = log2(matrix_fold_changes_s$fold_change)


matrix_fold_changes_s <- matrix_fold_changes_s %>% 
  filter(treatment == treatment1)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_fold_changes_per_treatment_signalling_state_markers_for_", treatment1, ".pdf", sep = ""),
      height = 15,
      width = 4)
}

# open the order for the factors antibody and cluster
order_abs_s <- read_tsv("C:/Users/niels/Documents/R/01_data_in_use/order_s_abs.tsv")

#transform into vector
order_abs_s <- as_vector(order_abs_s$order_s_abs)

# transform antibody and cluster into factor with right levels
matrix_fold_changes_s$antibody <- factor(matrix_fold_changes_s$antibody, levels = rev(order_abs_s))

b <- ggplot(data = matrix_fold_changes_s, 
       aes(y = antibody, 
           x = factor(cluster), 
           col = log2_fold_change, 
           size = minus_log10_p_value)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("#00008B", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "#8B0000"),
                         values = rescale(c(-1,0,1)),
                         limits = c(-2, 2)) +
  scale_size_continuous(limits = c(0,15)) +
  facet_wrap(facets = "treatment") +
  labs(title= "Signalling State Marker Antibodies",
       size = "-log10 p-value",
       colour = "log2 fold change") +
  xlab("Cluster") +
  ylab("") +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
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
  rremove("legend")

plot(b)

if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/dot_plots_fold_changes_per_treatment_per_cluster_for_", treatment1, ".pdf", sep = ""),
      height = 15,
      width = 4)
}

plot_grid(a, b, ncol = 1, rel_heights = c(20, 49))

#message at the end
if (save_figure == "yes") {
  dev.off()
}

matrix_fold_changes_ns$type <- "signalling_state_marker"
matrix_fold_changes_s$type <- "cell_state_marker"

matrix_fold_changes <- rbind(matrix_fold_changes_s, matrix_fold_changes_ns)

if (save_table == "yes") {
  write_tsv(matrix_fold_changes, path = paste(data_dir, "/matrix_fold_changes_iRSK_and_ip70S6K.tsv", sep = ""))
  print("table has been saved")
}

if (save_figure == "yes") {
  print("figures have been saved")
  # save figures to pdf
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/legend_dot_plots_fold_changes_per_treatment_cell_state_markers_for_", treatment1, ".pdf", sep = ""),
      height = 5,
      width = 8)
}

c <- ggplot(data = matrix_fold_changes_ns, 
            aes(y = antibody, 
                x = factor(cluster), 
                col = log2_fold_change, 
                size = minus_log10_p_value)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("#00008B", "#00008B", "#0000CC", "grey", "#FF0000", "#8B0000", "#8B0000"),
                         values = rescale(c(-1,0,1)),
                         limits = c(-2, 2)) +
  scale_size_continuous(limits = c(0,12)) +
  facet_wrap(facets = "treatment") +
  labs(title= "Cell State Marker Antibodies",
       size = "-log10 p-value",
       colour = "log2 fold change") +
  xlab("Cluster") +
  ylab("") +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
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

plot(c)

#message at the end
if (save_figure == "yes") {
  dev.off()
}

if (save_table == "yes") {
  write_tsv(matrix_fold_changes_ns, paste(data_dir, "/cluster_fold_changes_for_cell_state_markers_for_", treatment1, ".tsv", sep = ""))
  write_tsv(matrix_fold_changes_s, paste(data_dir, "/cluster_fold_changes_for_signalling_state_markers_for_", treatment1, ".tsv", sep = ""))
}
