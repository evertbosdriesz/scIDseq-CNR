## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# filter on treatment of interest
TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment1)


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

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/box_plots_", treatment1, "_normalised_per_cluster_on_EGF_log2_fold_change.pdf", sep = ""),
      width=4, height=3)
}


TMM_normalised_ns_long$log2_fold_change = log2(TMM_normalised_ns_long$intensity)


for (ab in sort(non_significant_abs)) {
  TMM_temp <- TMM_normalised_ns_long %>% 
    filter(antibody == ab)
    
  TMM_temp$cluster <- as_factor(TMM_temp$cluster)
  
    # What is the effect of the treatment on the value ?
    model=lm( TMM_temp$log2_fold_change ~ TMM_temp$cluster )
    ANOVA=aov(model)
    
    # Tukey test to study each pair of treatment :
    TUKEY <- TukeyHSD(x=ANOVA, 'TMM_temp$cluster', conf.level=0.95)
    
    LABELS <- generate_label_df(TUKEY , "TMM_temp$cluster")
  
  
    plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "log2_fold_change", color = "treatment", 
                        add = "jitter", size = 0.2,
                        add.params = list(size = 1.6)) +
      rotate_x_text(angle = 45)+
      theme_light() +
      labs(title = "Normalised antibody intensities per cluster",
           caption = "antibody has not been used for clustering") +
      ylab("log2 fold change") +
      xlab("Cluster") +
      facet_wrap(facets = vars(antibody), ncol = 1) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme_light() +
      theme(plot.title = element_blank(),
            plot.subtitle = element_text(size=10),
            plot.caption = element_blank(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9, colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank()) +
      theme(strip.text.x = element_text(size = 10, color = "black")) +
      scale_colour_manual(values = palette_plots[correct_colour]) +
      theme(legend.position = "bottom")

    # mirrors the y-axis when pos y is bigger
    if (max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change) >= abs(min(TMM_temp$log2_fold_change))) {
      plot_a <- plot_a + coord_cartesian(ylim = c(-(max(TMM_temp$log2_fold_change) + 
                                                      0.1 * max(TMM_temp$log2_fold_change)), 
                                                  max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change))) +
        geom_text(data = LABELS, aes(x = treatment, y = max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change)), label = LABELS$Letters)
    } 
    # mirrors the y-axis when neg y is bigger 
    if (max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change) <= abs(min(TMM_temp$log2_fold_change))) {
      plot_a <- plot_a + coord_cartesian(ylim = c(min(TMM_temp$log2_fold_change), 
                                                  abs(min(TMM_temp$log2_fold_change)))) +
        geom_text(data = LABELS, aes(x = treatment, y = abs(min(TMM_temp$log2_fold_change))), label = LABELS$Letters)
      
    }
  
    plot(plot_a)
}

TMM_normalised_s_long <- TMM_s %>% 
  gather(antibody, intensity, 5:length(TMM_s))

TMM_normalised_s_long$log2_fold_change = log2(TMM_normalised_s_long$intensity)


significant_abs <- as.vector(unique(TMM_normalised_s_long$antibody))

for (ab in sort(significant_abs)) {
  TMM_temp <- TMM_normalised_s_long %>% 
    filter(antibody == ab)
  
  TMM_temp$cluster <- as.factor(TMM_temp$cluster)
  
  # What is the effect of the treatment on the value ?
  model=lm( TMM_temp$log2_fold_change ~ TMM_temp$cluster )
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'TMM_temp$cluster', conf.level=0.95)
  
  LABELS <- generate_label_df(TUKEY , "TMM_temp$cluster")
  
  plot_a <- ggboxplot(TMM_temp, x = "cluster", y = "log2_fold_change", color = "treatment", 
                      add = "jitter", size = 0.2,
                      add.params = list(size = 1.6)) +
    rotate_x_text(angle = 45)+
    theme_light() +
    labs(title = "Normalised antibody intensities per cluster",
         caption = "antibody has not been used for clustering") +
    ylab("log2 fold change") +
    xlab("Cluster") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(size=10),
          plot.caption = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9, colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 10, color = "black")) +
    scale_colour_manual(values = palette_plots[correct_colour]) +
    theme(legend.position = "bottom")
  
  # mirrors the y-axis when pos y is bigger
  if (max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change) >= abs(min(TMM_temp$log2_fold_change))) {
    plot_a <- plot_a + coord_cartesian(ylim = c(-(max(TMM_temp$log2_fold_change) + 
                                                    0.1 * max(TMM_temp$log2_fold_change)), 
                                                max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change))) +
      geom_text(data = LABELS, aes(x = treatment, y = max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change)), label = LABELS$Letters)
  } 
  # mirrors the y-axis when neg y is bigger 
  if (max(TMM_temp$log2_fold_change) + 0.1 * max(TMM_temp$log2_fold_change) <= abs(min(TMM_temp$log2_fold_change))) {
    plot_a <- plot_a + coord_cartesian(ylim = c(min(TMM_temp$log2_fold_change), 
                                                abs(min(TMM_temp$log2_fold_change)))) +
      geom_text(data = LABELS, aes(x = treatment, y = abs(min(TMM_temp$log2_fold_change))), label = LABELS$Letters)
    
  }
     
    
  
  plot(plot_a)
}

#message at the end

if (save_figure == "yes") {
  print("figures have been saved")
  # save figures to pdf
  dev.off()
}