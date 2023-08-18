## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)

rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

# transform TMM_normalised into long format
TMM_normalised <- TMM_normalised %>% 
  pivot_longer(names_to = "antibody", values_to = "intensity", 4:length(TMM_normalised))

TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment_2 | treatment == treatment_1)



if(save_figure == "yes"){
  pdf(file = paste(figure_dir, "/violin_plots_", treatment_2, "_vs_", treatment_1, ".pdf", sep = ""),
      width = 2, height = 4)
}

# create violin plots of the antibody intensities
for (ab in sort(unique(TMM_normalised$antibody))) {
  TMM_temp <- TMM_normalised %>% 
    filter(antibody == ab)
  
  ymax1 <- max(TMM_temp$intensity)
  
  TMM_temp$treatment <- as.factor(TMM_temp$treatment)
  
  #create separate data frames for each treatment
  table_treatment_1 <- TMM_temp %>% 
    filter(treatment == treatment_1) %>% 
    pull(intensity)
  
  table_treatment_2 <- TMM_temp %>% 
    filter(treatment == treatment_2) %>% 
    pull(intensity)
  
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
  
  # create palette
  palette_right <- c(palette_plots[1], palette_plots[correct_colour])
  
  plot <- ggviolin(TMM_temp, x = "treatment", y = "intensity",
                    color = "treatment", palette = palette_right,
                    add = "jitter", alpha = 0.2) +
    geom_text(aes(x = 1.5, y = ymax1),
              label = sig,
              size = 10) +
    ggtitle(ab) +
    theme(legend.position = "bottom",
          axis.title.y.left = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    rremove("legend")
  
  print(plot)
}

if(save_figure == "yes") {
  dev.off()
} 

if(save_figure == "yes"){
  pdf(file = paste(figure_dir, "/violin_plots_", treatment_2, "_vs_", treatment_1, "_with_legend_and_axes.pdf", sep = ""),
      width = 2, height = 4)
}

plot <- ggviolin(TMM_temp, x = "treatment", y = "intensity",
                  color = "treatment", palette = palette_right,
                  add = "jitter", alpha = 0.2) +
  geom_text(aes(x = 1.5, y = ymax1),
            label = sig,
            size = 10) +
  ggtitle(ab) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot

if(save_figure == "yes") {
  dev.off()
} 