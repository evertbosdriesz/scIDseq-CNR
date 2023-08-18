## CODE

# load dataframe

TMM_normalised <- read_csv(file = input_file, col_names = TRUE)
rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

## select treatments
TMM_normalised <- TMM_normalised %>% 
  filter(treatment == "EGF" | treatment == "ip70S6K_EGF" | treatment == "iRSK_EGF")


# transform plate_number column into class factor in order for the 
# plate_numbers to be ordered correctly in the figure

TMM_normalised$plate_number <- factor(TMM_normalised$plate_number, 
                                   level = c("plate_2", "plate_3", 
                                             "plate_4", "plate_6", 
                                             "plate_7", "plate_8", 
                                             "plate_9", "plate_10", 
                                             "plate_11", "plate_12", 
                                             "plate_13"))

# transform dataframe to long format in order to facet, place the antibodies
# in the 'ab' column and the intensities for each ab in the
# 'intensity' column

TMM_normalised_long <- gather(TMM_normalised, key = "antibody", 
                              value = "intensity", -(1:3))

# alphabetise colum with antibodys in order for the facetting to take place in
# alphabetical order

TMM_normalised_long <- TMM_normalised_long %>%
  arrange(antibody)

# create the graph wih geom_density_ridges
# (https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html)

if (save_figures == "yes") {
  pdf(file = paste(figure_dir, "/kernel_density_estimates_treatment.pdf", sep = ""), paper = "a4",
      width = 11.7, height = 16.5)
}
  
# plot in kernel density estimates 


# creates the right colour palette
plot <- ggplot(data = TMM_normalised_long, mapping = aes(x = intensity, y = treatment)) +
  geom_density_ridges(aes(fill = treatment), 
                      position = "identity", 
                      scale = 2, #ensures that the graphs overlap
                      alpha=0.5,
                      size = 0.1,
                      show.legend = TRUE) +
  theme_light() +
  labs(title = "Kernel Density Estimates", fill = "Treatment") +
  ylab("Density") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette_plots) +
  theme_light() +
  theme(plot.title = element_text(hjust=0.5, size = 5)) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.5, "cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(facets = vars(antibody), ncol = 7) +
  theme(strip.text.x = element_text(size = 5, color = "black")) +
  theme(strip.background = element_blank()) +
  theme(legend.background = element_rect()) +
  coord_cartesian(xlim =c(-3, 3), ylim = c(1, 5)) + # ensures graphs are plotted correctly
  theme(plot.title = element_text(hjust = 0.5, size = 10))

print(plot)

if(save_figures == "yes"){
  dev.off()
}  

if (save_figures == "yes") {
  pdf(file = paste(figure_dir, "/kernel_density_estimates_plate_number.pdf", sep = ""), paper = "a4",
      width = 11.7, height = 16.5)
}

plot <- ggplot(data = TMM_normalised_long, mapping = aes(x = intensity, y = plate_number)) +
  geom_density_ridges(aes(fill = plate_number), 
                      position = "identity", 
                      alpha=0.5, 
                      scale = 3, #ensures that the graphs overlap
                      size = 0.1,
                      show.legend = TRUE) +
  theme_light() +
  labs(title = "Kernel Density Estimates", fill = "Plate number") +
  ylab("Density") +
  facet_wrap(facets = vars(antibody), ncol = 7) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_light() +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.5, "cm")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(strip.text.x = element_text(size = 5, color = "black")) +
  theme(strip.background = element_blank()) +
  scale_fill_manual(values = palette_plots) +
  theme(legend.background = element_rect()) +
  coord_cartesian(xlim =c(-3, 3), ylim = c(1, 14)) + # ensures graphs are plotted correctly
  theme(plot.title = element_text(hjust = 0.5, size = 10))

print(plot)

if(save_figures == "yes"){
  dev.off()
}  
