library(tidyverse)
library(umap)
library(yarrr)
library(ggpubr)
library(multcompView)     # for creating the significance labels
library(ComplexHeatmap)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$concensus_cluster=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$concensus_cluster) , ]
  return(Tukey.labels)
}

file1 <- "C:/Users/niels/Documents/R/leiden_heatmap/clean_data/2020-07-03/TMM_normalised_weighed_clusters_after_998_repeats_of_leiden_filtered.csv"
file2 <- "C:/Users/niels/Documents/R/leiden_heatmap/clean_data/2020-07-03/TMM_normalised_weighed_clusters_after_999_repeats_of_leiden_filtered.csv"
file3 <- "C:/Users/niels/Documents/R/leiden_heatmap/clean_data/2020-07-03/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered.csv"
file4 <- "C:/Users/niels/Documents/R/leiden_heatmap/clean_data/2020-07-03/TMM_normalised_weighed_clusters_after_1001_repeats_of_leiden_filtered.csv"
file5 <- "C:/Users/niels/Documents/R/leiden_heatmap/clean_data/2020-07-03/TMM_normalised_weighed_clusters_after_1002_repeats_of_leiden_filtered.csv"

abs_no_fold_change <- "R/01_data_in_use/non_significant_abs_kruskal_test.csv"

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/R/leiden_heatmap/figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

dat <- read.csv(file1)

dat$cluster_1 <- dat$cluster
dat$X <- NULL
dat$cluster <- NULL


new_dat <- read.csv(file2)
dat$cluster_2 <- new_dat$cluster

new_dat <- read.csv(file3)
dat$cluster_3 <- new_dat$cluster

new_dat <- read.csv(file4)

new_dat$cluster<- new_dat$cluster %>% 
  recode('1' = 1L, 
         '2' = 2L, 
         '3' = 3L, 
         '4' = 4L, 
         '5' = 5L, 
         '6' = 6L,
         '7' = 7L,
         '8' = 8L,
         '9' = 9L,
         '10' = 10L)

dat$cluster_4 <- new_dat$cluster

new_dat <- read.csv(file5)
dat$cluster_5 <- new_dat$cluster

ambigious <- vector()

true_cluster <- vector()

for (i in c(1:length(dat$plate_number))) {
  value <- dat[i,c((length(dat) - 4):length(dat))] %>% 
    as.vector() %>% 
    as.numeric
  
  true_cluster <- append(true_cluster, names(sort(table(value), 
                                                         decreasing = TRUE)[1]))
  
  value <- sd(value)

  if (value != 0) {
    ambigious <- append(ambigious, "yes")
  }
  if (value == 0) {
    ambigious <- append(ambigious, "no")}

}

dat$concensus_cluster <- as.numeric(true_cluster)
dat$ambigious <- ambigious



#place sample_id as row names
TMM_normalised <- dat %>%
  remove_rownames %>%
  column_to_rownames("sample_id")

#load dataframe that created list with shared non chaninging ABs

non_significant_abs <- read_csv(file = abs_no_fold_change)

# create factor out of the antibody column
non_significant_abs <- as.factor(non_significant_abs$antibody)


# create new data frame to only use the non significant the significant abs
TMM_normalised <- TMM_normalised %>%
  dplyr::select(all_of(non_significant_abs))

# make a matrix where you can compare each cell with all other cells

TMM_normalised <- as.matrix(TMM_normalised)


## Do UMAP and add to the visualisation data set
# set seed for umap (seed for leiden is not via R)
set.seed(1)

# # Run umap
umap_out <- umap(TMM_normalised)


umap_1 <- umap_out$layout[,1]
umap_2 <- umap_out$layout[,2]

dat$umap_1 <- umap_1
dat$umap_2 <- umap_2


TMM_final_cluster <- dat %>% 
  pivot_longer(names_to = "cluster_series", values_to = "cluster", 73:78)

pdf(file=paste(figure_dir, "/final_comparisons.pdf", sep = ""),
      width=10, height=10)


#FOR clusters on umap

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=ambigious), size = 4) + 
  ggtitle("") +
  labs(subtitle="", col = "Ambigious") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)


TMM_final_cluster <- TMM_final_cluster %>% 
  filter(cluster_series == "concensus_cluster")

TMM_final_cluster$cluster_series <- NULL
TMM_final_cluster$ambigious <- NULL
TMM_final_cluster$umap_1 <- NULL
TMM_final_cluster$umap_2 <- NULL

TMM_final_cluster <- TMM_final_cluster %>% 
  dplyr::select("sample_id", "plate_number", "treatment", "cluster", everything())

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(concensus_cluster)), size = 4) + 
  ggtitle("") +
  labs(subtitle="", col = "Cluster") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)

dev.off()


# without legend
pdf(file=paste(figure_dir, "/final_comparisons_without_legend.pdf", sep = ""),
    width=10, height=10)


#FOR clusters on umap

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=ambigious), size = 4, show.legend = FALSE) + 
  ggtitle("") +
  labs(subtitle="", col = "Ambigious") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(concensus_cluster)), size = 4, show.legend =  FALSE) + 
  ggtitle("") +
  labs(subtitle="", col = "Cluster") + 
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_colour_manual(values = palette_plots)

dev.off()

total_abs <- factor(colnames(TMM_final_cluster[5:length(TMM_final_cluster)]))

dat <- dat %>% 
  pivot_longer(names_to = "antibody", values_to = "intensity", 4:72)

pdf(file=paste(figure_dir, "/figures_feature_plots.pdf", sep = ""),
    width=5, height=5)

#create scatter plot of umap
for (ab in sort(unique(dat$antibody))) {
  temp_dat <- dat %>% 
    filter(antibody == ab)
  
  temp_dat$concensus_cluster <- as_factor(temp_dat$concensus_cluster)
  
  feature_plot <- ggplot(temp_dat) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity), size = 2) + 
    ggtitle(ab) +
    labs(subtitle="", col = "Z-score \nnormalised \nintensity") + 
    theme(plot.title = element_text(size=10, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))) - 0.05,
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))),
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))) + 0.05,
                                                     1))
  plot(feature_plot)
  
  # What is the effect of the treatment on the value ?
  model=lm(temp_dat$intensity ~ temp_dat$concensus_cluster)
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'temp_dat$concensus_cluster', conf.level=0.95)
  
  # Create the labels
  LABELS <- generate_label_df(TUKEY , "temp_dat$concensus_cluster")
  
  
  boxplot <- ggboxplot(temp_dat, x = "concensus_cluster", y = "intensity", color = "concensus_cluster", 
                       add = "jitter") +
    rotate_x_text(angle = 45)+
    theme_light() +
    ylab("Z-score normalised antibody intensity") +
    xlab("Cluster") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_blank(),
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

  # mirrors the y-axis when pos y is bigger
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) >= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(-(max(temp_dat$intensity) + 
                                                      0.1 * max(temp_dat$intensity)), 
                                                  max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity)), label = LABELS$Letters)
  } 
  # mirrors the y-axis when neg y is bigger 
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) <= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(min(temp_dat$intensity), 
                                                  abs(min(temp_dat$intensity)))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = abs(min(temp_dat$intensity))), label = LABELS$Letters)
    
  }
  
  plot(boxplot)
  
  temp_dat <- temp_dat %>% 
    filter(treatment == "EGF")
  
  # What is the effect of the treatment on the value ?
  model=lm(temp_dat$intensity ~ temp_dat$concensus_cluster)
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'temp_dat$concensus_cluster', conf.level=0.95)
  
  # Create the labels
  LABELS <- generate_label_df(TUKEY , "temp_dat$concensus_cluster")
  
  
  boxplot <- ggboxplot(temp_dat, x = "concensus_cluster", y = "intensity", color = "concensus_cluster", 
                       add = "jitter") +
    rotate_x_text(angle = 45)+
    theme_light() +
    ylab("Z-score normalised antibody intensity") +
    xlab("Cluster") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_blank(),
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
  
  # mirrors the y-axis when pos y is bigger
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) >= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(-(max(temp_dat$intensity) + 
                                                      0.1 * max(temp_dat$intensity)), 
                                                  max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity)), label = LABELS$Letters)
  } 
  # mirrors the y-axis when neg y is bigger 
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) <= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(min(temp_dat$intensity), 
                                                  abs(min(temp_dat$intensity)))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = abs(min(temp_dat$intensity))), label = LABELS$Letters)
    
  }
  
  plot(boxplot)

}

dev.off()

# without legends

pdf(file=paste(figure_dir, "/figures_feature_plots_without_legend.pdf", sep = ""),
    width=5, height=5)

#create scatter plot of umap
for (ab in sort(unique(dat$antibody))) {
  temp_dat <- dat %>% 
    filter(antibody == ab)
  
  temp_dat$concensus_cluster <- as_factor(temp_dat$concensus_cluster)
  
  feature_plot <- ggplot(temp_dat) + 
    geom_point(aes(x=umap_1, y=umap_2, color=intensity), size = 2, show.legend = FALSE) + 
    ggtitle(ab) +
    labs(subtitle="", col = "Z-score \nnormalised \nintensity") + 
    theme(plot.title = element_text(size=10, hjust = 0.5),
          plot.subtitle = element_text(size=10),
          plot.caption = element_text(size = 10),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    scale_colour_gradientn(colours = c("#0000CC", "white", "#FF0000"), 
                           space = "Lab", values = c(0,
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))) - 0.05,
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))),
                                                     1 - (max(temp_dat$intensity) * 1 / (max(temp_dat$intensity) - min(temp_dat$intensity))) + 0.05,
                                                     1))
  plot(feature_plot)
  
  # What is the effect of the treatment on the value ?
  model=lm(temp_dat$intensity ~ temp_dat$concensus_cluster)
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'temp_dat$concensus_cluster', conf.level=0.95)
  
  # Create the labels
  LABELS <- generate_label_df(TUKEY , "temp_dat$concensus_cluster")
  

  boxplot <- ggboxplot(temp_dat, x = "concensus_cluster", y = "intensity", color = "concensus_cluster", 
                       add = "jitter") +
    rotate_x_text(angle = 45)+
    theme_light() +
    ylab("Z-score normalised antibody intensity") +
    xlab("Cluster") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_blank(),
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
    scale_colour_manual(values = palette_plots) +
    rremove("legend")
  
  # mirrors the y-axis when pos y is bigger
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) >= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(-(max(temp_dat$intensity) + 
                                                    0.1 * max(temp_dat$intensity)), 
                                                max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity)), label = LABELS$Letters)
  } 
  # mirrors the y-axis when neg y is bigger 
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) <= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(min(temp_dat$intensity), 
                                                abs(min(temp_dat$intensity)))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = abs(min(temp_dat$intensity))), label = LABELS$Letters)
    
  }
  
  plot(boxplot)
  
  temp_dat <- temp_dat %>% 
    filter(treatment == "EGF")
  
  # What is the effect of the treatment on the value ?
  model=lm(temp_dat$intensity ~ temp_dat$concensus_cluster)
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'temp_dat$concensus_cluster', conf.level=0.95)
  
  # Create the labels
  LABELS <- generate_label_df(TUKEY , "temp_dat$concensus_cluster")
  
  
  boxplot <- ggboxplot(temp_dat, x = "concensus_cluster", y = "intensity", color = "concensus_cluster", 
                       add = "jitter") +
    rotate_x_text(angle = 45)+
    theme_light() +
    ylab("Z-score normalised antibody intensity") +
    xlab("Cluster") +
    facet_wrap(facets = vars(antibody), ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_light() +
    theme(plot.title = element_blank(),
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
    scale_colour_manual(values = palette_plots) +
    rremove("legend")
  
  # mirrors the y-axis when pos y is bigger
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) >= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(-(max(temp_dat$intensity) + 
                                                      0.1 * max(temp_dat$intensity)), 
                                                  max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity)), label = LABELS$Letters)
  } 
  # mirrors the y-axis when neg y is bigger 
  if (max(temp_dat$intensity) + 0.1 * max(temp_dat$intensity) <= abs(min(temp_dat$intensity))) {
    boxplot <- boxplot + coord_cartesian(ylim = c(min(temp_dat$intensity), 
                                                  abs(min(temp_dat$intensity)))) +
      geom_text(data = LABELS, aes(x = concensus_cluster, y = abs(min(temp_dat$intensity))), label = LABELS$Letters)
    
  }
  
  plot(boxplot)
  
}

dev.off()


# create new data frame to only use the non significant the significant abs
TMM_normalised_not_significant <- TMM_final_cluster %>%
  select(sample_id, plate_number, treatment, cluster, all_of(non_significant_abs))

TMM_normalised_significant <- TMM_final_cluster %>% 
  select(sample_id, plate_number, treatment, !all_of(non_significant_abs))

TMM_heatmap <- TMM_normalised_not_significant %>% 
  select(!one_of(c("cluster","plate_number", "treatment", "sample_id")))

rownames(TMM_heatmap) <- TMM_final_cluster$sample_id

TMM_heatmap <- data.matrix(TMM_heatmap)


pdf(file=paste(figure_dir, "/heatmap_antibody_intensities_cellular_level.pdf", sep = ""),
    width=11, height=8, paper = "a4r")

# create heatmap

# set the colour palette for the different plots
palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

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
  select(!one_of(c("cluster","plate_number", "treatment", "sample_id")))

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

dev.off()

