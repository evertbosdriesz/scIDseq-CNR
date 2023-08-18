## code

file1 <- "./clean_data/next_step/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered_1.csv"
file2 <- "./clean_data/next_step/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered_2.csv"
file3 <- "./clean_data/next_step/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered_3.csv"
file4 <- "./clean_data/next_step/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered_4.csv"
file5 <- "./clean_data/next_step/TMM_normalised_weighed_clusters_after_1000_repeats_of_leiden_filtered_5.csv"

abs_no_fold_change <- "../../../data/annotations/non_significant_abs_kruskal_test.csv"

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


# create new data frame to only use the non significant abs (cell state markers)
TMM_normalised <- TMM_normalised %>%
  select(all_of(non_significant_abs))

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

pdf(file=paste(figure_dir, "/comparisons_between_leiden_iterations.pdf", sep = ""),
      width=8, height=16/3)

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

#create scatter plot of umap
ggplot(TMM_final_cluster) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(cluster)), size = 2, show.legend = FALSE) + 
  facet_wrap(facets = "cluster_series") +
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

if (save_figure == "yes") {
  dev.off()
}

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_ambiguity.pdf", sep = ""),
    width=3, height=3)
}
  
#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=ambigious), size = 2) + 
  labs(subtitle="Clustered on cell state marker antibodies", col = "Ambigious") + 
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
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=ambigious), size = 2, show.legend = FALSE) + 
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
        strip.text = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}

TMM_final_cluster <- TMM_final_cluster %>% 
  filter(cluster_series == "concensus_cluster")

TMM_final_cluster$cluster_series <- NULL
TMM_final_cluster$ambigious <- NULL

# order correctly
TMM_final_cluster <- TMM_final_cluster %>% 
  select("sample_id", "plate_number", "treatment", "cluster", "umap_1", "umap_2", everything())

# save information with umaps
write_tsv(TMM_final_cluster, paste(data_dir, "/TMM_normalised_z_transformed_clusters_with_umap.tsv", sep = ""))


TMM_final_cluster$umap_1 <- NULL
TMM_final_cluster$umap_2 <- NULL

TMM_final_cluster <- TMM_final_cluster %>% 
  select("sample_id", "plate_number", "treatment", "cluster", everything())


if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/umap_concensus_clusters.pdf", sep = ""),
      width=3, height=3)
}

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(concensus_cluster)), size = 2) + 
  labs(col = "Cluster") + 
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

ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(concensus_cluster)), size = 2, show.legend = FALSE) + 
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
        strip.text = element_blank(),
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
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(treatment)), size = 2) + 
  labs(col = "Cluster") + 
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

ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(treatment)), size = 2, show.legend = FALSE) + 
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
        strip.text = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = palette_plots)

if (save_figure == "yes") {
  dev.off()
}

df_check <- TMM_final_cluster %>% 
  group_by(cluster, treatment) %>% 
  summarise(count = n())

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/distribution_plots.pdf", sep = ""),
      width=5, height=3)
}

ggplot(data = df_check, aes(x = factor(cluster), y = count, fill = treatment)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Fraction of cells") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = palette_plots)


ggplot(TMM_final_cluster, aes(x = factor(cluster), fill = treatment)) +
  geom_bar() +
  scale_fill_manual(values = palette_plots) +
  ylab("Cell count") +
  xlab("Cluster") +
  labs(fill = "Treatment") +
  theme(plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        plot.caption = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, colour = "black"),
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

write.csv(TMM_final_cluster, paste(data_dir, "/TMM_normalised_concencus_clusters.csv", sep = ""))