library(tidyverse)
library(umap)

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

pdf(file=paste(figure_dir, "/final_comparisons.pdf", sep = ""),
      width=10, height=10)


#FOR clusters on umap

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

#create scatter plot of umap
ggplot(TMM_final_cluster) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(cluster)), size = 2) + 
  facet_wrap(facets = "cluster_series") +
  ggtitle("UMAP on Cell State Markers") +
  labs(subtitle="Clustered on cell state marker antibodies",
       col = "Cluster") + 
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
  geom_point(aes(x=umap_1, y=umap_2, color=ambigious), size = 4) + 
  ggtitle("") +
  labs(subtitle="", col = "Ambigious") + 
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


TMM_final_cluster <- TMM_final_cluster %>% 
  filter(cluster_series == "concensus_cluster")

TMM_final_cluster$cluster_series <- NULL
TMM_final_cluster$ambigious <- NULL
TMM_final_cluster$umap_1 <- NULL
TMM_final_cluster$umap_2 <- NULL

TMM_final_cluster <- TMM_final_cluster %>% 
  select("sample_id", "plate_number", "treatment", "cluster", everything())

#create scatter plot of umap
ggplot(dat) + 
  geom_point(aes(x=umap_1, y=umap_2, color=as.factor(concensus_cluster)), size = 4) + 
  ggtitle("") +
  labs(subtitle="", col = "Cluster") + 
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

write.csv(TMM_final_cluster, paste(data_dir, "/TMM_normalised_concencus_clusters.csv", sep = "")) 

dev.off()
