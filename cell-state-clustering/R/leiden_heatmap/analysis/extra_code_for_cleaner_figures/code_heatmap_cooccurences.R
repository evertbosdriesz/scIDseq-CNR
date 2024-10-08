cor_clusters <- read_csv("R/leiden_heatmap/clean_data/2020-07-15/cooccurences_after_1000_repeats_of_leiden_filtered_1.csv")

cor_clusters$X1 <- NULL

cor_clusters <- as.matrix(cor_clusters)

if (save_figure == "yes") {
  pdf(file=paste(figure_dir, "/heatmap_cooccurences_after_", no_repeats, "_repeats_of_leiden_filtered.pdf", sep = ""),
      width=5, height=4)
}

ht2 <- Heatmap(cor_clusters, 
        name = "Co-occurence distance",
        column_title = "",
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        heatmap_legend_param = list(legend_direction = "horizontal")
)

draw(ht2, heatmap_legend_side = "bottom")

dev.off()
