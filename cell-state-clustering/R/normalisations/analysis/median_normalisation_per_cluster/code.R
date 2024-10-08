## CODE

# load comma separated file into dataframe

TMM_normalised <- read_csv(file = input_file)
rownames(TMM_normalised) <- TMM_normalised$X1
TMM_normalised$X1 <- NULL

TMM_normalised$cluster <- as_factor(TMM_normalised$cluster)

TMM_median_normalised <- tibble()

for (cluster_no in c(1:length(unique(TMM_normalised$cluster)))) {
  cluster_df <- TMM_normalised %>% 
    filter(cluster == cluster_no)
  cluster_df_long <- cluster_df %>% 
    gather(key = 'antibody', value = 'intensity', -sample_id, -plate_number, -treatment, -cluster)
  for (ab in sort(colnames(TMM_normalised[5:length(TMM_normalised)]))) {
    a <- cluster_df_long %>% 
      filter(cluster == cluster_no) %>% 
      filter(treatment == "EGF") %>% 
      filter(antibody == ab) %>% 
      summarise(median = median(intensity))
    a <- as.numeric(a)
    cluster_df[ab] <- cluster_df[ab] / a
  }
  TMM_median_normalised <- rbind(TMM_median_normalised, cluster_df)
}  

TMM_median_normalised_long <- TMM_median_normalised %>% 
  gather(key = 'antibody', value = 'intensity', -sample_id, -plate_number, -treatment, -cluster)

# save new data frame in the normalised data_folder
if(save_table == "yes") {
  write.csv(TMM_median_normalised, file = paste(data_dir, "/", output_file, sep = ""))
  print("Table has been saved in your normalised data folder")
} else {
  print("Table has not been saved")
}
