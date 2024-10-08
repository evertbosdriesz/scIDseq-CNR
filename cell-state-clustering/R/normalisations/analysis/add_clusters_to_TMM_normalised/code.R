## CODE

# load data frames
TMM_normalised <- read_tsv(file = input_file_TMM_normalisation)
TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment_1| treatment == treatment_2 | treatment == treatment_3)

TMM_clusters <- read_csv(file = input_file_clusters)
rownames(TMM_clusters) <- TMM_clusters$X1
TMM_clusters$X1 <- NULL

#check if sample_ids are still in the same order
if(unique(TMM_normalised$sample_id == TMM_clusters$sample_id) == TRUE) {
  print("rows are in the same order, proceed") } else {
    print("rows are not in the same order, ensure they are in the same order before proceeding") }

# place cluster in TMM_normalised
TMM_normalised$cluster <- TMM_clusters$cluster

TMM_normalised <- TMM_normalised %>% 
  select(sample_id, plate_number, treatment, cluster, everything())

# save new data frame in the normalised data_folder
if(save_table == "yes") {
  write.csv(TMM_normalised, file = paste(data_dir, "/", output_file, sep = ""))
  print("Table has been saved in your normalised data folder")
} else {
  print("Table has not been saved")
}