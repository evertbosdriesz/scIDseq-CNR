## LIBRARIES
library(tidyverse)

## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

input_file_clusters <- "R/01_data_in_use/TMM_normalised_z_transformed_concencus_clusters.csv" 

input_file_TMM_normalisation <- "R/01_data_in_use/TMM_normalised_matrix.tsv"

treatment_1 <- "ip70S6K_EGF"
treatment_2 <- "EGF"
treatment_3 <- "iRSK_EGF"

save_table = "yes"

output_file <- "TMM_normalised_with_clusters_matrix.csv"

## CREATE OUTPUT FOLDERS IF NEED BE

data_dir <- paste(paste(getwd(),'/R/normalisations/normalised_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}
